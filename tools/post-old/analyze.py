#!/usr/bin/python

import sys
import os
import re
import string
import copy
import math
from operator import sub
from sympy import Symbol, ask, Q
from sympy.parsing.sympy_parser import parse_expr
import pprint

from params import doSymSubs, isSpeciesArray
from common import options, numIters, diff, prod, int_types, float_types

import parser
import grapher
import spreadsheet

pp = pprint.PrettyPrinter(width=126)

def rangeSize(x):
  return x[1] - x[0] + 1

def getMember(x, tag):
  for item in x:
    if item['name'] == tag:
      return item

def getR(arrays):
  return filter(lambda x: x['accesstype'] == 'readonly', arrays)
def getRW(arrays):
  return filter(lambda x: x['accesstype'] == 'readwrite', arrays)
def getW(arrays):
  return filter(lambda x: x['accesstype'] == 'writeonly' or \
                          x['accesstype'] == 'writeread', arrays)

def getRNames(arrays):
  return set(map(lambda x: x['name'], getR(arrays)))
def getRWNames(arrays):
  return set(map(lambda x: x['name'], getRW(arrays)))
def getWNames(arrays):
  return set(map(lambda x: x['name'], getW(arrays)))

def printLoops(info, style='full'):

  t = {'readonly': 'R', 'readwrite': 'RW', 'writeonly': 'W', 'writeread': 'WR'}

  def printLoop(info, indent=0):
    if style == 'simple':
      for linfo in info['loops']:
        print ' '*(indent+1), linfo['linenum'], linfo['loopvar'], linfo['range']
        if 'loops' in linfo:
          printLoop(linfo, indent+2)
    elif style == 'full':
      for linfo in info['loops']:
        print ' '*(indent+1), linfo['linenum'], linfo['sweeps'], linfo['loopvars'], linfo['ranges']
        print ' '*(indent+3), 'FLOPS: adds: %s, muls: %s, divs: %s, spes: %s' % \
                              (linfo['flops']['adds'], linfo['flops']['multiplies'], \
                               linfo['flops']['divides'], linfo['flops']['specials'])
        print ' '*(indent+3), 'ARRAYS:'
        for ainfo in linfo['arrays']:
          print ' '*(indent+5), ainfo['name'], ainfo['copies'], t[ainfo['accesstype']]
          pp.pprint(ainfo['access'])
        print ' '*(indent+3), 'STATE ARRAYS:'
        for ainfo in linfo['stateArrays']:
          print ' '*(indent+5), ainfo['name'], ainfo['copies'], ainfo['access']
        print ' '*(indent+3), 'SCALARS:'
        print ' '*(indent+5), map(lambda x: x['name'], linfo['scalars'])
        if 'loops' in linfo:
          printLoop(linfo, indent+2)
    else:
      pass

  for (fname, finfo) in info.iteritems():
    print fname
    printLoop(finfo)


class StaticAnalysis:

  def __init__(self, filename):
    self.info = parser.XMLParser(filename).info()
    if options.flag_debug:
      pp.pprint(self.info)
      printLoops(self.info, 'simple')
    self.analyzeInfo()
    if options.flag_debug:
      printLoops(self.info, 'full')

  def analyzeInfo(self):

    def getRegisters(scalars, stateArrays, arrays):
      ints = filter(lambda x: x['datatype'] in int_types, scalars)
      floats = filter(lambda x: x['datatype'] in float_types, scalars)
      constInts = filter(lambda x: x['const'] == 'true', ints)
      constFloats = filter(lambda x: x['const'] == 'true', floats)
      SAInts = filter(lambda x: x['datatype'] == 'int', stateArrays)
      SAFloats = filter(lambda x: x['datatype'] == 'double', stateArrays)
      numSAIntRegs = sum(map(lambda x: x['copies'] * len(x['access']), SAInts))
      numSAFloatRegs = sum(map(lambda x: x['copies'] * len(x['access']), SAFloats))
      arrayNames = set(map(lambda x: os.path.splitext(x['name'])[0], arrays))
      return {'ints': len(ints) + numSAIntRegs, 'floats': len(floats) + numSAFloatRegs, \
              'const': {'ints': len(constInts), 'floats': len(constFloats)}, \
              'ptrs': len(arrayNames)}

    # compute access types
    def getAccessTypes(loop):
      def numNonzero(x):
        return sum(map(lambda y: y != 0, x))
      for array in loop['arrays']:
        types = ['cell', 'faces', 'edges', 'all']
        array['stenciltype'] = types[max(map(numNonzero, array['access'].keys()))]

    # compute number of cells, pencils, and planes accessed by stencil pattern
    def getPencilsAndPlanes(loop):

      # raw number of pencils and planes accessed during a corresponding sweep
      def numBWCells(array):
        return array['copies'] * len(set(array['access'].keys()))
      def numBWPencils(array):
        return array['copies'] * len(set(map(lambda x:x[1:3], array['access'].keys())))
      def numBWPlanes(array):
        return array['copies'] * len(set(map(lambda x:x[2], array['access'].keys())))
      def numBWCopies(array):
        return array['copies']

      def maxDiff(x):
        if len(x) == 1:
          return 1
        return max(diff(sorted(x)))

      # return the set of cells in each pencil
      def cellsByPencil(array):
        d = {}
        for access in array['access'].keys():
          if access[1:3] not in d:
            d[access[1:3]] = set()
          d[access[1:3]].add(access[0])
        return d

      # return the set of pencils in each plane
      def pencilsByPlane(array):
        d = {}
        for access in array['access'].keys():
          if access[2] not in d:
            d[access[2]] = set()
          d[access[2]].add(access[1])
        return d

      # maximum gap size between cells/pencils/planes for working set calculation
      def maxCellGap(array):
        return max(map(maxDiff, cellsByPencil(array).values()))
      def maxPencilGap(array):
        return max(map(maxDiff, pencilsByPlane(array).values()))
      def maxPlaneGap(array):
        return maxDiff(map(lambda x:x[2], array['access'].keys()))

      # WS is determined by largest gap between elements and the size of the
      # span of elements: WS = Max - Min + MaxGap.
      def computeWS(pattern, gap = None):
        if options.flag_ignore_gap:
          return len(pattern)
        else:
          pattern = sorted(pattern)
          if not gap:
            gap = maxDiff(pattern)
          return pattern[-1] - pattern[0] + gap

      def numWSPlanes(array, gap = None):
        result = computeWS(set(map(lambda x:x[2], array['access'].keys())), gap)
        return array['copies'] * result

      # Only return non-zero if array has reuse between either planes or pencils
      def numWSReusePlanes(array, gap = None):
        result = numWSPlanes(array, gap) / array['copies']
        if numWSReusePencils(array) == 0 and \
           (gap and result == gap or not gap and result == 1):
          return 0
        return array['copies'] * result

      # WS is determined by largest gap in Y-dim across all Z-planes,
      # then apply same formula as numWSPlanes for each plane and sum.
      def numWSPencils(array, gap = None):
        result = sum(map(lambda x:computeWS(x, gap), pencilsByPlane(array).values()))
        return array['copies'] * result

      # Number of pencils needed for a working set in planes with reuse
      # Differentiate between Z-planes of an array with/without reuse
      def numWSReusePencils(array, gap = None):
        # no reuse if only one pencil in a plane
        pencilsWithReuse = filter(lambda x: len(x) > 1, pencilsByPlane(array).values())
        result = sum(map(lambda x:computeWS(x, gap), pencilsWithReuse))
        return array['copies'] * result

      # TODO: implement these (for now, just use BWCells)
      #       take into account cache line size?
      def numWSCells(array, gap = None):
        return numBWCells(array)
      def numWSReuseCells(array, gap = None):
        return numBWCells(array)

      def getGhostZone(array):
        def minmax(x, y):
          return [min(x[0], y), max(x[1], y)]
        ghost = [[0,0]] * options.dim
        for access in array['access'].keys():
          ghost = map(minmax, ghost, access)
        return ghost

      # get lists of arrays by access type
      readArrays = getR(loop['arrays'])
      readWriteArrays = getRW(loop['arrays'])
      writeArrays = getW(loop['arrays'])

      # compute gap sizes needed to be covered by working sets
      cellGap = 1
      pencilGap = 1
      planeGap = 1
      if readArrays:
        cellGap = max(map(maxCellGap, readArrays))
        pencilGap = max(map(maxPencilGap, readArrays))
        planeGap = max(map(maxPlaneGap, readArrays))
      loop['cellGap'] = cellGap
      loop['pencilGap'] = pencilGap
      loop['planeGap'] = planeGap

      # compute ghost zone, bandwidth, and working set info per array
      for array in loop['arrays']:
        array['ghost'] = getGhostZone(array)
        array['BW'] = {}
        array['BW']['numCells'] = numBWCells(array)
        array['BW']['numPencils'] = numBWPencils(array)
        array['BW']['numPlanes'] = numBWPlanes(array)
        array['BW']['numCopies'] = numBWCopies(array)
        array['WS'] = {}
        array['WS']['numCells'] = numWSCells(array, cellGap)
        array['WS']['numReuseCells'] = numWSReuseCells(array, cellGap)
        array['WS']['numPencils'] = numWSPencils(array, pencilGap)
        array['WS']['numReusePencils'] = numWSReusePencils(array, pencilGap)
        array['WS']['numPlanes'] = numWSPlanes(array, planeGap)
        array['WS']['numReusePlanes'] = numWSReusePlanes(array, planeGap)

      # compute aggregate bandwidth and working set info by access type
      loop['BW'] = {'numCells': {}, 'numPencils': {}, 'numPlanes': {}, 'numArrays': {}}
      loop['WS'] = {'numCells': {}, 'numPencils': {}, 'numPlanes': {}, 'numArrays': {}}

      loop['BW']['numCells']['R'] = sum(map(numBWCells, readArrays))
      loop['BW']['numCells']['RW'] = sum(map(numBWCells, readWriteArrays))
      loop['BW']['numCells']['W'] = sum(map(numBWCells, writeArrays))

      loop['BW']['numPencils']['R'] = sum(map(numBWPencils, readArrays))
      loop['BW']['numPencils']['RW'] = sum(map(numBWPencils, readWriteArrays))
      loop['BW']['numPencils']['W'] = sum(map(numBWPencils, writeArrays))

      loop['BW']['numPlanes']['R'] = sum(map(numBWPlanes, readArrays))
      loop['BW']['numPlanes']['RW'] = sum(map(numBWPlanes, readWriteArrays))
      loop['BW']['numPlanes']['W'] = sum(map(numBWPlanes, writeArrays))

      loop['BW']['numArrays']['R'] = sum(map(lambda x: x['copies'], readArrays))
      loop['BW']['numArrays']['RW'] = sum(map(lambda x: x['copies'], readWriteArrays))
      loop['BW']['numArrays']['W'] = sum(map(lambda x: x['copies'], writeArrays))

      loop['WS']['numCells']['RR'] = sum(map(lambda x:numWSReuseCells(x, cellGap), readArrays))
      loop['WS']['numCells']['R'] = sum(map(lambda x:numWSCells(x, cellGap), readArrays))
      loop['WS']['numCells']['RW'] = sum(map(lambda x:numWSCells(x, cellGap), readWriteArrays))
      loop['WS']['numCells']['W'] = sum(map(lambda x:numWSCells(x, cellGap), writeArrays))

      loop['WS']['numPencils']['RR'] = sum(map(lambda x:numWSReusePencils(x, pencilGap), readArrays))
      loop['WS']['numPencils']['R'] = sum(map(lambda x:numWSPencils(x, pencilGap), readArrays))
      loop['WS']['numPencils']['RW'] = sum(map(lambda x:numWSPencils(x, pencilGap), readWriteArrays))
      loop['WS']['numPencils']['W'] = sum(map(lambda x:numWSPencils(x, pencilGap), writeArrays))

      loop['WS']['numPlanes']['RR'] = sum(map(lambda x:numWSReusePlanes(x, planeGap), readArrays))
      loop['WS']['numPlanes']['R'] = sum(map(lambda x:numWSPlanes(x, planeGap), readArrays))
      loop['WS']['numPlanes']['RW'] = sum(map(lambda x:numWSPlanes(x, planeGap), readWriteArrays))
      loop['WS']['numPlanes']['W'] = sum(map(lambda x:numWSPlanes(x, planeGap), writeArrays))

    # set memTableFlag to true if we want input arrays to be present from
    # beginning of the function and output arrays until the end
    # set memTableFlag to false if arrays can be pulled into memory
    # and written out on an as needed basis (e.g. like with a cache)
    def getLiveArrays(fname, finfo, memTableFlag = True):

      # first, some helper functions to compute and manipulate access ranges

      def computeAccessRange(ranges, ghost):
        return map(lambda x, y: map(lambda w, z: doSymSubs(w) + z, x, y), ranges, ghost)

      def mySymMin(x, y):
        # HACK: for now, assume ng equals 4 so we can compare access ranges
        if options.flag_warn and not 'ng==4' in options.warned:
          options.warned.add('ng==4')
          print >> sys.stderr, 'Warning: assuming ng == 4 for PW execution model'
        x = x.subs('ng', 4)
        y = y.subs('ng', 4)
        testExpr = ask(~Q.positive(x - y), Q.positive(Symbol('ng')))
        if type(testExpr) == type(True):
          return x if testExpr else y
        else:
          return parse_expr('min(%s, %s)' % (x, y))

      def mySymMax(x, y):
        # HACK: for now, assume ng equals 4 so we can compare access ranges
        if options.flag_warn and not 'ng==4' in options.warned:
          options.warned.add('ng==4')
          print >> sys.stderr, 'WARNING: assuming ng == 4 for PW execution model'
        x = x.subs('ng', 4)
        y = y.subs('ng', 4)
        testExpr = ask(~Q.negative(x - y), Q.positive(Symbol('ng')))
        if type(testExpr) == type(True):
          return x if testExpr else y
        else:
          return parse_expr('max(%s, %s)' % (x, y))

      def maxAccessRange(x, y):
        # unpack range and stencil type from inputs
        (x, xt) = x
        (y, yt) = y
        maxRange = map(lambda w, z: (mySymMin(w[0], z[0]), mySymMax(w[1], z[1])), x, y)
        if xt == 'all' or yt == 'all':
          maxStencilType = 'all'
        elif xt == 'faces' or yt == 'faces':
          maxStencilType = 'faces'
        elif xt == 'edges' or yt == 'edges':
          maxStencilType = 'edges'
        else:
          maxStencilType = 'cell'
        return (maxRange, maxStencilType)

      # if no custom execution order, use default
      numLoops = len(finfo['loops'])
      if 'execorder' in finfo:
        if len(finfo['execorder']) != numLoops or \
           set(finfo['execorder']) != set(range(0, numLoops)):
          raise Exception('custom execution order must be a permutation of range(0, len(loops))')
      else:
        finfo['execorder'] = range(0, numLoops)

      # gather RW information into table
      table = []
      arrays = {}
      loops_execorder = [finfo['loops'][i] for i in finfo['execorder']]
      for idx, linfo in enumerate(loops_execorder):
        lname = "%s.%03d" % (fname, linfo['linenum'])
        table.append({})
        for ainfo in linfo['arrays']:

          aname = ainfo['name']
          accesstype = ainfo['accesstype']
          copies = ainfo['copies']

          if accesstype == 'readonly':
            table[idx][aname] = ('R', copies)
          elif accesstype == 'writeonly':
            table[idx][aname] = ('W', copies)
          elif accesstype == 'readwrite':
            table[idx][aname] = ('RW', copies)
          elif accesstype == 'writeread':
            table[idx][aname] = ('WR', copies)

          if accesstype == 'readonly':
            dirty = False
          else:
            dirty = True

          if aname not in arrays:
            # compute access range and initialize array analysis data
            accessRange = (computeAccessRange(linfo['ranges'], ainfo['ghost']), ainfo['stenciltype'])
            arrays[aname] = {'firstidx': idx, 'firstaccess': accesstype, \
                             'copies': copies, 'accessrange': accessRange}
          else:
            # update max access range
            accessRange = (computeAccessRange(linfo['ranges'], ainfo['ghost']), ainfo['stenciltype'])
            arrays[aname]['accessrange'] = maxAccessRange(arrays[aname]['accessrange'], accessRange)

            # do some checks
            lastidx = arrays[aname]['lastidx']
            lastaccess = arrays[aname]['lastaccess']
            dirty = dirty or arrays[aname]['dirty']

            if arrays[aname]['copies'] != copies and options.flag_warn:
                print >> sys.stderr, 'Warning: number of copies mismatch: (%s, %s, %s)' % \
                      (aname, arrays[aname]['copies'], copies)

            if (lastaccess == 'writeonly' or lastaccess == 'readwrite') and \
               (accesstype == 'writeonly' or accesstype == 'writeread') and options.flag_warn:
              print >> sys.stderr, 'Warning: value not used between two writes: (%s, %s.%d, %s)' % \
                    (aname, fname, loops_execorder[lastidx]['linenum'], lname)

            # check if array needs to be present in cache since last access
            if accesstype == 'readonly' or accesstype == 'readwrite':
              for i in xrange(lastidx+1, idx):
                table[i][aname] = ('P', copies)

          # update array analysis data
          arrays[aname]['lastidx'] = idx
          arrays[aname]['lastaccess'] = accesstype
          arrays[aname]['dirty'] = dirty

      # compute block working set sizes using computed array max access ranges
      for (idx, linfo) in enumerate(loops_execorder):
        linfo['WS'][ 'numBlocks'] = 0
        linfo['WS']['sizeBlocks'] = 0
        # for each array accessed (R/W/RW) or present (P) during this loop
        for aname in table[idx].keys():
          array = arrays[aname]
          if 'size' not in array:
            box = (Symbol('__BoxX__'), Symbol('__BoxY__'), Symbol('__BoxZ__'))
            acc = map(doSymSubs, map(rangeSize, array['accessrange'][0]))
            # check stencil shape to see if we can trim edges and corners from the box 
            if array['accessrange'][1] == 'faces':
              array['size'] = array['copies'] * (acc[0]*box[1]*box[2] + box[0]*acc[1]*box[2] + \
                                                 box[0]*box[1]*acc[2] - 2*box[0]*box[1]*box[2])
            elif array['accessrange'][1] == 'edges':
              array['size'] = array['copies'] * (prod(acc) - prod(map(sub, acc, box)))
            else:
              # this works for 'all' or 'cell'
              array['size'] = array['copies'] * prod(acc)
            # since this working set is for blocked execution, substitute block sizes for box sizes
            array['size'] = array['size'].subs((('__BoxX__', '__BlockX__'), \
                                                ('__BoxY__', '__BlockY__'), \
                                                ('__BoxZ__', '__BlockZ__')))
          linfo['WS'][ 'numBlocks'] += array['copies']
          linfo['WS']['sizeBlocks'] += array['size']

      # initialize bandwidth numbers to 0
      for linfo in loops_execorder:
        linfo['BW'][ 'numBlocks'] = {'R': 0, 'RW': 0, 'W': 0}
        linfo['BW']['sizeBlocks'] = {'R': 0, 'RW': 0, 'W': 0}

      # for each array, figure out if/when it needs to be streamed to/from memory
      for aname in arrays:
        base = string.split(aname, '.')[0]
        firstidx = arrays[aname]['firstidx']
        lastidx = arrays[aname]['lastidx']

        if base in finfo['nonlocal']:
          streamin = (arrays[aname]['firstaccess'] == 'readonly' or \
                      arrays[aname]['firstaccess'] == 'readwrite')
          streamout = arrays[aname]['dirty']
          if streamin:
            loops_execorder[firstidx]['BW']['numBlocks']['R'] += arrays[aname]['copies']
            loops_execorder[firstidx]['BW']['sizeBlocks']['R'] += arrays[aname]['size']
            if memTableFlag:
              for i in xrange(0, firstidx):
                table[i][aname] = ('P', arrays[aname]['copies'])
          if streamout:
            if streamin:
              loops_execorder[lastidx]['BW']['numBlocks']['RW'] += arrays[aname]['copies']
              loops_execorder[lastidx]['BW']['sizeBlocks']['RW'] += arrays[aname]['size']
              # since we counted a read earlier, don't overcount the extra read
              loops_execorder[lastidx]['BW']['numBlocks']['R'] -= arrays[aname]['copies']
              loops_execorder[lastidx]['BW']['sizeBlocks']['R'] -= arrays[aname]['size']
            else:
              loops_execorder[lastidx]['BW']['numBlocks']['W'] += arrays[aname]['copies']
              loops_execorder[lastidx]['BW']['sizeBlocks']['W'] += arrays[aname]['size']
            if memTableFlag:
              for i in xrange(lastidx + 1, len(table)):
                table[i][aname] = ('P', arrays[aname]['copies'])

      finfo['liveness'] = table

    # translate a hierarchical loop nest into a flat struct that
    #   represents a 3D loop nest
    # NOTE: This function relies on some syntax specific to the CNS/SMC
    #   codes.  Currently only handles certain types of loop nests
    #   with the number of chemical species being looped over in the
    #   outermost or innermost loop of a quadruply nested loop, and the
    #   X,Y,Z spatial loops being perfectly nested within one another.
    #   The loop over the number of chemical species does not need to
    #   be perfectly nested with respect to the X,Y,Z loops.
    def flattenLoopNest(linfo):

      def getRanges(linfo):
        result = [[linfo['loopvar']], [linfo['range']], [linfo['stride']]]
        if len(linfo['loops']) > 0:
          temp = getRanges(linfo['loops'][0])
          result[0].extend(temp[0])
          result[1].extend(temp[1])
          result[2].extend(temp[2])
        return result

      def isSpatial(x):
        return re.search('d?lo\(\d\)', x[0]) and \
               re.search('d?hi\(\d\)', x[1])

      # triply nested spatial loop
      def isType0(ranges):
        return len(ranges) == 3 and all(map(isSpatial, ranges))

      # quadruply nested loop with spatial loops 2, 3, and 4
      def isType1(ranges):
        return len(ranges) == 4 and all(map(isSpatial, ranges[1:]))

      # quadruply nested loop with spatial loops 1, 2, and 3
      def isType2(ranges):
        return len(ranges) == 4 and all(map(isSpatial, ranges[:-1]))

      # merge info from child loop linfo2 into parent loop linfo1
      def mergeLoops(linfo1, linfo2):

        def mergeAccess(array1, array2):
          access1 = array1['access']
          access2 = array2['access']
          # HACK: if there's a numCopies mismatch, try to recognize if there's an implicit
          # knowledge of the number of elements in the array
          if array1['copies'] == 1 and array2['copies'] != 1:
            maxIndex = max(map(lambda x: x[0], access1.keys()))
            print >> sys.stderr, 'Warning: During merge of %s: assuming %s equals %s' % \
                                 (array2['name'], array2['copies'], maxIndex)
            if access2.keys() != [(0,)]:
              raise Exception('unsupported access merge')
            for idx in xrange(1, maxIndex + 1):
              if (idx,) in access1:
                access1[(idx,)]['reads' ] += access2[(0,)]['reads' ]
                access1[(idx,)]['writes'] += access2[(0,)]['writes']
              else:
                access1[(idx,)] = access2[(0,)]
          else:
            for index in access2:
              if index in access1:
                access1[index]['reads' ] += access2[index]['reads' ]
                access1[index]['writes'] += access2[index]['writes']
              else:
                access1[index] = access2[index]
          return access1

        iters = numIters(linfo2['range'])
        if linfo2['loops']:
          raise Exception('child loop cannot have nested loops of its own')

        # process flops
        for floptype in ['adds', 'multiplies', 'divides', 'specials']:
          linfo1['flops'][floptype] = linfo1['flops'][floptype] + iters * linfo2['flops'][floptype]

        # process scalars
        for scalar2 in linfo2['scalars']:
          scalar2['reads'] *= iters
          scalar2['writes'] *= iters
          scalar1 = getMember(linfo1['scalars'], scalar2['name'])
          if scalar1:
            scalar1['reads'] += scalar2['reads']
            scalar1['writes'] += scalar2['writes']
          else:
            linfo1['scalars'].append(scalar2)

        # process spatial arrays and state arrays
        for category in ['arrays', 'stateArrays']:
          for array2 in linfo2[category]:
            if category == 'arrays' and isSpeciesArray(array2['name']) or \
               category == 'stateArrays' and array2['arraytype'] == 'relStateArray':
              array2['copies'] = iters * array2['copies']
            else:
              for access in array2['access'].values():
                access['reads'] *= iters
                access['writes'] *= iters
            array1 = getMember(linfo1[category], array2['name'])
            if array1:
              array1['access'] = mergeAccess(array1, array2)
            else:
              linfo1[category].append(array2)

      [loopvars, ranges, strides] = getRanges(linfo)
      if not (isType0(ranges) or isType1(ranges) or isType2(ranges)):
        if options.flag_warn:
          print >> sys.stderr, "Warning: ignoring loop nest: %s" % str(ranges)
        return

      # HACK: figure out loop type based on ranges
      sweeps = 1
      if isType1(ranges):
        # species sweep in outermost loop
        sweeps = numIters(linfo['range'])
        linfo = linfo['loops'][0]
        loopvars = loopvars[1:]
        ranges = ranges[1:]
        strides = strides[1:]
      elif isType2(ranges):
        # species sweep(s) inside third spatial loop level
        loopvars = loopvars[0:-1]
        ranges = ranges[0:-1]
        strides = strides[0:-1]
        linfo3 = linfo['loops'][0]['loops'][0]
        for linfo4 in linfo3['loops']:
          mergeLoops(linfo3, linfo4)
        del linfo3['loops']

      # use data in innermost loop
      linfo['loops'][0]['loops'][0]['linenum'] = linfo['linenum']
      linfo = linfo['loops'][0]['loops'][0]
      del linfo['loopvar']
      del linfo['range']
      del linfo['stride']
      loopvars.reverse()
      ranges.reverse()
      strides.reverse()
      linfo['loopvars'] = loopvars
      linfo['ranges'] = ranges
      linfo['strides'] = strides
      linfo['sweeps'] = sweeps

      return linfo

    # main analyzeInfo loop
    for (fname, finfo) in self.info.items():
      newLoops = []
      for (idx, linfo) in enumerate(finfo['loops']):
        linfo = flattenLoopNest(linfo)
        if not linfo:
          continue
        linfo['registers'] = getRegisters(linfo['scalars'], linfo['stateArrays'], linfo['arrays'])
        getAccessTypes(linfo)
        getPencilsAndPlanes(linfo)
        newLoops.append(linfo)
      finfo['loops'] = newLoops
      getLiveArrays(fname, finfo, memTableFlag = True)


  def printInfo(self):
    for (fname, finfo) in self.info.iteritems():
      print fname
      pp.pprint(self.info[fname])


def main(args):
  if len(args) < 3:
    print "Usage: <me> <staticxml-input> <statictsv-output>"
    sys.exit(1)

  (staticxml, statictsv) = args[1:3]

  sa = StaticAnalysis(staticxml)
  if options.flag_debug:
    sa.printInfo()

  grapher.Grapher(sa).generateDepGraphs()
  spreadsheet.SpreadsheetWriter(sa).exportSpreadsheet(statictsv)
 
if __name__ == '__main__':
  main(sys.argv)
