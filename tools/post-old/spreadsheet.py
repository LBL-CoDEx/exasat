import os
import re

from common import numIters, diff, numNonZero
import params
import analyze
import parser

TEMPFILE = 'output.temp.csv'

def getRef(ref):
  if type(ref) == type(''):
    return '$REF{%s}' % ref
  elif type(ref) == type(()):
    return '$REF(%d; %d)' % ref
  elif type(ref) == type(0):
    return '$REF(0; %d)' % ref

def getRefs(*args):
  return tuple(map(getRef, args))

def doRefSubs(expr):
  return params.doRefSubs_l(expr, getRef)

def excelRef(i, j, refType = 'N'):
  rowPrefix = ''
  colPrefix = ''
  if refType == 'E' or refType == 'C':
    colPrefix = '$'
  if refType == 'E' or refType == 'R':
    rowPrefix = '$'

  if i < 0 or j < 0:
    raise Exception('Spreadsheet reference out of range: (%d, %d)' % (i, j))
  (letter1, letter2) = divmod(j, 26)
  if letter1 == 0:
    return '%s%s%s%d' % (colPrefix, chr(ord('A') + letter2), rowPrefix, i+1)
  else:
    return '%s%s%s%d' % (colPrefix, chr(ord('A') + letter1 - 1) + chr(ord('A') + letter2), rowPrefix, i+1)

def processRefs(outfile, infile):
  
  fin = open(infile, 'rt')
  lines = fin.readlines()
  fin.close()

  # first pass: scan for labels
  labels = {}
  for (i, line) in enumerate(lines):
    cells = line.strip('\n\r').split('|')
    for (j, val) in enumerate(cells):
      m = re.search('\#\$([ECR])LABEL\{(.*?)\}', val)
      if m:
        (ltype, label) = m.groups()
        if label in labels and \
           (ltype == 'E' or \
            ltype == 'C' and j != labels[label][1] or \
            ltype == 'R' and i != labels[label][0]):
          raise Exception('duplicate label detected: %s at (%d, %d) and (%d, %d)' % \
                          (label, labels[label][0], labels[label][1], i, j))
        labels[label] = (i, j, ltype)

  # second pass: process and write file
  fout = open(outfile, 'wt')
  for (i, line) in enumerate(lines):
    cells = line.strip('\n\r').split('|')
    for (j, val) in enumerate(cells):

      # remove labels
      m = re.search('\#\$([ECR])LABEL\{(.*?)\}', val)
      if m:
        val = re.sub(re.escape(m.group(0)), '', val)

      # process relative references
      m = re.search('\$REF\((-?\d+);\s*(-?\d+)\)', val)
      while m:
        val = re.sub(re.escape(m.group(0)), \
                     excelRef(i + int(m.group(1)), j + int(m.group(2)), 'N'), \
                     val)
        m = re.search('\$REF\((-?\d+);\s*(-?\d+)\)', val)

      # process labeled references
      m = re.search('\$REF\{(.*?)\}', val)
      while m:
        ref = labels[m.group(1)]
        (r, c) = ref[0:2]
        if ref[2] == 'C':
          r = i
        elif ref[2] == 'R':
          c = j
        val = re.sub(re.escape(m.group(0)), excelRef(r, c, ref[2]), val)
        m = re.search('\$REF\{(.*?)\}', val)

      fout.write(val + '\t') # switch from pipes to tabs
    fout.write('\n')
  fout.close()


class SpreadsheetWriter(object):

  def __init__(self, sa):
    self.info = sa.info

  def exportSpreadsheet(self, outfile):

    def writeParameters(f, my_params):
      f.write('Parameters:\n')
      for (key, _) in params.model_params:
        value = my_params[key]
        f.write('%s|%s#$ELABEL{%s}\n' % (key, str(value), key))
      f.write('\n')

    def writeComm(f):
      f.write('COMM ANALYSIS\n')
      f.write('|||Iteration Space\n')
      f.write('function|line|array#$CLABEL{array}|' + \
              'X#$CLABEL{iterx}|Y#$CLABEL{itery}|Z#$CLABEL{iterz}|' + \
              'Total Off-Node (MB)#$CLABEL{totalOffNode}|Num Messages#$CLABEL{numMsgs}|' + \
              'BW Cost (s)#$CLABEL{BWCost}|Latency Cost (s)#$CLABEL{latencyCost}|Total Cost (s)\n')
      numArrays = 0
      for (fname, finfo) in self.info.iteritems():
        for comm in finfo['comms']:
          for array in comm['arrays']:
            numArrays += 1
            f.write('%s|%s|%s|' % (fname, comm['linenum'], array['name']))
            f.write('=%s|=%s|=%s|' % getRefs('X Problem Size', 'Y Problem Size', 'Z Problem Size'))

            numComps = array['numComps']
            f.write('=%s*%s*(%d*%s*%s+%d*%s*%s+%d*%s*%s)/2^20|' % (doRefSubs(numComps), getRef('Word Size'), \
              diff(array['ghost'][0])[0], getRef('itery'), getRef('iterz'), \
              diff(array['ghost'][1])[0], getRef('iterx'), getRef('iterz'), \
              diff(array['ghost'][2])[0], getRef('iterx'), getRef('itery')))
            f.write('%d|' % numNonZero(array['ghost']))
            f.write('=%s*2^20/(%s*2^30)|' % getRefs('totalOffNode', 'NIC BW (GB/s)'))
            f.write('=%s*%s/10^6|' % getRefs('numMsgs', 'NIC Latency (us)'))
            f.write('=%s+%s|' % getRefs('BWCost', 'latencyCost'))

            f.write('\n')

      # totals
      f.write('Totals' + '|'*6)
      f.write('=sum(%s:%s)|' % getRefs((-numArrays, 0), (-1, 0)) * 5)
      f.write('\n\n')

    def writeSummary(f):
      # write loop summary header with column labels
      f.write('LOOP ANALYSIS\n')
      f.write('|||Iteration Space|||Block Iteration Space|||||Flops/cell/sweep||||Registers||RW/W Arrays||RW/W Working Set|||Working set (naive)|||Working set (streaming writes)|||Working set (reuse only)|||Working set (actual)|||Memory Traffic (GBytes)|||||Computation (Gflops)|||Estimated execution times (s)\n')
      f.write('function|line|sweeps#$CLABEL{sweeps}|' + \
              'X#$CLABEL{iterx}|Y#$CLABEL{itery}|Z#$CLABEL{iterz}|' + \
              'X#$CLABEL{blockx}|Y#$CLABEL{blocky}|Z#$CLABEL{blockz}|' + \
              'X (cache line)#$CLABEL{blockxcl}|Num Blocks#$CLABEL{numBlocks}|' + \
              'add#$CLABEL{add}|mul#$CLABEL{mul}|div#$CLABEL{div}|special#$CLABEL{spec}|GP Regs|FP Regs|' + \
              'RW Arrays#$CLABEL{RWArrays}|W Arrays#$CLABEL{WArrays}|' + \
              'RW/W WS Planes#$CLABEL{WSWPlanes}|RW/W WS Pencils#$CLABEL{WSWPencils}|RW/W WS Cells#$CLABEL{WSWCells}|' + \
              'WS/plane/core (kB)#$CLABEL{WSPlane}|' + \
              'WS/pencil/core (kB)#$CLABEL{WSPencil}|' + \
              'WS/cell/core (kB)#$CLABEL{WSCell}|' + \
              'WS/plane/core (kB)#$CLABEL{WSStreamPlane}|' + \
              'WS/pencil/core (kB)#$CLABEL{WSStreamPencil}|' + \
              'WS/cell/core (kB)#$CLABEL{WSStreamCell}|' + \
              'WS/plane/core (kB)#$CLABEL{WSReusePlane}|' + \
              'WS/pencil/core (kB)#$CLABEL{WSReusePencil}|' + \
              'WS/cell/core (kB)#$CLABEL{WSReuseCell}|' + \
              'WS/plane/core (kB)#$CLABEL{WSActualPlane}|' + \
              'WS/pencil/core (kB)#$CLABEL{WSActualPencil}|' + \
              'WS/cell/core (kB)#$CLABEL{WSActualCell}|' + \
              'Reuse between planes#$CLABEL{blockGbytes}|' + \
              'Reuse between pencils#$CLABEL{planeGbytes}|' + \
              'Reuse within pencils#$CLABEL{pencilGbytes}|' + \
              'No reuse within pencils#$CLABEL{cellGbytes}|' + \
              'Actual#$CLABEL{aGbytes}|' + \
              'Gflops performed#$CLABEL{gflops}|Weighted Gflops#$CLABEL{wGflops}|B/F ratio|' + \
              'time (CPU)#$CLABEL{timeCPU}|time (DRAM)#$CLABEL{timeRAM}|time (CPU and DRAM)|\n')

      # write summary info
      for fname in sorted(self.info.keys()):
        function = self.info[fname]
        for loop in function['loops']:
          flops = loop['flops']
          loopid = '%s.%d' % (fname, loop['linenum'])

          f.write('%s|%d|=%s' % (fname, loop['linenum'], doRefSubs(loop['sweeps'])))

          # iteration space
          f.write('|=%s' % doRefSubs(numIters(loop['ranges'][0])))
          f.write('|=%s' % doRefSubs(numIters(loop['ranges'][1])))
          f.write('|=%s' % doRefSubs(numIters(loop['ranges'][2])))

          # block iteration space
          f.write('|=%s' % getRefs('X Block Size'))
          f.write('|=%s' % getRefs('Y Block Size'))
          f.write('|=%s' % getRefs('Z Block Size'))

          # block X rounded up to nearest cache line and num blocks
          f.write('|=ceiling(%s/(%s/%s),1)*(%s/%s)' % getRefs('blockx', 'Cache Line Size', 'Word Size', \
                                                                        'Cache Line Size', 'Word Size'))
          f.write('|=%s/%s*%s/%s*%s/%s' % getRefs('iterx', 'blockx', 'itery', 'blocky', 'iterz', 'blockz'))

          # flops
          f.write('|=%s|=%s|=%s|=%s' % tuple(map(lambda x: doRefSubs(flops[x]), \
                                                 ['adds', 'multiplies', 'divides', 'specials'])))

          # registers
          f.write('|=%s|=%s' % (doRefSubs(loop['registers']['ints'] + loop['registers']['ptrs']), \
                                doRefSubs(loop['registers']['floats'])))

          # read/write and write arrays
          numArrays = loop['BW']['numArrays']
          f.write('|=%s|=%s' % (doRefSubs(numArrays['RW']), doRefSubs(numArrays['W'])))

          # working set planes, pencils, and cells
          numPlanes = loop['WS']['numPlanes']
          numPencils = loop['WS']['numPencils']
          numCells = loop['WS']['numCells']
          f.write('|=%s' % doRefSubs(numPlanes['RW'] + numPlanes['W']))
          f.write('|=%s' % doRefSubs(numPencils['RW'] + numPencils['W']))
          f.write('|=%s' % doRefSubs(numCells['RW'] + numCells['W']))

          # working set sizes (kB)
          f.write('|=%s+%s*%s*%s*%s/2^10' % getRefs('%s.WSPlane' % loopid, 'Word Size', 'WSWPlanes', 'blockxcl', 'blocky'))
          f.write('|=%s+%s*%s*%s/2^10' % getRefs('%s.WSPencil' % loopid, 'Word Size', 'WSWPencils', 'blockxcl'))
          f.write('|=%s+%s*%s/2^10' % getRefs('%s.WSCell' % loopid, 'Word Size', 'WSWCells'))
          f.write('|=%s' % getRefs('%s.WSPlane' % loopid))
          f.write('|=%s' % getRefs('%s.WSPencil' % loopid))
          f.write('|=%s' % getRefs('%s.WSCell' % loopid))
          f.write('|=%s' % getRefs('%s.WSReusePlane' % loopid))
          f.write('|=%s' % getRefs('%s.WSReusePencil' % loopid))
          f.write('|=%s' % getRefs('%s.WSReuseCell' % loopid))
          f.write('|=IF(%s=%s,%s,IF(%s, %s, IF(%s, %s, %s)))' % \
                  getRefs('blockGbytes', 'planeGbytes', 'WSActualPencil', \
                          'NTA Hints', 'WSReusePlane', 'Streaming Writes', 'WSStreamPlane', 'WSPlane'))
          f.write('|=IF(%s=%s,%s,IF(%s, %s, IF(%s, %s, %s)))' % \
                  getRefs('planeGbytes', 'pencilGbytes', 'WSActualCell', \
                          'NTA Hints', 'WSReusePencil', 'Streaming Writes', 'WSStreamPencil', 'WSPencil'))
          f.write('|=IF(%s, %s, IF(%s, %s, %s))' % \
                  getRefs('NTA Hints', 'WSReuseCell', 'Streaming Writes', 'WSStreamCell', 'WSCell'))

          # GBytes transferred (use results from detailed read-only array analysis)
          f.write('|=%s*(%s+(%s*%s+%s*%s)*%s*%s*%s*%s/2^30)' % \
                  getRefs('sweeps', '%s.blockGbytes' % loopid, \
                          'RWArrays', 'RW cost', 'WArrays', 'W cost',
                          'numBlocks', 'blockxcl', 'blocky', 'blockz'))
          f.write('|=%s*(%s+(%s*%s+%s*%s)*%s*%s*%s*%s/2^30)' % \
                  getRefs('sweeps', '%s.planeGbytes' % loopid, \
                          'RWArrays', 'RW cost', 'WArrays', 'W cost',
                          'numBlocks', 'blockxcl', 'blocky', 'blockz'))
          f.write('|=%s*(%s+(%s*%s+%s*%s)*%s*%s*%s*%s/2^30)' % \
                  getRefs('sweeps', '%s.pencilGbytes' % loopid, \
                          'RWArrays', 'RW cost', 'WArrays', 'W cost',
                          'numBlocks', 'blockxcl', 'blocky', 'blockz'))
          f.write('|=%s*(%s+(%s*%s+%s*%s)*%s*%s*%s*%s/2^30)' % \
                  getRefs('sweeps', '%s.cellGbytes' % loopid, \
                          'RWArrays', 'RW cost', 'WArrays', 'W cost',
                          'numBlocks', 'blockxcl', 'blocky', 'blockz'))
          f.write('|=IF(%s<%s,MIN(%s:%s),IF(%s<%s,MIN(%s:%s),IF(%s<%s,MIN(%s:%s),%s)))' % \
                  getRefs('WSActualPlane', '$/thread group (kB)', 'blockGbytes', 'cellGbytes', \
                          'WSActualPencil', '$/thread group (kB)', 'planeGbytes', 'cellGbytes', \
                          'WSActualCell', '$/thread group (kB)', 'pencilGbytes', 'cellGbytes', \
                          'cellGbytes'))

          # Gflops and arithmetic intensity
          f.write('|=%s*(%s+%s+%s+%s)*%s*%s*%s/1e9' % getRefs('sweeps', 'add', 'mul', 'div', 'spec', \
                                                              'iterx', 'itery', 'iterz'))
          f.write('|=%s*(%s+%s+%s*%s+%s*%s)*%s*%s*%s/1e9' % \
                  getRefs('sweeps', 'add', 'mul', 'Division Cost', 'div', 'Special Cost', 'spec', \
                          'iterx', 'itery', 'iterz'))
          f.write('|=(%s*2^30)/(%s*10^9)' % getRefs('aGbytes', 'wGflops'))

          # estimated execution times
          f.write('|=%s/(%s*%s)' % getRefs('wGflops', 'Gflop/s/thread', 'Threads'))
          f.write('|=%s/(%s*%s)' % getRefs('aGbytes', 'GB/s/thread', 'Threads'))
          f.write('|=max(%s:%s)' % getRefs('timeCPU', 'timeRAM'))

          f.write('\n')

      # totals
      numLoops = sum([len(func['loops']) for func in self.info.values()])
      f.write('Total/Max' + '|'*30)
      f.write('|=max(%s:%s)' % getRefs((-numLoops, 0), (-1, 0)) * 3)
      f.write('|=sum(%s:%s)' % getRefs((-numLoops, 0), (-1, 0)) * 7)
      f.write('|=(%s*2^30)/(%s*10^9)' % getRefs('aGbytes', 'wGflops'))
      f.write('|=sum(%s:%s)' % getRefs((-numLoops, 0), (-1, 0)) * 3)

      f.write('\n\n')

    def writePWSummary(f):
      # write loop summary header with column labels
      f.write('WOODWARD LOOP ANALYSIS\n')
      f.write('|||Iteration Space|||Block Iteration Space|||||Flops/cell/sweep||||Registers||Num Blocks||||Working Set|Bandwidth||||Computation (Gflops)|||Estimated execution times (s)\n')
      f.write('function|line|sweeps#$CLABEL{sweeps}|' + \
              'X#$CLABEL{iterx}|Y#$CLABEL{itery}|Z#$CLABEL{iterz}|' + \
              'X#$CLABEL{blockx}|Y#$CLABEL{blocky}|Z#$CLABEL{blockz}|' + \
              'X (cache line)#$CLABEL{blockxcl}|Num Blocks#$CLABEL{numBlocks}|' + \
              'add#$CLABEL{add}|mul#$CLABEL{mul}|div#$CLABEL{div}|special#$CLABEL{spec}|GP Regs|FP Regs|' + \
              'Resident#$CLABEL{ResidentBlocks}|' + \
              'R#$CLABEL{RBlocks}|' + \
              'W#$CLABEL{WBlocks}|' + \
              'RW#$CLABEL{RWBlocks}|' + \
              'kB#$CLABEL{WSBlock}|' + \
              'R (GB)#$CLABEL{BWBlockR}|' + \
              'W (GB)#$CLABEL{BWBlockW}|' + \
              'RW (GB)#$CLABEL{BWBlockRW}|' + \
              'Total (GB)#$CLABEL{BWBlock}|' + \
              'Gflops performed#$CLABEL{gflopsPW}|Weighted Gflops#$CLABEL{wGflopsPW}|B/F ratio|' + \
              'time (CPU)#$CLABEL{timeCPUPW}|time (DRAM)#$CLABEL{timeRAMPW}|time (CPU and DRAM)|\n')

      # write summary info
      for fname in sorted(self.info.keys()):
        function = self.info[fname]
        for loop in function['loops']:
          flops = loop['flops']
          loopid = '%s.%d' % (fname, loop['linenum'])

          f.write('%s|%d|=%s' % (fname, loop['linenum'], doRefSubs(loop['sweeps'])))

          # iteration space
          f.write('|=%s' % doRefSubs(numIters(loop['ranges'][0])))
          f.write('|=%s' % doRefSubs(numIters(loop['ranges'][1])))
          f.write('|=%s' % doRefSubs(numIters(loop['ranges'][2])))

          # expanded block iteration space with overlapped ghost regions
          f.write('|=%s+(%s-%s)' % getRefs('X Block Size', 'iterx', 'X Problem Size'))
          f.write('|=%s+(%s-%s)' % getRefs('Y Block Size', 'itery', 'Y Problem Size'))
          f.write('|=%s+(%s-%s)' % getRefs('Z Block Size', 'iterz', 'Z Problem Size'))

          # block X rounded up to nearest cache line and num blocks
          f.write('|=ceiling(%s/(%s/%s),1)*(%s/%s)' % getRefs('blockx', 'Cache Line Size', 'Word Size', \
                                                                        'Cache Line Size', 'Word Size'))
          f.write('|=%s/%s*%s/%s*%s/%s' % getRefs('X Problem Size', 'X Block Size', \
                                                  'Y Problem Size', 'Y Block Size', \
                                                  'Z Problem Size', 'Z Block Size'))

          # flops
          f.write('|=%s|=%s|=%s|=%s' % tuple(map(lambda x: doRefSubs(flops[x]), \
                                                 ['adds', 'multiplies', 'divides', 'specials'])))

          # registers
          f.write('|=%s|=%s' % (doRefSubs(loop['registers']['ints'] + loop['registers']['ptrs']), \
                                doRefSubs(loop['registers']['floats'])))

          # number of blocks resident, read, and written
          f.write('|=%s'*4 % tuple(map(doRefSubs, [loop['WS']['numBlocks'], loop['BW']['numBlocks']['R'], \
                                                   loop['BW']['numBlocks']['W'], loop['BW']['numBlocks']['RW']])))

          # working set and bandwidth estimates
          f.write('|=%s*(%s)/2^10' % (getRef('Word Size'), doRefSubs(loop['WS']['sizeBlocks'])))
          f.write('|=IF(%s<%s, %s*%s*(%s)/2^30, "N/A")' % \
            (getRefs('WSBlock', '$/thread group (kB)', 'R cost', 'numBlocks') + \
            (doRefSubs(loop['BW']['sizeBlocks']['R']),)))
          f.write('|=IF(%s<%s, %s*%s*(%s)/2^30, "N/A")' % \
            (getRefs('WSBlock', '$/thread group (kB)', 'W cost', 'numBlocks') + \
            (doRefSubs(loop['BW']['sizeBlocks']['W']),)))
          f.write('|=IF(%s<%s, %s*%s*(%s)/2^30, "N/A")' % \
            (getRefs('WSBlock', '$/thread group (kB)', 'RW cost', 'numBlocks') + \
            (doRefSubs(loop['BW']['sizeBlocks']['RW']),)))
          f.write('|=%s+%s+%s' % getRefs('BWBlockR', 'BWBlockW', 'BWBlockRW'))

          # Gflops and arithmetic intensity
          f.write('|=%s*(%s+%s+%s+%s)*%s*%s*%s*%s/1e9' % getRefs('sweeps', 'add', 'mul', 'div', 'spec', \
                                                              'blockx', 'blocky', 'blockz', 'numBlocks'))
          f.write('|=%s*(%s+%s+%s*%s+%s*%s)*%s*%s*%s*%s/1e9' % \
                  getRefs('sweeps', 'add', 'mul', 'Division Cost', 'div', 'Special Cost', 'spec', \
                          'blockx', 'blocky', 'blockz', 'numBlocks'))
          f.write('|=(%s*2^30)/(%s*10^9)' % getRefs('BWBlock', 'wGflopsPW'))

          # estimated execution times
          f.write('|=%s/(%s*%s)' % getRefs('wGflopsPW', 'Gflop/s/thread', 'Threads'))
          f.write('|=%s/(%s*%s)' % getRefs('BWBlock', 'GB/s/thread', 'Threads'))
          f.write('|=max(%s:%s)' % getRefs('timeCPUPW', 'timeRAMPW'))

          f.write('\n')

      # totals
      numLoops = sum([len(func['loops']) for func in self.info.values()])
      f.write('Total/Max' + '|'*20)
      f.write('|=max(%s:%s)' % getRefs((-numLoops, 0), (-1, 0)))
      f.write('|=sum(%s:%s)' % getRefs((-numLoops, 0), (-1, 0)) * 6)
      f.write('|=(%s*2^30)/(%s*10^9)' % getRefs('BWBlock', 'wGflopsPW'))
      f.write('|=sum(%s:%s)' % getRefs((-numLoops, 0), (-1, 0)) * 3)

      f.write('\n\n')

    def writeDetails(f):
      f.write('DETAILED READ-ONLY ARRAY INFO\n\n')
      for fname in sorted(self.info.keys()):
        f.write('Function:|%s\n' % (fname))
        function = self.info[fname]
        for loop in function['loops']:

          # print header
          loopid = '%s.%d' % (fname, loop['linenum'])
          f.write('Loop line num:|%d\n' % (loop['linenum']))
          f.write('|||Iteration Space|||Block Iteration Space|||||Block Access Space|||||Bandwidth|||' + \
                  'Working Set|||Reuse WS|||WS (all reads)|||WS (reuse only)|||GBytes transferred/sweep\n')
          f.write('||Name|X|Y|Z|X|Y|Z|X (cache line)|Num Blocks|' + \
                  'X#$CLABEL{accessx}|Y#$CLABEL{accessy}|Z#$CLABEL{accessz}|X (cache line)#$CLABEL{accessxcl}|' + \
                  'Copies#$CLABEL{aBWCopies}|Planes#$CLABEL{aBWPlanes}|' + \
                  'Pencils#$CLABEL{aBWPencils}|Cells#$CLABEL{aBWCells}|' + \
                  'Planes#$CLABEL{aWSPlanes}|Pencils#$CLABEL{aWSPencils}|Cells#$CLABEL{aWSCells}|' + \
                  'Planes#$CLABEL{aWSReusePlanes}|Pencils#$CLABEL{aWSReusePencils}|Cells#$CLABEL{aWSReuseCells}|' + \
                  'WS/plane/core (kB)#$CLABEL{aWSPlane}|' + \
                  'WS/pencil/core (kB)#$CLABEL{aWSPencil}|' + \
                  'WS/cell/core (kB)#$CLABEL{aWSCell}|' + \
                  'WS/plane/core (kB)#$CLABEL{aWSReusePlane}|' + \
                  'WS/pencil/core (kB)#$CLABEL{aWSReusePencil}|' + \
                  'WS/cell/core (kB)#$CLABEL{aWSReuseCell}|' + \
                  'Reuse between planes#$CLABEL{aBlockGbytes}|' + \
                  'Reuse between pencils#$CLABEL{aPlaneGbytes}|' + \
                  'Reuse within pencils#$CLABEL{aPencilGbytes}|' + \
                  'No reuse within pencils#$CLABEL{aCellGbytes}|' + \
                  '\n')

          readArrays = sorted(analyze.getR(loop['arrays']), key=lambda x:x['name'])
          for array in readArrays:

            # name
            f.write('||%s' % (array['name']))

            # iteration space
            f.write('|=%s' % doRefSubs(numIters(loop['ranges'][0])))
            f.write('|=%s' % doRefSubs(numIters(loop['ranges'][1])))
            f.write('|=%s' % doRefSubs(numIters(loop['ranges'][2])))

            # block iteration space (with also X rounded up to nearest cache line) and num blocks
            f.write('|=%s' % getRefs('X Block Size'))
            f.write('|=%s' % getRefs('Y Block Size'))
            f.write('|=%s' % getRefs('Z Block Size'))
            f.write('|=ceiling(%s/(%s/%s),1)*(%s/%s)' % getRefs('blockx', 'Cache Line Size', 'Word Size', \
                                                                          'Cache Line Size', 'Word Size'))
            f.write('|=%s/%s*%s/%s*%s/%s' % getRefs('iterx', 'blockx', 'itery', 'blocky', 'iterz', 'blockz'))

            # block access space (with also X rounded up to nearest cache line)
            f.write('|=%s+%d' % (getRef('blockx'), diff(array['ghost'][0])[0]))
            f.write('|=%s+%d' % (getRef('blocky'), diff(array['ghost'][1])[0]))
            f.write('|=%s+%d' % (getRef('blockz'), diff(array['ghost'][2])[0]))
            f.write('|=ceiling(%s/(%s/%s),1)*(%s/%s)' % getRefs('accessx', 'Cache Line Size', 'Word Size', \
                                                                           'Cache Line Size', 'Word Size'))

            # bandwidth and working set figures from access pattern analysis
            f.write('|=%s'*10 % tuple(map(doRefSubs, [array['copies'], array['BW']['numPlanes'], \
                    array['BW']['numPencils'], array['BW']['numCells'], array['WS']['numPlanes'], \
                    array['WS']['numPencils'], array['WS']['numCells'], array['WS']['numReusePlanes'], \
                    array['WS']['numReusePencils'], array['WS']['numReuseCells']])))

            # working set size formulas (kB)

            if array['stenciltype'] == 'faces':
              f.write('|=%s*%s*((%s-1)*(%s*%s)+1*(%s*%s+%s*%s-%s*%s))/2^10' % \
                  getRefs('Word Size', 'aBWCopies', \
                          'aWSPlanes', 'blockxcl', 'blocky', \
                          'blockxcl', 'accessy', \
                          'accessxcl', 'blocky', \
                          'blockxcl', 'blocky'))
            else:
              f.write('|=%s*%s*%s*%s*%s/2^10' % \
                  getRefs('Word Size', 'aBWCopies', 'aWSPlanes', 'accessxcl', 'accessy'))

            if array['stenciltype'] == 'faces':
              f.write('|=%s*%s*((%s-1)*(%s)+1*(%s))/2^10' % \
                  getRefs('Word Size', 'aBWCopies', 'aWSPencils', 'blockxcl', 'accessxcl'))
            else:
              f.write('|=%s*%s*%s*%s/2^10' % getRefs('Word Size', 'aBWCopies', 'aWSPencils', 'accessxcl'))

            f.write('|=%s*%s*%s/2^10' % getRefs('Word Size', 'aBWCopies', 'aWSCells'))

            # working set size formulas (reuse only) (kB)

            if array['stenciltype'] == 'faces':
              f.write('|=%s*%s*((%s-1)*(%s*%s)+1*(%s*%s+%s*%s-%s*%s))/2^10' % \
                  getRefs('Word Size', 'aBWCopies', \
                          'aWSReusePlanes', 'blockxcl', 'blocky', \
                          'blockxcl', 'accessy', \
                          'accessxcl', 'blocky', \
                          'blockxcl', 'blocky'))
            else:
              f.write('|=%s*%s*%s*%s*%s/2^10' % \
                  getRefs('Word Size', 'aBWCopies', 'aWSReusePlanes', 'accessxcl', 'accessy'))

            if array['stenciltype'] == 'faces':
              f.write('|=%s*%s*((%s-1)*(%s)+1*(%s))/2^10' % \
                  getRefs('Word Size', 'aBWCopies', 'aWSReusePencils', 'blockxcl', 'accessxcl'))
            else:
              f.write('|=%s*%s*%s*%s/2^10' % getRefs('Word Size', 'aBWCopies', 'aWSReusePencils', 'accessxcl'))

            f.write('|=%s*%s*%s/2^10' % getRefs('Word Size', 'aBWCopies', 'aWSReuseCells'))

            # memory traffic formulas (GBytes per sweep)
            if array['stenciltype'] == 'faces':
              f.write('|=%s*%s*%s*(%s*%s*%s+%s*%s*%s+%s*%s*%s-2*%s*%s*%s)/2^30' % \
                  getRefs('numBlocks', 'aBWCopies', 'R cost', 'accessxcl', 'blocky', 'blockz', 'blockxcl', \
                          'accessy', 'blockz', 'blockxcl', 'blocky', 'accessz', 'blockxcl', 'blocky', 'blockz'))
              f.write('|=%s*%s*%s*((%s-1)*(%s*%s)+1*(%s*%s+%s*%s-%s*%s))*%s/2^30' % \
                  getRefs('numBlocks', 'aBWCopies', 'R cost', \
                          'aBWPlanes', 'blockxcl', 'blocky', \
                          'accessxcl', 'blocky', \
                          'blockxcl', 'accessy', \
                          'blockxcl', 'blocky', \
                          'blockz'))
              f.write('|=%s*%s*%s*((%s-1)*%s+1*%s)*%s*%s/2^30' % \
                  getRefs('numBlocks', 'aBWCopies', 'R cost', \
                          'aBWPencils', 'blockxcl', \
                          'accessxcl', \
                          'blocky', 'blockz'))
            else:
              f.write('|=%s*%s*%s*%s*%s*%s/2^30' % getRefs('numBlocks', 'aBWCopies', 'R cost', \
                                                           'accessxcl', 'accessy', 'accessz'))
              f.write('|=%s*%s*%s*%s*%s*%s*%s/2^30' % getRefs('numBlocks', 'aBWCopies', 'aBWPlanes', 'R cost', \
                                                              'accessxcl', 'accessy', 'blockz'))
              f.write('|=%s*%s*%s*%s*%s*%s*%s/2^30' % getRefs('numBlocks', 'aBWCopies', 'aBWPencils', 'R cost', \
                                                              'accessxcl', 'blocky', 'blockz'))
            f.write('|=%s*%s*%s*%s*%s*%s*%s/2^30' % getRefs('numBlocks', 'aBWCopies', 'aBWCells', 'R cost', \
                                                            'blockx', 'blocky', 'blockz'))

            f.write('\n')

          def addFname(*args):
            return tuple(map(lambda x: '%s.%s' % (loopid, x), args))

          temp = '|=sum(%%s:%%s)#$ELABEL{%s}' * 10 % \
                 addFname('WSPlane', 'WSPencil', 'WSCell', \
                          'WSReusePlane', 'WSReusePencil', 'WSReuseCell', \
                          'blockGbytes', 'planeGbytes', 'pencilGbytes', 'cellGbytes')

          f.write('||Total' + '|'*22 + temp % (getRefs((-max(len(readArrays), 1), 0), (-1, 0)) * 10))
          f.write('\n')

        f.write('\n')

    def writeRWTable(f):

      f.write('ARRAY READ-WRITE TABLES\n\n')
      for (fname, finfo) in self.info.iteritems():
        # gather array names
        arrays = {}
        for liveness in finfo['liveness']:
          for (aname, (x, copies)) in liveness.iteritems():
            arrays[aname] = copies
        numArrays = len(arrays)

        # write header
        f.write('%s|Copies#$CLABEL{rwCopies}|' % fname)
        for linfo in [finfo['loops'][i] for i in finfo['execorder']]:
          f.write('%d|' % linfo['linenum'])
        f.write('Reads|Writes|Present|\n')

        # write table
        numLoops = len(finfo['liveness'])
        for aname in sorted(arrays.keys()):
          f.write('%s|=%s|' % (aname, doRefSubs(arrays[aname])))
          for liveness in finfo['liveness']:
            if aname in liveness:
              f.write(liveness[aname][0])
            f.write('|')
          f.write('=%s*(COUNTIF(%s:%s,"R")+COUNTIF(%s:%s,"RW"))|' % \
                  getRefs('rwCopies', -numLoops, -1, -numLoops, -1))
          f.write('=%s*(COUNTIF(%s:%s,"W")+COUNTIF(%s:%s,"RW")+COUNTIF(%s:%s,"WR"))|' % \
                  getRefs('rwCopies', -numLoops-1, -2, -numLoops-1, -2, -numLoops-1, -2))
          f.write('=%s*COUNTIF(%s:%s,"P")|' % getRefs('rwCopies', -numLoops-2, -3))
          f.write('\n')

        # write footprints
        f.write('Footprint||')
        for idx in xrange(0, numLoops):
          f.write('=SUMIF(%s:%s,"<>",%s:%s)|' % getRefs((-numArrays,0), (-1,0), (-numArrays,-idx-1), (-1,-idx-1)))
        f.write('\n\n')

    def writeSVTable(f):

      def idxToStr(aname, idx):
        return '%s%s' % (ainfo['name'], idx)
      def strToIdx(aname):
        (aname, idx) = aname.split('(')
        return (aname, tuple(eval('(' + idx)))

      f.write('STATE VARIABLE TABLES\n\n')
      for (fname, finfo) in self.info.iteritems():
        # gather array names
        SVs = set()
        copies = {}
        for linfo in finfo['loops']:
          for sinfo in linfo['scalars']:
            SVs.add(sinfo['name'])
            copies[sinfo['name']] = 1
          for ainfo in linfo['stateArrays']:
            if ainfo['arraytype'] == 'relStateArray':
              aname = ainfo['name']
              SVs.add(aname)
              if aname in copies and copies[aname] != ainfo['copies']:
                raise Exception('numCopies mismatch')
              copies[aname] = ainfo['copies']
            else:
              for idx in ainfo['access'].keys():
                aname = idxToStr(ainfo['name'], idx)
                SVs.add(aname)
                copies[aname] = 1
        SVs = sorted(SVs)

        # write header
        f.write('%s|copies|' % fname)
        for linfo in [finfo['loops'][i] for i in finfo['execorder']]:
          f.write('%d|' % linfo['linenum'])
        f.write('\n')

        # write table
        for vname in SVs:
          f.write('%s|=%s|' % (vname, doRefSubs(copies[vname])))
          if re.search('\(', vname):
            (vname, idx) = strToIdx(vname)
          for linfo in [finfo['loops'][i] for i in finfo['execorder']]:
            # check if variable is a scalar
            scalar = analyze.getMember(linfo['scalars'], vname)
            if scalar:
              f.write('=%s' % doRefSubs(scalar['reads'] + scalar['writes']))
            # check if variable is a state array
            array = analyze.getMember(linfo['stateArrays'], vname)
            if array and array['arraytype'] == 'relStateArray':
              f.write('=%s' % doRefSubs(array['access'][(0,)]['reads'] + array['access'][(0,)]['writes']))
            if array and array['arraytype'] == 'absStateArray' and idx in array['access']:
              f.write('=%s' % doRefSubs(array['access'][idx]['reads'] + array['access'][idx]['writes']))
            f.write('|')
          f.write('\n')
        f.write('\n')

    # load parameters
    problem_xml = os.getenv('problem_xml')
    machine_xml = os.getenv('machine_xml')
    my_params = params.defaultParams()
    if problem_xml:
      params.customizeParams(my_params, parser.MachineXMLParser(problem_xml).params)
    if machine_xml:
      params.customizeParams(my_params, parser.MachineXMLParser(machine_xml).params)

    f = open(TEMPFILE, 'wt')
    writeParameters(f, my_params)
    writeComm(f)
    writeSummary(f)
    writePWSummary(f)
    writeDetails(f)
    writeRWTable(f)
    writeSVTable(f)
    f.close()

    processRefs(outfile, TEMPFILE)
    os.remove(TEMPFILE)
