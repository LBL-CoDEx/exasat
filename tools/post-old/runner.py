import string

import os
import math
import params
import parser
from params import model_params, doParamSubs
from common import numIters, diff, mymax, numNonZero, prod
import analyze

def writeRegAllocTables(f, result, dataType):
  f.write('REGISTER ALLOCATION TABLE\n\n')
  for (fname, func) in result['functions'].iteritems():
    regAlloc = func['regAlloc'][dataType]
    f.write('%s\n' % fname)
    f.write('\tHits\tMisses\tTotal\n')
    f.write(('State \t'+'%d\t'*3+'\n') % (regAlloc['state']['hits'], regAlloc['state']['misses'], \
                                          regAlloc['state']['hits'] + regAlloc['state']['misses']))
    f.write(('Stream\t'+'%d\t'*3+'\n') % (regAlloc['stream']['hits'], regAlloc['stream']['misses'], \
                                          regAlloc['stream']['hits'] + regAlloc['stream']['misses']))
    f.write('\n\n')
    for loop in func['loops']:
      regAlloc = loop['regAlloc'][dataType]
      f.write('Line:\t%d\n' % (loop['linenum']))
      f.write('\tIn Registers\tIn Cache\tTotal\tHits\tMisses\tTotal\n')
      f.write(('State \t'+'%d\t'*6+'\n') % (regAlloc['state']['inRegs'], regAlloc['state']['inCache'], \
                                            regAlloc['state']['inRegs'] + regAlloc['state']['inCache'], \
                                            regAlloc['state']['hits'], regAlloc['state']['misses'], \
                                            regAlloc['state']['hits'] + regAlloc['state']['misses']))
      f.write(('Stream\t'+'%d\t'*6+'\n') % (regAlloc['stream']['inRegs'], regAlloc['stream']['inCache'], \
                                            regAlloc['stream']['inRegs'] + regAlloc['stream']['inCache'], \
                                            regAlloc['stream']['hits'], regAlloc['stream']['misses'], \
                                            regAlloc['stream']['hits'] + regAlloc['stream']['misses']))
      f.write('\nAccess Count Scatterplot Graph\n')
      table = regAlloc['state']['table']
      if table:
        f.write(     '\t' + string.join(map(lambda x: str(x[0]), table), '\t') + '\n')
        f.write('State\t' + string.join(map(lambda x: str(x[1]), table), '\t') + '\n')
      table = regAlloc['stream']['table']
      if table:
        f.write(      '\t' + string.join(map(lambda x: str(x[0]), table), '\t') + '\n')
        f.write('Stream\t' + string.join(map(lambda x: str(x[1]), table), '\t') + '\n')
      f.write('\n')
    f.write('\n')
  f.write('\n')


class Runner(object):

  def __init__(self, sa):
    self.info = sa.info

  def runModel(self, param_overrides = None):

    # compute state and streaming variable register allocation
    def regAllocModel(linfo, my_params):

      # split available registers between state and stream variables
      def splitRegs(table, numRegs):
        numStateVars = len(table['state'])
        numStreamVars = len(table['stream'])
        numStateRegs = 0
        numStreamRegs = 0
        if numStateVars == 0:
          numStreamRegs = min(numStreamVars, numRegs)
        elif numStreamVars == 0:
          numStateRegs = min(numStateVars, numRegs)
        else:
          while numStateRegs + numStreamRegs < numRegs and \
                numStateRegs + numStreamRegs < numStateVars + numStreamVars:
            if numStateRegs >= numStateVars or \
               table['stream'][numStreamRegs] > table['state'][numStateRegs]:
              numStreamRegs += 1
            else:
              numStateRegs += 1
        return {'state': numStateRegs, 'stream': numStreamRegs}

      # count accesses of each variable type
      result = {'double': {'state': {}, 'stream': {}}, 'int': {'state': {}, 'stream': {}}}
      for (dataType, numRegs) in [('double', my_params['FP Regs/thread']), ('int', my_params['Int Regs/thread'])]:
        table = {'state': [], 'stream': []}
        for sinfo in filter(lambda x: x['datatype'] == dataType, linfo['scalars']):
          table['state'].append(doParamSubs(sinfo['reads'] + sinfo['writes'], my_params))
        for (dest, source) in [(table['state' ], linfo['stateArrays']), \
                               (table['stream'], linfo['arrays'     ])]:
          for ainfo in filter(lambda x: x['datatype'] == dataType, source):
            for idx in ainfo['access'].keys():
              # TODO: merge paramsubs
              dest.extend([doParamSubs(ainfo['access'][idx]['reads'] + \
                                       ainfo['access'][idx]['writes'], my_params)] * \
                          int(doParamSubs(ainfo['copies'], my_params)))
        table['state'] = sorted(table['state'], reverse=True)
        table['stream'] = sorted(table['stream'], reverse=True)

        # determine register allocation between state and streaming variables
        numRegs = splitRegs(table, numRegs)

        # count number of register "hits" and "misses" (i.e. register vs. L1 cache accesses)
        # and summarize access table for graphing
        for varType in ['state', 'stream']:
          t = table[varType]
          nr = numRegs[varType]
          nv = len(t)
          hits = 0
          misses = 0
          if t:
            hits   = sum(t[:nr])
            misses = sum(t[ nr:])
            t = [(1, t[0])] + \
                [(i+1, t[i]) for i in xrange(1,nv-1) if \
                 t[i-1] != t[i] or t[i] != t[i+1]] + \
                [(nv, t[-1])]
          result[dataType][varType]['inRegs' ] = nr
          result[dataType][varType]['inCache'] = nv - nr
          result[dataType][varType]['hits'   ] = hits
          result[dataType][varType]['misses' ] = misses
          result[dataType][varType]['table'  ] = t

        # must load streaming values into registers on each iteration
        # first access is a miss, not a hit
        result[dataType]['stream']['misses'] += result[dataType]['stream']['inRegs']
        result[dataType]['stream']['hits'  ] -= result[dataType]['stream']['inRegs']

      return result

    def runLoopModel(linfo, my_params):
      result = {}

      result['linenum'] = linfo['linenum']
      result['range'] = (doParamSubs(numIters(linfo['ranges'][0]), my_params), \
                         doParamSubs(numIters(linfo['ranges'][1]), my_params), \
                         doParamSubs(numIters(linfo['ranges'][2]), my_params))
      blockx = my_params['X Block Size']
      blocky = my_params['Y Block Size']
      blockz = my_params['Z Block Size']
      result['block'] = [blockx, blocky, blockz]
      result['numBlocks'] = float(prod(result['range'])) / prod(result['block'])

      # round blockx up to nearest cache line multiple
      CLWords = my_params['Cache Line Size'] / my_params['Word Size']
      blockx = math.ceil(float(blockx) / CLWords) * CLWords

      # registers
      result['GPRegs'] = doParamSubs(linfo['registers']['ints'] + linfo['registers']['ptrs'], my_params)
      result['FPRegs'] = doParamSubs(linfo['registers']['floats'], my_params)
      result['regAlloc'] = regAllocModel(linfo, my_params)

      # Compute working sets and memory traffic for read-only arrays
      arrays = []
      for ainfo in analyze.getR(linfo['arrays']):
        array = {}
        array['name'] = ainfo['name']
        array['access'] = map(lambda x, y:x + diff(y)[0], result['block'], ainfo['ghost'])
        accessx = array['access'][0]
        accessy = array['access'][1]
        accessz = array['access'][2]

        # round accessx up to nearest cache line multiple
        accessx = math.ceil(float(accessx) / CLWords) * CLWords

        # WS calculation for generic case
        array['WS'] = {'all'  : {'plane' : my_params['Word Size'] * ainfo['WS']['numPlanes'] * \
                                           accessx * accessy, \
                                 'pencil': my_params['Word Size'] * ainfo['WS']['numPencils'] * \
                                           accessx, \
                                 'cell'  : my_params['Word Size'] * ainfo['WS']['numCells'], \
                                }, \
                       'reuse': {'plane' : my_params['Word Size'] * ainfo['WS']['numReusePlanes'] * \
                                           accessx * accessy, \
                                 'pencil': my_params['Word Size'] * ainfo['WS']['numReusePencils'] * \
                                           accessx, \
                                 'cell'  : my_params['Word Size'] * ainfo['WS']['numReuseCells'], \
                                }, \
                      }
        # fix the plane WS for faces-only stencils
        if ainfo['stenciltype'] == 'faces':

          array['WS']['all']['plane'] = my_params['Word Size'] * \
              ((ainfo['WS']['numPlanes'] - 1) * (blockx * blocky) + \
               1 * (blockx * accessy + accessx * blocky - blockx * blocky))
          array['WS']['all']['pencil'] = my_params['Word Size'] * \
              ((ainfo['WS']['numPencils'] - 1) * blockx + 1 * accessx)

          array['WS']['reuse']['plane'] = my_params['Word Size'] * \
              ((ainfo['WS']['numReusePlanes'] - 1) * (blockx * blocky) + \
               1 * (blockx * accessy + accessx * blocky - blockx * blocky))
          array['WS']['reuse']['pencil'] = my_params['Word Size'] * \
              ((ainfo['WS']['numReusePencils'] - 1) * blockx + 1 * accessx)

        # BW calculation for generic case
        array['BW'] = {'block' : result['numBlocks'] * ainfo['BW']['numCopies'] * my_params['R cost'] * \
                                 accessx * accessy * accessz, \
                       'plane' : result['numBlocks'] * ainfo['BW']['numPlanes'] * my_params['R cost'] * \
                                 accessx * accessy * blockz, \
                       'pencil': result['numBlocks'] * ainfo['BW']['numPencils'] * my_params['R cost'] * \
                                 accessx * blocky * blockz, \
                       'cell'  : result['numBlocks'] * ainfo['BW']['numCells'] * my_params['R cost'] * \
                                 blockx * blocky * blockz, \
                      }
        # fix the block and plane BW for faces-only stencils
        if ainfo['stenciltype'] == 'faces':
          array['BW']['block'] = result['numBlocks'] * ainfo['BW']['numCopies'] * my_params['R cost'] * \
              (accessx * blocky * blockz + blockx * accessy * blockz + \
               blockx * blocky * accessz - 2 * blockx * blocky * blockz)
          array['BW']['plane'] = result['numBlocks'] * my_params['R cost'] * \
              ((ainfo['BW']['numPlanes'] - 1) * (blockx * blocky) + \
               1 * (accessx * blocky + blockx * accessy - blockx * blocky)) * blockz
          array['BW']['pencil'] = result['numBlocks'] * my_params['R cost'] * \
              ((ainfo['BW']['numPencils'] - 1) * blockx + 1 * accessx) * blocky * blockz

        arrays.append(array)

      sumWSAllPlane = sum(map(lambda x: x['WS']['all']['plane'], arrays))
      sumWSAllPencil = sum(map(lambda x: x['WS']['all']['pencil'], arrays))
      sumWSAllCell = sum(map(lambda x: x['WS']['all']['cell'], arrays))
      sumWSReusePlane = sum(map(lambda x: x['WS']['reuse']['plane'], arrays))
      sumWSReusePencil = sum(map(lambda x: x['WS']['reuse']['pencil'], arrays))
      sumWSReuseCell = sum(map(lambda x: x['WS']['reuse']['cell'], arrays))

      sumBWBlock = sum(map(lambda x: x['BW']['block'], arrays))
      sumBWPlane = sum(map(lambda x: x['BW']['plane'], arrays))
      sumBWPencil = sum(map(lambda x: x['BW']['pencil'], arrays))
      sumBWCell = sum(map(lambda x: x['BW']['cell'], arrays))

      result['arrays'] = arrays

      # Compute working sets for different scenarios
      result['WS'] = {'all'   : {'plane' : sumWSAllPlane + my_params['Word Size'] * \
                                           (linfo['WS']['numPlanes']['RW'] + linfo['WS']['numPlanes']['W']) * \
                                           blockx * blocky, \
                                 'pencil': sumWSAllPencil + my_params['Word Size'] * \
                                           (linfo['WS']['numPencils']['RW'] + linfo['WS']['numPencils']['W']) * \
                                           blockx, \
                                 'cell'  : sumWSAllCell + my_params['Word Size'] * \
                                           (linfo['WS']['numCells']['RW'] + linfo['WS']['numCells']['W']), \
                                }, \
                      'stream': {'plane' : sumWSAllPlane, \
                                 'pencil': sumWSAllPencil, \
                                 'cell'  : sumWSAllCell, \
                                }, \
                      'reuse' : {'plane' : sumWSReusePlane, \
                                 'pencil': sumWSReusePencil, \
                                 'cell'  : sumWSReuseCell, \
                                }, \
                     }

      # Determine actual working sets based on cache utilization policy
      if my_params['NTA Hints']:
        result['WS']['actual'] = result['WS']['reuse']
      elif my_params['Streaming Writes']:
        result['WS']['actual'] = result['WS']['stream']
      else:
        result['WS']['actual'] = result['WS']['all']

      # Compute memory traffic for different reuse scenarios
      numSweeps = doParamSubs(linfo['sweeps'], my_params)
      RWWBW = (linfo['BW']['numArrays']['RW'] * my_params['RW cost'] + \
               linfo['BW']['numArrays']['W' ] * my_params['W cost' ]) * result['numBlocks'] * blockx * blocky * blockz
      result['BW'] = {'block' : numSweeps * (sumBWBlock  + RWWBW), \
                      'plane' : numSweeps * (sumBWPlane  + RWWBW), \
                      'pencil': numSweeps * (sumBWPencil + RWWBW), \
                      'cell'  : numSweeps * (sumBWCell   + RWWBW), \
                     }

      # do symbolic parameter substitutions
      for x in result['WS'].values():
        for (y, z) in x.iteritems():
          x[y] = doParamSubs(z, my_params)
      for (y, z) in result['BW'].iteritems():
        result['BW'][y] = doParamSubs(z, my_params)

      # bandwidth should be no worse than model prediction for cases with worse reuse,
      #   but model approximations cause different inaccuracies for different cases
      result['BW']['pencil'] = min(result['BW']['pencil'], result['BW']['cell'])
      result['BW']['plane'] = min(result['BW']['plane'], result['BW']['pencil'])
      result['BW']['block'] = min(result['BW']['block'], result['BW']['plane'])

      # if there's no difference between memory traffic, working set is effectively reduced
      if result['BW']['pencil'] == result['BW']['cell']:
        result['WS']['actual']['cell'] = 0
      if result['BW']['plane'] == result['BW']['pencil']:
        result['WS']['actual']['pencil'] = result['WS']['actual']['cell']
      if result['BW']['block'] == result['BW']['plane']:
        result['WS']['actual']['plane'] = result['WS']['actual']['pencil']

      # Compute "actual" memory traffic based on type of reuse given available cache
      if result['WS']['actual']['plane'] <= my_params['$/thread group (kB)'] * 2**10:
        result['BW']['actual'] = result['BW']['block']
      elif result['WS']['actual']['pencil'] <= my_params['$/thread group (kB)'] * 2**10:
        result['BW']['actual'] = result['BW']['plane']
      elif result['WS']['actual']['cell'] <= my_params['$/thread group (kB)'] * 2**10:
        result['BW']['actual'] = result['BW']['pencil']
      else:
        result['BW']['actual'] = result['BW']['cell']

      # save final WS and BW
      result['WSFinal'] = result['WS']['actual']['plane']
      result['BWFinal'] = result['BW']['actual']

      # number of flops and weighted flops
      result['adds'      ] = numSweeps * prod(result['range']) * doParamSubs(linfo['flops']['adds'      ], my_params)
      result['multiplies'] = numSweeps * prod(result['range']) * doParamSubs(linfo['flops']['multiplies'], my_params)
      result['divides'   ] = numSweeps * prod(result['range']) * doParamSubs(linfo['flops']['divides'   ], my_params)
      result['specials'  ] = numSweeps * prod(result['range']) * doParamSubs(linfo['flops']['specials'  ], my_params)
      result['flops'] = result['adds'] + result['multiplies'] + result['divides'] + result['specials']
      result['wflops'] = result['adds'] + result['multiplies'] + my_params['Division Cost'] * result['divides'] + \
                                                                 my_params['Special Cost'] * result['specials']

      # arithmetic intensity
      if result['wflops'] != 0:
        result['BF'] = float(result['BWFinal']) / result['wflops']
      else:
        result['BF'] = float('nan')

      # execution time
      result['cputime'] = float(result['wflops']) / \
                          (my_params['Gflop/s/thread'] * my_params['Threads'] * 10**9)
      result['ramtime'] = float(result['BW']['actual']) / \
                          (my_params['GB/s/thread'] * my_params['Threads'] * 2**30)
      # assume perfect overlap
      if result['cputime'] > result['ramtime']:
        result['ramtime'] = 0
      else:
        result['cputime'] = 0
      result['time'] = max(result['cputime'], result['ramtime'])

      return result

    def runPWLoopModel(linfo, my_params):
      result = {}

      result['linenum'] = linfo['linenum']
      numSweeps = doParamSubs(linfo['sweeps'], my_params)
      result['range'] = (doParamSubs(numIters(linfo['ranges'][0]), my_params), \
                         doParamSubs(numIters(linfo['ranges'][1]), my_params), \
                         doParamSubs(numIters(linfo['ranges'][2]), my_params))
      blockx = my_params['X Block Size'] + (result['range'][0] - my_params['X Problem Size'])
      blocky = my_params['Y Block Size'] + (result['range'][1] - my_params['Y Problem Size'])
      blockz = my_params['Z Block Size'] + (result['range'][2] - my_params['Z Problem Size'])
      result['block'] = [blockx, blocky, blockz]
      numBlocks = (my_params['X Problem Size'] * my_params['Y Problem Size'] * my_params['Z Problem Size']) / \
                   (my_params['X Block Size'] * my_params['Y Block Size'] * my_params['Z Block Size'])
      result['numBlocks'] = numBlocks

      # registers
      result['GPRegs'] = doParamSubs(linfo['registers']['ints'] + linfo['registers']['ptrs'], my_params)
      result['FPRegs'] = doParamSubs(linfo['registers']['floats'], my_params)
      result['regAlloc'] = regAllocModel(linfo, my_params)

      # save working set and memory traffic
      result['WSFinal'] = doParamSubs(my_params['Word Size'] * linfo['WS']['sizeBlocks'], my_params)
      if result['WSFinal'] <= my_params['$/thread group (kB)'] * 2**10:
        result['BWFinal'] = numBlocks * doParamSubs(my_params[ 'R cost'] * linfo['BW']['sizeBlocks']['R'] + \
                                                    my_params[ 'W cost'] * linfo['BW']['sizeBlocks']['W'] + \
                                                    my_params['RW cost'] * linfo['BW']['sizeBlocks']['RW'], my_params)
      else:
        # if not enough cache, punt on this method
        result['BWFinal'] = float('inf')

      # number of flops and weighted flops
      numCellIters = numSweeps * numBlocks * prod(result['block'])
      result['adds'      ] = numCellIters * doParamSubs(linfo['flops']['adds'      ], my_params)
      result['multiplies'] = numCellIters * doParamSubs(linfo['flops']['multiplies'], my_params)
      result['divides'   ] = numCellIters * doParamSubs(linfo['flops']['divides'   ], my_params)
      result['specials'  ] = numCellIters * doParamSubs(linfo['flops']['specials'  ], my_params)
      result['flops'] = result['adds'] + result['multiplies'] + result['divides'] + result['specials']
      result['wflops'] = result['adds'] + result['multiplies'] + my_params['Division Cost'] * result['divides'] + \
                                                                 my_params['Special Cost'] * result['specials']

      # arithmetic intensity
      if result['wflops'] != 0:
        result['BF'] = float(result['BWFinal']) / result['wflops']
      else:
        result['BF'] = float('nan')

      # execution time
      result['cputime'] = float(result['wflops']) / \
                          (my_params['Gflop/s/thread'] * my_params['Threads'] * 10**9)
      result['ramtime'] = float(result['BWFinal']) / \
                          (my_params['GB/s/thread'] * my_params['Threads'] * 2**30)
      # assume perfect overlap
      if result['cputime'] > result['ramtime']:
        result['ramtime'] = 0
      else:
        result['cputime'] = 0
      result['time'] = max(result['cputime'], result['ramtime'])

      return result

    def runCommModel(comm, my_params):

      def computeBW(ainfo, my_params):
        xy = my_params['X Problem Size'] * my_params['Y Problem Size']
        xz = my_params['X Problem Size'] * my_params['Z Problem Size']
        yz = my_params['Y Problem Size'] * my_params['Z Problem Size']
        xyGhost = diff(ainfo['ghost'][0])[0]
        xzGhost = diff(ainfo['ghost'][1])[0]
        yzGhost = diff(ainfo['ghost'][2])[0]
        return my_params['Word Size'] * doParamSubs(ainfo['numComps'], my_params) * (xyGhost * xy + xzGhost * xz + yzGhost * yz)

      result = {}
      arrays = []
      for ainfo in filter(lambda x: x['type'] == 'ghost', comm['arrays']):
        array = {'name': ainfo['name'], 'BW': computeBW(ainfo, my_params), 'messages': numNonZero(ainfo['ghost'])}
        arrays.append(array)
      result['arrays'] = arrays
      result['BW'] = sum(map(lambda x: x['BW'], arrays))
      result['messages'] = sum(map(lambda x: x['messages'], arrays))
      result['time'] = float(result['BW']) / (my_params['NIC BW (GB/s)'] * 2**30) + \
                       float(result['messages']) * my_params['NIC Latency (us)'] * 1e-6
      return result

    def runFuncModel(func, my_params):
      comms = []
      for cinfo in func['comms']:
        comms.append(runCommModel(cinfo, my_params))
      loops = []
      for linfo in func['loops']:
        if 'PWFlag' in my_params and my_params['PWFlag']:
          loops.append(runPWLoopModel(linfo, my_params))
        else:
          loops.append(runLoopModel(linfo, my_params))
      result = {'comms': comms, \
                'loops': loops, \
                'cputime' : sum(map(lambda x: x['cputime'], loops)), \
                'ramtime' : sum(map(lambda x: x['ramtime'], loops)), \
                'time' : sum(map(lambda x: x['time'], loops)), \
                'adds': sum(map(lambda x: x['adds'], loops)), \
                'multiplies': sum(map(lambda x: x['multiplies'], loops)), \
                'divides': sum(map(lambda x: x['divides'], loops)), \
                'specials': sum(map(lambda x: x['specials'], loops)), \
                'flops': sum(map(lambda x: x['flops'], loops)), \
                'wflops': sum(map(lambda x: x['wflops'], loops)), \
                'BW'   : sum(map(lambda x: x['BWFinal'], loops)), \
                'WS'   : mymax(map(lambda x: x['WSFinal'], loops)), \
                'GPRegs': mymax(map(lambda x: x['GPRegs'], loops)), \
                'FPRegs': mymax(map(lambda x: x['FPRegs'], loops)), \
                'regAlloc': {}, \
                'cTime': sum(map(lambda x: x['time'], comms)), \
                'cBW'  : sum(map(lambda x: x['BW'], comms)), \
                'cMsgs': sum(map(lambda x: x['messages'], comms)), \
               }

      # summarize register hits and misses
      for dataType in ['int', 'double']:
        result['regAlloc'][dataType] = {}
        for varType in ['state', 'stream']:
          result['regAlloc'][dataType][varType] = {}
          for tag in ['hits', 'misses']:
            result['regAlloc'][dataType][varType][tag] = sum(map(lambda x: x['regAlloc'][dataType][varType][tag], loops))

      return result

    # load parameters
    problem_xml = os.getenv('problem_xml')
    machine_xml = os.getenv('machine_xml')
    my_params = params.defaultParams()
    if problem_xml:
      params.customizeParams(my_params, parser.MachineXMLParser(problem_xml).params)
    if machine_xml:
      params.customizeParams(my_params, parser.MachineXMLParser(machine_xml).params)
    params.customizeParams(my_params, param_overrides)

    functions = {}
    for (fname, func) in self.info.iteritems():
      functions[fname] = runFuncModel(func, my_params)
    result = {'params': my_params, \
              'functions': functions, \
              'cputime' : sum(map(lambda x: x['cputime'], functions.values())), \
              'ramtime' : sum(map(lambda x: x['ramtime'], functions.values())), \
              'time' : sum(map(lambda x: x['time'], functions.values())), \
              'adds': sum(map(lambda x: x['adds'], functions.values())), \
              'multiplies': sum(map(lambda x: x['multiplies'], functions.values())), \
              'divides': sum(map(lambda x: x['divides'], functions.values())), \
              'specials': sum(map(lambda x: x['specials'], functions.values())), \
              'flops': sum(map(lambda x: x['flops'], functions.values())), \
              'wflops': sum(map(lambda x: x['wflops'], functions.values())), \
              'BW'   : sum(map(lambda x: x['BW'], functions.values())), \
              'WS'   : mymax(map(lambda x: x['WS'], functions.values())), \
              'GPRegs': mymax(map(lambda x: x['GPRegs'], functions.values())), \
              'FPRegs': mymax(map(lambda x: x['FPRegs'], functions.values())), \
              'regAlloc': {}, \
              'cTime': sum(map(lambda x: x['cTime'], functions.values())), \
              'cBW'  : sum(map(lambda x: x['cBW'], functions.values())), \
              'cMsgs': sum(map(lambda x: x['cMsgs'], functions.values())), \
             }

    # summarize register hits and misses
    for dataType in ['int', 'double']:
      result['regAlloc'][dataType] = {}
      for varType in ['state', 'stream']:
        result['regAlloc'][dataType][varType] = {}
        for tag in ['hits', 'misses']:
          result['regAlloc'][dataType][varType][tag] = sum(map(lambda x: x['regAlloc'][dataType][varType][tag], \
                                                               functions.values()))
      f = open('RegAllocTables.%s.temp.tsv' % dataType, 'wt')
      writeRegAllocTables(f, result, dataType)
      f.close()

    self.result = result
    return result
