#!/usr/bin/python

from math import sqrt
import sys
sys.path.append('/home/dunat/Projects/exasat/tools/post')
from analyzeXML import StaticAnalysis
import pprint

def runCompute(filename, chemfile, config, problemSize): 

  sa = StaticAnalysis(filename)
  
  r = sa.runModel(config)
  numflops = (float(r['wflops'])/10**9)
  onnodetraffic = ((float(r['BW']))/ 2**30)
  offnodetraffic = ((float(r['cBW']))* 2**(-30))

  cputime = (numflops / ((r['params']['Gflop/s/thread'] * r['params']['Threads'])))
  dramtime = (onnodetraffic / (r['params']['GB/s/thread']* r['params']['Threads']))
  comptime = r['time']
  NICbwtime = ((offnodetraffic / r['params']['NIC BW (GB/s)']))
  NIClattime = (float(r['cMsgs']) * r['params']['NIC Latency (us)'] * 1e-6)
  commtime = r['cTime']

#running chemistry 
  sa = StaticAnalysis(chemfile)
  r = sa.runModel(config)

  config['X Block Size'] = problemSize
  config['Y Block Size'] = problemSize
  config['Z Block Size'] = problemSize

  chemnumflops = (float(r['wflops'])/10**9)
  chemonnodetraffic = ((float(r['BW']))/ 2**30)
  chemcputime = (chemnumflops / ((r['params']['Gflop/s/thread'] * r['params']['Threads'])))
  chemdramtime = (chemonnodetraffic / (r['params']['GB/s/thread']* r['params']['Threads']))
  chemcomptime = r['time']

  totaltime = commtime + comptime + chemcomptime 
  gflops = ((numflops + chemnumflops)/ totaltime)
  mpoints = (float((problemSize * problemSize * problemSize )) / totaltime /10**6)

#finite diff
  print "%g " % numflops,
  print "%g " % onnodetraffic, 
  print "%g " % cputime, 
  print "%g " % dramtime,
  print "%g " % comptime,

#chemistry 
  print "%g " % chemnumflops,
  print "%g " % chemonnodetraffic, 
  print "%g " % chemcputime, 
  print "%g " % chemdramtime,
  print "%g " % chemcomptime,

#communication
  print "%g " % offnodetraffic, 
  print "%g " % NICbwtime,
  print "%g " % NIClattime,
  print "%g " % commtime,

  print "%g " % totaltime,
  print "%g " % gflops,
  print "%g " % mpoints,  
  print 

def runExaNode1():
#exaNode1 Settings 
  filename = ['/home/dunat/Projects/exasat/examples/xmls/advance-smc-modified.xml',
              '/home/dunat/Projects/exasat/examples/xmls/comm/advance-smc-modified-s3d-comm.xml']

  versions = ['baseline-compact', 'baseline-s3d']

  blockSize = [16, 16,  16, 8,  8, 8, # for 9 , 21, 53   
               8, 4,  8, 4,  8, 4,   # 60, 70, 80
               8,4, 8,4]

#  problemSize = [32, 64, 128, 256, 512]
  species = [9, 21, 53, 60, 70, 80, 90, 100]
  problemSize = 64
  #naive

  print "Results for exaNode1"
  print "SMC Problem Size Info Finite Difference Part , , Chemistry Part , , , Communication Part"
  print "Version BoxSize BlockSize NSpecies",
  print "NumberofGFlops On-NodeTraffic(GB) ",
  print "CPUTime(s) DRAMTime(s) ComputeTime(s)",
  print "NumberofGFlops On-NodeTraffic(GB) ",
  print "CPUTime(s) DRAMTime(s) ComputeTime(s)",
  print "Off-NodeTraffic(GB)",
  print "NICBWTime(s) NICLatTime(s) CommTime(s)",
  print "TotalTime(s) Gflops/s MegaPoints/s"

  ver = 0 
  for version in versions:

    fname = filename[ver]

    for nspec in species:
      dummy = 0
      print version, problemSize, "%dX%d" % (blockSize[dummy], blockSize[dummy + 1]) , 
      print nspec, 

      config = {}
      config['Num Species'] = nspec
      config['GB/s/thread'] = 1
      config['Gflop/s/thread'] = 10 
      config['Threads']= 1024 
      config['$/thread group (kB)']= 1024
      config['NIC BW (GB/s)']=100
      config['NIC Latency (us)']=0.4
      
      config['X Problem Size'] = problemSize
      config['Y Problem Size'] = problemSize
      config['Z Problem Size'] = problemSize
      config['X Block Size'] = blockSize[dummy]
      config['Y Block Size'] = blockSize[dummy + 1]
      config['Z Block Size'] = problemSize

      if nspec == 9 :
        chemfilename = '/home/dunat/Projects/exasat/examples/xmls/lidryer.xml'

      if nspec == 21 :
        chemfilename = '/home/dunat/Projects/exasat/examples/xmls/drm19.xml'
          
      if nspec >= 53 : 
        chemfilename = '/home/dunat/Projects/exasat/examples/xmls/grimech30.xml'

      runCompute(fname, chemfilename, config, problemSize)

      dummy = dummy + 2
      
  ver  = ver + 1
  print

#program starts here

runExaNode1()

#runDemo()
#printBlockCacheTables()

