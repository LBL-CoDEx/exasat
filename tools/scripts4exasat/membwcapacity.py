#!/usr/bin/python

from math import sqrt
import sys
sys.path.append('/home/dunat/Projects/exasat/tools/post')

from analyzeXML import StaticAnalysis
import pprint

def runCompute(filename, chemfile, config, problemSize): 

  sa = StaticAnalysis(filename)
  
  r = sa.runModel(config)
  numflops = (float(r['flops'])/10**9)
  wnumflops = (float(r['wflops'])/10**9)
  onnodetraffic = ((float(r['BW']))/ 2**30)
  offnodetraffic = ((float(r['cBW']))* 2**(-30))

  cputime = (wnumflops / ((r['params']['Gflop/s/thread'] * r['params']['Threads'])))
  dramtime = (onnodetraffic / (r['params']['GB/s/thread']* r['params']['Threads']))
  comptime = r['time']
  NICbwtime = ((offnodetraffic / r['params']['NIC BW (GB/s)']))
  NIClattime = (float(r['cMsgs']) * r['params']['NIC Latency (us)'] * 1e-6)
  commtime = r['cTime']

#running chemistry 
  config['X Block Size'] = problemSize
  config['Y Block Size'] = problemSize
  config['Z Block Size'] = problemSize

  sa = StaticAnalysis(chemfile)
  r = sa.runModel(config)

  chemnumflops = (float(r['flops'])/10**9)
  wchemnumflops = (float(r['wflops'])/10**9)
  chemonnodetraffic = ((float(r['BW']))/ 2**30)
  chemcputime = (chemnumflops / ((r['params']['Gflop/s/thread'] * r['params']['Threads'])))
  chemdramtime = (chemonnodetraffic / (r['params']['GB/s/thread']* r['params']['Threads']))
  chemcomptime = r['time']

  totaltime = commtime + comptime + chemcomptime 
  teraflops = ((numflops + chemnumflops)/ totaltime/ 1000)
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
  print "%g " % teraflops,
  print "%g " % mpoints,  
  print 

def runExaNode1(version):
#exaNode1 Settings 

  xmlpath = '/home/dunat/Projects/exasat/examples/xmls/'
  filename = xmlpath + "advance-smc-modified.xml"

  pSize = [128, 128, 176, 176, 176, 224, 224, 224, 272, 272, 272, 352, 352, 352, 448, 448, 448, 576, 576, 576]
  membw = [4,   2,   4,   2,   1,   2,   1, 0.5, 2, 1, 0.5, 1, 0.5, 0.25, 0.5, 0.25,0.125,0.5,0.25,0.125]
#  membw = [4,    4,   2,   4,   2,   1,   2,   1,  0.5, 1, 0.5, 0.5]
#  pSize = [576,448,448,352,352,352,272,272,272,224,224,176]
#  memcapacity = [1024, 512, 512, 256, 256, 128, 128, 128, 64, 64, 32]
#  membw = [0.5, 0.5, 1, 0.5, 1, 2, 1, 2, 4, 2, 4, 4]
  species = 53

  #naive
  blockSize = [8, 8] #[16, 16,  16, 8,  8, 8, # for 9 , 21, 53   
              # 8, 4,  8, 4,  8, 4,   # 60, 70, 80
              # 8,4, 8,4]
  print "Results for exaNode1"
  print "SMC Version %s" % version 

  print "Problem Size Info Finite Difference Part , , Chemistry Part , , , Communication Part"
  print "BoxSize BlockSize NSpecies",
  print "NumberofGFlops On-NodeTraffic(GB) ",
  print "CPUTime(s) DRAMTime(s) ComputeTime(s)",
  print "NumberofGFlops On-NodeTraffic(GB) ",
  print "CPUTime(s) DRAMTime(s) ComputeTime(s)",
  print "Off-NodeTraffic(GB)",
  print "NICBWTime(s) NICLatTime(s) CommTime(s)",
  print "TotalTime(s) Teraflop/s MegaPoints/s"

  p = 0 
  dummy = 0
#  pim = True 

  for bw in membw:

    nspec = species 

    problemSize = pSize[p]
    print problemSize, "%dX%d" % (blockSize[dummy], blockSize[dummy + 1]) , 
    print nspec, 

    config = {}
    config['Num Species'] = nspec
    config['Gflop/s/thread'] = 10 
    config['Threads']= 1024 
    config['GB/s/thread'] = float(bw) * 1024/ config['Threads']
    config['$/thread group (kB)']= 1024
    config['NIC BW (GB/s)']=100
    config['NIC Latency (us)']=0.4

    if pim == True : 
      config['Gflop/s/thread'] = 10 
      config['Threads']= 64
      #config['$/thread group (kB)']= 0

    config['Division Cost'] = 20
    config['Special Cost']= 42
    config['X Problem Size'] = problemSize
    config['Y Problem Size'] = problemSize
    config['Z Problem Size'] = problemSize
    config['X Block Size'] = blockSize[dummy]
    config['Y Block Size'] = blockSize[dummy + 1]
    config['Z Block Size'] = problemSize

    if nspec == 9 :
      chemfilename = xmlpath + 'lidryer.xml'

    if nspec == 21 :
      chemfilename = xmlpath + 'drm19.xml'
      
    if nspec >= 53 : 
      chemfilename = xmlpath + 'grimech30.xml'

    runCompute(filename, chemfilename, config, problemSize)

    #dummy = dummy + 2
    p = p + 1

  
#program starts here



version = "baseline"
runExaNode1(version)

#runDemo()
#printBlockCacheTables()

