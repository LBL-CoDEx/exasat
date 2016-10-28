#!/usr/bin/python
from math import sqrt
import sys
sys.path.append('/home/dunat/Projects/exasat/tools/post')

from analyzeXML import StaticAnalysis
import pprint

def runCompute(filename, config, problemSize): 

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

  totaltime = commtime + comptime 
  teraflops = ((numflops)/ totaltime / 1000)
  mpoints = (float((problemSize * problemSize * problemSize )) / totaltime /10**6)

#finite diff
  print "%g " % numflops,
  print "%g " % onnodetraffic, 
  print "%g " % cputime, 
  print "%g " % dramtime,
  print "%g " % comptime,

#communication
  print "%g " % offnodetraffic, 
  print "%g " % NICbwtime,
  print "%g " % NIClattime,
  print "%g " % commtime,

  print "%g " % totaltime,
  print "%g " % teraflops,
  print "%g " % mpoints,  
  print 

def runExaNode1():
#exaNode1 Settings 

  xmlpath = '/home/dunat/Projects/exasat/examples/xmls/'
  s3d = True

  filename = xmlpath + "diffterm_compact.xml"#"advance-smc-modified.xml"
  version = "diffterm-compact"

  if s3d:
    filename = xmlpath + "diffterm_s3d.xml"
    version = "diffterm-S3D"

  pSize = 128
  species = [9, 21, 53, 71, 107]

  #naive
  blockSize = [16, 16,  16, 8,  8, 8, # for 9 , 21, 53   
               8, 4, 8, 4 ] # for 71 and 107

  #s3d
  if s3d:
    blockSize = [32, 32, 32, 32, 
                 32, 32, 
                 32, 32, 
                 32, 32]
                                 
             
  print "Results for exaNode1"

  print "SMC Problem Size Info Finite Difference Part , , Communication Part"
  print "Version BoxSize BlockSize NSpecies",
  print "NumberofGFlops On-NodeTraffic(GB) ",
  print "CPUTime(s) DRAMTime(s) ComputeTime(s)",
  print "Off-NodeTraffic(GB)",
  print "NICBWTime(s) NICLatTime(s) CommTime(s)",
  print "TotalTime(s) Teraflops/s MegaPoints/s"

  dummy = 0        
  for nspec in species:

    problemSize = pSize
    print version, problemSize, "%dX%d" % (blockSize[dummy], blockSize[dummy + 1]) , 
    print nspec, 
      
    config = {}
    config['Num Species'] = nspec
    config['Gflop/s/thread'] = 10 
    config['Threads']= 1024 
    config['GB/s/thread'] = 1
    config['$/thread group (kB)']= 1024
    config['NIC BW (GB/s)']=100
    config['NIC Latency (us)']=0.4

    config['Division Cost'] = 40
    config['Special Cost'] = 125
    config['X Problem Size'] = problemSize
    config['Y Problem Size'] = problemSize
    config['Z Problem Size'] = problemSize
    config['X Block Size'] = blockSize[dummy]
    config['Y Block Size'] = blockSize[dummy + 1]
    config['Z Block Size'] = problemSize
    
    runCompute(filename, config, problemSize)
    dummy = dummy + 2

#program starts here


runExaNode1()

#runDemo()
#printBlockCacheTables()

