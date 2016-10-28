#!/usr/bin/python

from math import sqrt
from analyzeXML import StaticAnalysis
import pprint


def runCompute(filename, config, problemSize, BS): 

  sa = StaticAnalysis(filename)

  config['X Problem Size'] = problemSize
  config['Y Problem Size'] = problemSize
  config['Z Problem Size'] = problemSize
  config['X Blocking Factor'] = BS
  config['Y Blocking Factor'] = BS
  config['Z Blocking Factor'] = 1
  
  r = sa.runModel(config)
#  print "Params: %s" % str(r['params'])
  print "Numberof Flops %g" % (float(r['flops'])/10**9)
  print "On-node Traffic %g" % ((float(r['BW']))/ 2**30) 
  print "Off-node Traffic %g" % ((float(r['cBW']))* 2**(-30))
#  print "Perf Gflops/s %g" % ((float(r['flops']) / 10**9)/r['time'])
#  print "Perf MPoints/s %g" % (float((problemSize * problemSize * problemSize )) /r['time'] /10**6)
  print "CPU Time(s) %g " % (float(r['flops']) / ((r['params']['Gflop/s/core'] * r['params']['Cores/chip'])) / 10**9)
  print "DRAM Time(s) %g " % (((float(r['BW']) / r['params']['GB/s/chip'])) * 2**(-30))
  print "Comp Time(s) %g " % r['time']
  print "NICBW Time(s) %g " % ((float(r['cBW']) / r['params']['NIC BW (GB/s)']) * 2**(-30))
  print "NICLat Time(s) %g " % (float(r['cMsgs']) * r['params']['NIC Latency (us)'] * 1e-6)
  print "Comm Time(s) %g " % r['cTime']
  print "Total Time(s) %g " % (float(r['time']) + r['cTime'])
  print 

def runExaNode1(filename, version):
#exaNode1 Settings 

  problemSize = [32, 64, 128, 256, 512]
  #naive
  BS=[1, 1, 1, 1, 1]
  #BF = [3, 6, 16, 32]
  #simple fuse
  #BF = [2, 3, 6, 12]

  print "Results for exaNode1"
  print "SMC Version %s" % version 

  dummy = 0
                               
  for pSize in problemSize:
    print "Box Size", pSize
    print "Blocking Factor", BF[dummy]
    config = {}
    config['GB/s/thread'] = 1
    config['Gflop/s/thread'] = 10 
    config['Threads']= 1024 
    config['$/thread group (kB)']= 1024
    config['NIC BW (GB/s)']=100
    config['NIC Latency (us)']=0.4
    runCompute(filename, config, pSize, BF[dummy])
    dummy = dummy + 1

#blocking factors are for 128KB cache per core 
#  runNodeConfig(filenameFused, config, 64, 1, 'fused-simple')
#  runNodeConfig(filenameNaive, config, 128, 11, 'naive')
#  runNodeConfig(filenameFused, config, 128, 1, 'fused-simple')
  
  
#program starts here

filename = '../../examples/xmls/advance-smc-modified.xml'
#filename = '../../examples/xmls/advance-fused-nomod.xml'

version = "baseline"
runExaNode1(filename, version)

#version = "simple-fusion"
#BF=1
#runExaNode1(filenameFused, version, BF)


#runDemo()
#printBlockCacheTables()

