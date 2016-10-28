#!/usr/bin/python

from math import sqrt
from analyzeXML import StaticAnalysis
import pprint


def runCompute(filename, config, problemSize, BF): 

  sa = StaticAnalysis(filename)

  config['X Problem Size'] = problemSize
  config['Y Problem Size'] = problemSize
  config['Z Problem Size'] = problemSize
  config['X Blocking Factor'] = BF
  config['Y Blocking Factor'] = BF
  config['Z Blocking Factor'] = 1
  
  r = sa.runModel(config)
#  print "Params: %s" % str(r['params'])
  print "Numberof Flops %g" % (float(r['flops'])/10**9)
#  print "Perf Gflops/s %g" % ((float(r['flops']) / 10**9)/r['time'])
#  print "Perf MPoints/s %g" % (float((problemSize * problemSize * problemSize )) /r['time'] /10**6)
  print "CPU Time(s) %g " % (float(r['flops']) / ((r['params']['Gflop/s/core'] * r['params']['Cores/chip'])) / 10**9)
  print "DRAM Time(s) %g " % (((float(r['BW']) / r['params']['GB/s/chip'])) * 2**(-30))
  print "Comp Time(s) %g " % r['time']
#  print "NICBW Time(s) %g " % ((float(r['cBW']) / r['params']['NIC BW (GB/s)']) * 2**(-30))
#  print "NICLat Time(s) %g " % (float(r['cMsgs']) * r['params']['NIC Latency (us)'] * 1e-6)
#  print "Comm Time(s) %g " % r['cTime']
#  print "Total Time(s) %g " % (float(r['time']) + r['cTime'])
#  print 

def runComm(filename, config, problemSize, topo):

  sa = StaticAnalysis(filename)

  config['X Problem Size'] = problemSize * topo[0]
  config['Y Problem Size'] = problemSize * topo[1]
  config['Z Problem Size'] = problemSize * topo[2]
  
  r = sa.runModel(config)
  print "NICBW Time(s) %g " % ((float(r['cBW']) / r['params']['NIC BW (GB/s)']) * 2**(-30))
  print "NICLat Time(s) %g " % (float(r['cMsgs']) * r['params']['NIC Latency (us)'] * 1e-6)
  print "Comm Time(s) %g " % r['cTime']
  #  print "Total Time(s) %g " % (float(r['time']) + r['cTime'])
  print 

  
def runExaNode1(filename, version, BF):
#exaNode1 Settings 
  numBoxes = [1,8,16,32,64]

  topology = [[1,1,1], [2,2,2], [4,2,2], [4,4,2], [4,4,4]]                             

  problemSize = 64

  print "Results for exaNode1"
  print "Box Size %g" % problemSize 
  print "CNS Version %s" % version 

  dummy = 0
                               
  for nBoxes in numBoxes:
    print "Num Boxes ", nBoxes 
    print "Box Topology ", topology[dummy]
    config = {}
    config['GB/s/chip'] = 4096 / nBoxes
    config['Gflop/s/core'] = 10 
    config['Cores/chip']= 1024 / nBoxes
    config['$/core (MB)']= 1
    config['NIC BW (GB/s)']=400
    config['NIC Latency (us)']=0.02
    runCompute(filename, config, problemSize, BF)                               
    runComm(filename, config, problemSize, topology[dummy])
    dummy = dummy + 1

#blocking factors are for 128KB cache per core 
#  runNodeConfig(filenameFused, config, 64, 1, 'fused-simple')
#  runNodeConfig(filenameNaive, config, 128, 11, 'naive')
#  runNodeConfig(filenameFused, config, 128, 1, 'fused-simple')
  
  
#program starts here

filenameNaive = '../../examples/xmls/advance-flat.xml'
filenameFused = '../../examples/xmls/advance-fused-nomod.xml'

version = "naive"
BF=6
runExaNode1(filenameNaive, version, BF)

version = "simple-fusion"
BF=1
runExaNode1(filenameFused, version, BF)



#runDemo()
#printBlockCacheTables()

