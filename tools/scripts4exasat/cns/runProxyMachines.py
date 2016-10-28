#!/usr/bin/python

from math import sqrt
from analyzeXML import StaticAnalysis
import pprint


def runNodeConfig(filename, config, problemSize, BF, version): 

  sa = StaticAnalysis(filename)

  config['X Problem Size'] = problemSize
  config['Y Problem Size'] = problemSize
  config['Z Problem Size'] = problemSize
  config['X Blocking Factor'] = BF
  config['Y Blocking Factor'] = BF
  config['Z Blocking Factor'] = 1
  
  r = sa.runModel(config)
#  print "Params: %s" % str(r['params'])
  print "Perf Gflops/s %g" % ((float(r['flops']) / 10**9)/(r['time']+r['cTime']))
  print "Perf MPoints/s %g" % (float((problemSize * problemSize * problemSize )) /(r['time']+r['cTime']) /10**6)
  print "CPU Time(s) %g " % (float(r['flops']) / ((r['params']['Gflop/s/core'] * r['params']['Cores/chip'])) / 10**9)
  print "DRAM Time(s) %g " % (((float(r['BW']) / r['params']['GB/s/chip'])) * 2**(-30))
  print "Comp Time(s) %g " % r['time']
  print "NICBW Time(s) %g " % ((float(r['cBW']) / r['params']['NIC BW (GB/s)']) * 2**(-30))
  print "NICLat Time(s) %g " % (float(r['cMsgs']) * r['params']['NIC Latency (us)'] * 1e-6)
  print "Comm Time(s) %g " % r['cTime']
  print "Total Time(s) %g " % (float(r['time']) + r['cTime'])
  print "Box Size %g" % problemSize 
  print "CNS Version %s" % version 
  print 

def runHopper(filenameNaive, filenameFused):
#Hopper Settings 
  config = {}
  config['NIC BW (GB/s)']=1
  
  print "Results for Hopper"
  runNodeConfig(filenameNaive, config, 64, 2, 'naive')
  runNodeConfig(filenameFused, config, 64, 4, 'fused-simple')
  runNodeConfig(filenameNaive, config, 128, 4, 'naive')
  runNodeConfig(filenameFused, config, 128, 8, 'fused-simple')
  
def runExaNode1(filenameNaive, filenameFused):
#exaNode1 Settings 
  config = {}
  config['GB/s/chip'] = 1024
  config['Gflop/s/core'] = 10 
  config['Cores/chip']= 1024
  config['$/core (MB)']= 128
  config['NIC BW (GB/s)']=100
  config['NIC Latency (us)']=0.4
  
  print "Results for exaNode1"
#blocking factors are for 128KB cache per core 
  runNodeConfig(filenameNaive, config, 64, 6, 'naive')
  runNodeConfig(filenameFused, config, 64, 1, 'fused-simple')
  runNodeConfig(filenameNaive, config, 128, 11, 'naive')
  runNodeConfig(filenameFused, config, 128, 1, 'fused-simple')
  
def runExaNode2(filenameNaive, filenameFused):
#exaNode2 Settings 
  config = {}
  config['GB/s/chip'] = 4096
  config['Gflop/s/core'] = 10 
  config['Cores/chip']= 1024
  config['$/core (MB)']= 128
  config['NIC BW (GB/s)']=400
  config['NIC Latency (us)']=0.02
  
  print "Results for exaNode2"

#blocking factors are for 128KB cache per core 
  runNodeConfig(filenameNaive, config, 64, 6, 'naive')
  runNodeConfig(filenameFused, config, 64, 1, 'fused-simple')
  runNodeConfig(filenameNaive, config, 128, 11, 'naive')
  runNodeConfig(filenameFused, config, 128, 1, 'fused-simple')

def runPIM(filenameNaive, filenameFused):
#PIM Settings 
  config = {}
  config['GB/s/chip'] = 1024
  config['Gflop/s/core'] = 10 
  config['Cores/chip']= 64
  config['$/core (MB)']= 0
  config['NIC BW (GB/s)'] = 400
  config['NIC Latency (us)'] = 0.02
  
  print "Results for PIM"
#blocking factors are for 128KB cache per core
#PIM is working on a 1/16 of the problem 
  runNodeConfig(filenameNaive, config, 64, 1, 'naive')
  runNodeConfig(filenameFused, config, 64, 1, 'fused-simple')
  runNodeConfig(filenameNaive, config, 128, 1, 'naive')
  runNodeConfig(filenameFused, config, 128, 1, 'fused-simple')
  
#program starts here

filenameNaive = '../../examples/xmls/advance-flat.xml'
filenameFused = '../../examples/xmls/advance-fused-nomod.xml'
  
runHopper(filenameNaive, filenameFused)
runExaNode1(filenameNaive, filenameFused)
runExaNode2(filenameNaive, filenameFused)
runPIM(filenameNaive, filenameFused)


#runDemo()
#printBlockCacheTables()

