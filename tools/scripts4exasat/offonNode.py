#!/usr/bin/python
from math import sqrt
import sys
sys.path.append('/home/dunat/Projects/exasat/tools/post')

from analyzeXML import StaticAnalysis
import pprint

def runCompute(filename, chemfile, config, problemSize): 

  sa = StaticAnalysis(filename)
  
  r = sa.runModel(config)

  onnodetraffic = ((float(r['BW']))/ 2**30)
  offnodetraffic = ((float(r['cBW']))* 2**(-30))

#running chemistry 
  config['X Block Size'] = problemSize
  config['Y Block Size'] = problemSize
  config['Z Block Size'] = problemSize

  sa = StaticAnalysis(chemfile)
  r = sa.runModel(config)

  chemonnodetraffic = ((float(r['BW']))/ 2**30)

#finite diff
  print "%g " % onnodetraffic, 

#chemistry 
  print "%g " % chemonnodetraffic, 

#total 
  print "%g " % (onnodetraffic + chemonnodetraffic), 

#communication

  print "%g " % offnodetraffic


def runExaNode1():
#exaNode1 Settings 

  xmlpath = '/home/dunat/Projects/exasat/examples/xmls/'
  light = True
  simplefusion = False

  filename = xmlpath + "advance-smc-modified.xml"
  version = ['cache-blocking']

  pSize = 128
  #naive
#noblocking
#  blockSize = [pSize, pSize, pSize, pSize, pSize, pSize, 
#               pSize, pSize, pSize, pSize]
  blockSize = [16, 16,  16, 8,  8, 8, # for 9 , 21, 53   
               8, 4, 8, 4 ] # for 71 and 107

  #no simd on div and exp, simd only div, simd both div and exp
  divcost = [39]
  specialcost = [125] 

  species = [9, 21, 53, 71, 107]

  if simplefusion:
    filename = xmlpath + "advance-smc-fused.xml"#"advance-smc-modified.xml"
    version = ['simple-fusion']
    blockSize = [16, 8, 8, 4, pSize, pSize, pSize, pSize, pSize, pSize]

  if light:
    filename = xmlpath + "advance-smc-fused-light.xml"#"advance-smc-modified.xml"
    version = ['light-fusion']
    blockSize = [16, 16, 16, 8, 8, 4, 8, 4, 8, 2]

  print "Results for exaNode1"
#  print "SMC Version %s" % version 

  print "SMC Problem Size Info Finite Chemistry On-Node Off-Node"
  print "Version BoxSize BlockSize NSpecies",
  print "On-NodeTraffic(GB) ",
  print "On-NodeTraffic(GB) ",
  print "On-NodeTraffic(GB) ",
  print "Off-NodeTraffic(GB)"

  cost = 0 
  for div in divcost:
    p = 0 
    dummy = 0        
    for nspec in species:

      problemSize = pSize
      print version[cost], problemSize, "%dX%d" % (blockSize[dummy], blockSize[dummy + 1]) , 
      print nspec, 
      
      config = {}
      config['Num Species'] = nspec
      config['Gflop/s/thread'] = 10 
      config['Threads']= 1024 
      config['GB/s/thread'] = 1
      config['$/thread group (kB)']= 1024
      config['NIC BW (GB/s)']=100
      config['NIC Latency (us)']=0.4

      config['Division Cost'] = divcost[cost] 
      config['Special Cost'] = specialcost[cost]
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
      
      if nspec == 53 : 
        chemfilename = xmlpath + 'grimech30.xml'

      if nspec == 71 : 
        chemfilename = xmlpath + 'hai.xml'

      if nspec == 107 : 
        chemfilename = xmlpath + 'prf_ethanol.xml'

      runCompute(filename, chemfilename, config, problemSize)

      dummy = dummy + 2
      p = p + 1  

    cost = cost + 1
    print   
#program starts here


runExaNode1()

#runDemo()
#printBlockCacheTables()

