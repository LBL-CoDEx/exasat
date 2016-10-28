#!/usr/bin/python
from math import sqrt
import sys
sys.path.append('/home/dunat/Projects/exasat/tools/post')

from analyzeXML import StaticAnalysis
import pprint

def runCompute(filename, chemfile, config, problemSize): 

  sa = StaticAnalysis(filename)
  
  r = sa.runModel(config)

  domainSize = problemSize * problemSize * problemSize
  adds = r['adds']/ domainSize
  muls = r['multiplies'] /domainSize
  divs = r['divides'] /domainSize
  specials = r['specials'] /domainSize

#running chemistry 
  config['X Block Size'] = problemSize
  config['Y Block Size'] = problemSize
  config['Z Block Size'] = problemSize

  sa = StaticAnalysis(chemfile)
  r = sa.runModel(config)
  chemadds = r['adds'] / domainSize
  chemmuls = r['multiplies'] /domainSize
  chemdivs = r['divides'] /domainSize
  chemspecials = r['specials'] /domainSize

#  mpoints = (float((problemSize * problemSize * problemSize )) / totaltime /10**6)


#finite diff
  print "Dynamics", problemSize, config['Num Species'],
  print "%g " % adds, 
  print "%g " % muls,
  print "%g " % (divs),
  print "%g " % (specials),
  print 
#chemistry 
  print "Chemistry", problemSize, config['Num Species'],
  print "%g " % chemadds, 
  print "%g " % chemmuls,
  print "%g " % (chemdivs),
  print "%g " % (chemspecials),

  print 

def runExaNode1():
#exaNode1 Settings 

  xmlpath = '/home/dunat/Projects/exasat/examples/xmls/'
  filename = xmlpath + "advance-smc-modified.xml"

  #no simd on div and exp, simd only div, simd both div and exp
  #divcost = [40, 20, 20]
  #specialcost = [125, 125, 42] 
  #version = ['baseline', 'simd-div', 'simd-exp']

  pSize = 128
  species = [9, 21, 53, 71, 107]

  #naive
  blockSize = [16, 16,  16, 8,  8, 8, # for 9 , 21, 53   
               8, 4, 8, 4 ] # for 71 and 107
              # 8, 4,  8, 4,  8, 4,   # 60, 70, 80
              # 8,4, 8,4] # 90, 100

  print "SMC FP instruction mix per grid point"
  print "Version BoxSize NSpecies Add Mul Div Specials"

  dummy = 0 #outside of species loop       
  for nspec in species:

      problemSize = pSize

      config = {}
      config['Num Species'] = nspec
      config['Gflop/s/thread'] = 10 
      config['Threads']= 1024 
      config['GB/s/thread'] = 1
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
       
#program starts here


runExaNode1()

#runDemo()
#printBlockCacheTables()

