#!/usr/bin/python

from analyzeXML import StaticAnalysis

sa = StaticAnalysis('../../examples/xmls/advance.xml')

for config in [{'Problem Size' : 64, 'Blocking Factor': 1, 'GB/s/chip': 1000, 'Gflop/s/core': 10, 'Cores/chip': 1000, '$/core (MB)': 0.0625}, 
               {'Problem Size' : 8, 'Blocking Factor': 2, 'GB/s/chip': 1000, 'Gflop/s/core': 10, 'Cores/chip': 1000, '$/core (MB)': 0.0625}, 
               {'Problem Size' : 8, 'Blocking Factor': 4, 'GB/s/chip': 1000, 'Gflop/s/core': 10, 'Cores/chip': 1000, '$/core (MB)': 0.0625}, 
               {'Problem Size' : 8, 'Blocking Factor': 8, 'GB/s/chip': 1000, 'Gflop/s/core': 10, 'Cores/chip': 1000, '$/core (MB)': 0.0625}, 
               ]:
  r = sa.runModel(config)
  print "Params: %s" % str(r['params'])
  print "Comp: %g Gflops" % (float(r['flops']) / 10**9)
  print "Data: %g GB" % (float(r['BW']) / 2**30)
  print "Time: %g sec" % (r['time'] )
  print
