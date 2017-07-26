import operator
from params import doSymSubs

class options:
  dim = 3
  flag_graph = False
  flag_debug = False
  flag_warn = True
  flag_ignore_gap = True # ignore the gaps in stencil pattern when calculating working set
  warned = set()

def numIters(x):
  return doSymSubs('-(%s) + (%s) + 1' % x)

def diff(x):
  t = type(x)
  return t([x[i] - x[i-1] for i in xrange(1,len(x))])

def mymax(x):
  if len(x) > 0:
    return max(x)
  return 0

def numNonZero(x):
  if type(x) == type(0) or type(x) == type(0.):
    return int(x != 0)
  else:
    return sum(map(numNonZero, x))

def prod(x):
  return reduce(operator.mul, x)

int_types = set(['int'])
float_types = set(['float', 'double'])
