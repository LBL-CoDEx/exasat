import operator
from params import doSymSubs

class options:
  dim = 3
  flag_debug = True
  flag_warn = True
  flag_ignore_gap = True # ignore the gaps in stencil pattern when calculating working set
  warned = set()

def numIters(x):
  # NOTE: this is for inclusive upper bound, make sure consistent across C and Fortran
  if type(x[0]) == int and type(x[1]) == int:
    return -x[0] + x[1] + 1
  return doSymSubs('-(%s) + (%s) + 1' % x)

def diff(x):
  t = type(x)
  return t([x[i] - x[i-1] for i in xrange(1,len(x))])

def mymax(x):
  if len(x) > 0:
    return max(x)
  return 0

def numNonZero(x):
  if type(x) == int or type(x) == float:
    return int(x != 0)
  else:
    return sum(map(numNonZero, x))

def prod(x):
  return reduce(operator.mul, x)

def collapse(x):
  return reduce(operator.add, x, [])

# internal helpers

def rstd(i):
  if type(i) != tuple or len(i) == 1:
    return (i, i)
  return i

def rdisj(i1, i2):
  return i1[1] < i2[0] or i1[0] > i2[1]

def rmerge(i1, i2):
  return (min(i1[0], i2[0]), max(i1[1], i2[1]))

# externally facing

def rangeDisjoint(i1, i2):
  return rdisj(rstd(i1), rstd(i2))

def rangeMerge(i1, i2):
  return rmerge(rstd(i1), rstd(i2))
