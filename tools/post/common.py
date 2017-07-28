#!/usr/bin/env python

""" Some common ExaSAT helper functions. """
__author__ = "Cy Chan"
__copyright__ = "Copyright 2016, The Regents of the University of California, through Lawrence Berkeley National Laboratory"
__credits__ = ["Cy Chan"]
__license__ = "Modified BSD License (see LICENSE file)"
__version__ = "2.0"
__maintainer__ = "Cy Chan"
__email__ = "cychan@lbl.gov"
__status__ = "Production"

import operator

class options:
  flag_debug = True
  flag_warn = True
  flag_ignore_gap = True # ignore the gaps in stencil pattern when calculating working set
  flag_ignore_conds = False # conservatively assume *all* code within *all* conditional blocks are executed
  flag_verbose_conditionals = True
  flag_verbose_parser = False
  warned = set()

def numIters(x):
  # TODO: this is for inclusive upper bound, make sure consistent across C and Fortran
  return -x[0] + x[1] + 1

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
