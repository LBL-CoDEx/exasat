#!/usr/bin/python

import sys
import os
import re
from copy import deepcopy

from parser import XMLParser, Collection, Flops, Scalar, Array, ArrayAccess, Conditional
from box import Box
from common import options

# Helper functions

def accessSize(accs):
  """Returns the number of points in a list of (possibly overlapping) accesses using Box library."""
  def makeBox(acc):
    # add 1 since upper bound is inclusive
    return Box(intervals = map(lambda x: (x,x+1) if type(x)==int else (x[0], x[1]+1), acc.index))
  # take the union across boxes to eliminate overlaps
  return reduce(lambda x,y: x.union(y), map(makeBox, accs)).size()

_reported = set()
def reportReuse(linenum, ws, cache):
  if linenum not in _reported:
    _reported.add(linenum)
    print "Reuse report: loop %d WS: %d %s %d" % \
          (linenum, ws, '<=' if ws <= cache else '>', cache)

# helper class for checking conditionals
class TableCondsChecker(object):
  __slots__ = ['table']
  def __init__(self, conds_table):
    if options.flag_ignore_conds:
      print "WARNING: flag_ignore_conds is True, including all conditionally executed blocks!"
    self.table = conds_table
  def check_cond(self, cond):
    if cond.condition not in self.table:
      print cond.condition
      print self.table
      raise Exception("Did not find conditional in table")
    return self.table[cond.condition] == cond.when
  def check_conds(self, conds):
    conds_sat = all(map(lambda x: self.check_cond(x), conds))
    if options.flag_verbose_conditionals and len(conds) > 0:
      print "Conditions", map(lambda x: x.condition, conds), conds_sat
    return conds_sat
  def __call__(self, x):
    # x is either a Conditional or a list of Conditionals
    if options.flag_ignore_conds:
      return True
    if type(x) == Conditional:
      return self.check_cond(x)
    else:
      assert type(x) == list
      return self.check_conds(x)

# Analysis classes

class FlopCount(object):
  slots = ['name', 'adds', 'multiplies', 'divides', 'specials']
  def __init__(self, flops):
    self.name = ''
    self.adds = flops.adds
    self.multiplies = flops.multiplies
    self.divides = flops.divides
    self.specials = flops.specials
  def loop(self, loop, siblings):
    return self * loop.iter_n()
  def __str__(self):
    return "FlopCount A=%s, M=%s, D=%s, S=%s" % \
           (self.adds, self.multiplies, self.divides, self.specials)
  def __iadd__(self, other):
    self.adds += other.adds
    self.multiplies += other.multiplies
    self.divides += other.divides
    self.specials += other.specials
    return self
  def __mul__(self, n):
    result = FlopCount(self)
    result.adds *= n
    result.multiplies *= n
    result.divides *= n
    result.specials *= n
    return result

  @staticmethod
  def visitor(conds_chk):
    # capture condition checker
    def f(arg, conds):
      assert type(arg) == Flops or type(arg) == Scalar or type(arg) == Array
      if type(arg) == Flops and conds_chk(conds):
        return Collection([FlopCount(arg)])
      return Collection()
    return f


class StateVar(object):
  slots = ['name', 'type', 'reads', 'writes']
  def __init__(self, sv=None, array=None, access=None):
    if sv:
      self.name = sv.name
      self.type = sv.type
      self.reads = sv.reads
      self.writes = sv.writes
    elif array and not access:
      # array pointer
      self.name = array.name.split('.')[0]
      self.type = array.type + ' *'
      self.reads = sum(map(lambda x: x.reads + x.writes, array.accesses))
      self.writes = 0
    elif array and access:
      # array element (non-loop varying)
      self.name = array.name + str(access.index)
      self.type = array.type
      self.reads = access.reads
      self.writes = access.writes
    else:
      raise Exception("unsupported constructor call")
  def loop(self, loop, siblings):
    return self * loop.iter_n()
  def __str__(self):
    return "SV %s %s, R=%s, W=%s" % (self.type, self.name, self.reads, self.writes)
  def __iadd__(self, other):
    assert other.name == self.name and other.type == self.type
    self.reads += other.reads
    self.writes += other.writes
    return self
  def __mul__(self, n):
    result = StateVar(self)
    result.reads *= n
    result.writes *= n
    return result

  @staticmethod
  def visitor(conds_chk):
    def f(arg, conds):
      assert type(arg) == Flops or type(arg) == Scalar or type(arg) == Array
      if type(arg) == Scalar and conds_chk(conds):
        return Collection([StateVar(sv=arg)])
      elif type(arg) == Array and conds_chk(conds):
        return Collection([StateVar(array=arg)] + \
                          map(lambda ac: StateVar(array=arg, access=ac),
                              filter(ArrayAccess.isStateVar, arg.accesses)))
      return Collection()
    return f


class ArrayVar(object):
  slots = ['name', 'type', 'reads', 'writes']
  def __init__(self, array):
    if type(array) == ArrayVar:
      # copy constructor
      self.name = array.name
      self.type = array.type
      self.reads = array.reads
      self.writes = array.writes
    elif type(array == Array):
      # from parser.Array object
      self.name = array.name
      self.type = array.type
      self.reads = sum(map(lambda x: 0 if x.isStateVar() else x.reads, array.accesses))
      self.writes = sum(map(lambda x: 0 if x.isStateVar() else x.writes, array.accesses))
    else:
      raise Exception("unknown arg type")
  def loop(self, loop, siblings):
    return self * loop.iter_n()
  def __str__(self):
    return "AR %s %s, L=%s, S=%s" % (self.type, self.name, self.reads, self.writes)
  def __iadd__(self, other):
    assert other.name == self.name and other.type == self.type
    self.reads += other.reads
    self.writes += other.writes
    return self
  def __mul__(self, n):
    result = ArrayVar(self)
    result.reads *= n
    result.writes *= n
    return result

  @staticmethod
  def visitor(conds_chk):
    def f(arg, conds):
      assert type(arg) == Flops or type(arg) == Scalar or type(arg) == Array
      if type(arg) == Array and not arg.onlyStateVars() and conds_chk(conds):
        return Collection([ArrayVar(arg)])
      return Collection()
    return f


class WorkingSet(object):
  """The working set associated with an array in a code region."""
  slots = ['name', 'type', 'accesses', 'size_']
  def __init__(self, array, clearAccs = False, params = None):
    self.name = array.name
    self.type = array.type
    self.accesses = []
    if not clearAccs:
      self.accesses = deepcopy(filter(lambda x: not x.isStateVar(), array.accesses))
      if params:
        self.accesses = map(lambda x: x.subParams(params), self.accesses)
    self.size_ = None
  def __str__(self):
    s = "WS %s %s\n" % (self.type, self.name)
    for ac in self.accesses:
      s += str(ac) + '\n'
    return s
  def __iadd__(self, other):
    assert other.name == self.name and other.type == self.type
    map(self.consume, other.accesses)
    return self
  def size(self):
    # memoize
    if not self.size_:
      self.size_ = accessSize(self.accesses)
    return self.size_
  def bytes(self):
    return self.size() * bytesPerWord
  def consume(self, acc):
    # access regions can overlap
    self.accesses.append(acc)
    self.size_ = None
  def loop(self, loop, siblings = None):
    """Unroll along loop variable, combining accesses along loop dimension."""

    def tupleDel(t, i):
      return tuple(t[:i] + t[i+1:])
    def tupleIns(t, i, v):
      return tuple(t[:i] + (v,) + t[i:])
    def appendToDict(d, k, v):
      if k not in d:
        d[k] = []
      d[k].append(v)

    result = WorkingSet(self, clearAccs = True)
    buckets = {}
    for acc in self.accesses:
      if loop.loopvar not in acc.loopvars:
        result.consume(acc) # nothing to unroll
      else:
        # sort accesses into buckets to be merged
        i = acc.loopvars.index(loop.loopvar) # which index corresponds to the loopvar
        key = (i, tupleDel(acc.index, i), tupleDel(acc.loopvars, i)) # remove the index corresponding to the loopvar
        appendToDict(buckets, key, acc.index[i]) # add the index to the key's bucket

    # each bucket now contains a set of accesses that differ only in the dimension of the loop
    # these accesses can be merged, conservatively, by taking the minimum and maximum offset
    for ((i, idx, lvars), offsets) in buckets.iteritems():
      # insert an access interval covering the looped over accesses
      (lb, ub) = loop.range
      r = (lb + min(offsets), ub + max(offsets))
      idx = tupleIns(idx, i, r)
      lvars = tupleIns(lvars, i, '')
      result.consume(ArrayAccess(index=idx, loopvars=lvars))

    return result

  @staticmethod
  def visitor(conds_chk):
    def f(arg, conds):
      assert type(arg) == Flops or type(arg) == Scalar or type(arg) == Array
      if type(arg) == Array and not arg.onlyStateVars() and conds_chk(conds):
        return Collection([WorkingSet(arg)])
      return Collection()
    return f


class TrafficRegion(object):
  """The memory traffic associated with a list of regions and a count for how
     many times those regions are accessed."""
  slots = ['accesses', 'count', 'size_']
  def __init__(self, accesses, count = 1):
    self.accesses = accesses
    self.count = count
    self.size_ = self.count * accessSize(self.accesses)
  def __str__(self):
    s = " TrafficRegion, count=%g:\n" % self.count
    for acc in self.accesses:
      s += str(acc) + '\n'
    return s
  def __imul__(self, n):
    self.count *= n
    self.size_ *= n
    return self
  def size(self):
    return self.size_


class Traffic(object):
  """The memory traffic required to execute a code region.
    
     Contains several TrafficRegions and the working set to determine reuse in loops."""
  slots = ['name', 'type', 'regions', 'ws', 'ws_block_n', 'params', 'blockParams', 'cache']
  def __init__(self, array=None, params=None, blockParams=None, cache=None, copy=None):
    if copy:
      # make a copy (excluding traffic regions and ws)
      self.name = copy.name
      self.type = copy.type
      self.params = copy.params           # for problem size and others
      self.blockParams = copy.blockParams # for cache blocking only
      self.cache = copy.cache             # size of cache in machine model
    else:
      # copy from an Array object
      self.name = array.name
      self.type = array.type

      # copy non-state var accesses, do parameter substitution
      accesses = filter(lambda x: not x.isStateVar(), array.accesses)
      accesses = map(lambda x: x.subParams(params), accesses)

      # traffic regions accessed
      self.regions = [TrafficRegion(accesses)]

      # need to track working set to compute data reuse
      self.ws = WorkingSet(array, params=params)
      self.ws_block_n = 1 # single iteration

      # other parameters
      self.params = params
      self.blockParams = blockParams
      self.cache = cache
  def __str__(self):
    s = "MT %s %s, words=%d, bytes=%d\n" % (self.type, self.name, self.size(), self.bytes())
    for region in self.regions:
      s += str(region)
    return s
  def __iadd__(self, other):
    assert other.name == self.name and other.type == self.type
    self.ws += other.ws # add other's working set to ours
    map(self.consume, other.regions)
    return self
  def __imul__(self, n):
    for region in self.regions:
      region *= n
    return self
  def consume(self, region):
    # these regions may overlap
    self.regions.append(region)
  def bytes(self):
    return self.size() * bytesPerWord
  def size(self):
    return sum(map(lambda x: x.size(), self.regions))
  def loop(self, origLoop, siblings):
    """If enough cache for reuse within a block, traffic equals the blocked working set times the number of blocks.

       Otherwise, multiply traffic per iteration by number of iterations."""

    # total working set for all arrays touched in the loop (across all siblings) for one iteration
    # siblings are all different arrays, so assume no aliasing (overlap) between them
    loopWS = sum(map(lambda x: x.ws.bytes(), siblings))

     # only need to substitute params for the loop bounds
    loop = origLoop.subParams(self.params, shallow=True)
    blockLoop = origLoop.subParams(self.blockParams, shallow=True)

    # number of blocks in the current loop dimension
    numBlocks = float(loop.range[1]      - loop.range[0]+1) / \
                (blockLoop.range[1] - blockLoop.range[0]+1)

    result = Traffic(copy = self)
    result.ws = self.ws.loop(blockLoop) # working set of the blocked loop
    result.ws_block_n = self.ws_block_n * numBlocks # number of blocks total

    reportReuse(loop.linenum, loopWS, self.cache)
    if loopWS <= self.cache:
      # WS of all siblings fit in cache
      # traffic equals blocked working set times number of total blocks
      # assumes no reuse between blocks (i.e. problem blocked correctly for cache)
      result.regions = [TrafficRegion(result.ws.accesses, count=result.ws_block_n)]
    else:
      # no reuse across iterations, multiply traffic by iteration count
      result.regions = self.regions
      result *= loop.iter_n()
    return result

  @staticmethod
  def visitor(conds_chk, params, blockParams, cache):
    # capture problem size, blocking params, and machine model cache
    def f(arg, conds):
      assert type(arg) == Flops or type(arg) == Scalar or type(arg) == Array
      if type(arg) == Array and not arg.onlyStateVars() and conds_chk(conds):
        return Collection([Traffic(arg, params, blockParams, cache)])
      return Collection()
    return f

class StaticAnalysis(object):

  slots = ['functions']

  def __init__(self, filename):
    self.functions = XMLParser(filename).functions

  def dump(self, params, blockParams, cache, conds_chk):
    for function in self.functions:
      print "************"
      print "* %s *" % function.name
      print "************"
      for loop in function.body.loops:
        print
        print "*************"
        print "* Loop %4d *" % loop.linenum
        print "*************"
        print

        if flagSubParams:
          loop = loop.subParams(params)

        print "Floating Point Ops (A/S/M/D):"
        print loop.visit(FlopCount.visitor(conds_chk))

        print "State Variables (R/W):"
        print loop.visit(StateVar.visitor(conds_chk))

        print "Array Variables (L/S):"
        print loop.visit(ArrayVar.visitor(conds_chk))

        print "Working Set:"
        print loop.visit(WorkingSet.visitor(conds_chk))

        mt = loop.visit(Traffic.visitor(conds_chk, params, blockParams, cache))
        print
        print "Memory Traffic (L/S) using cache model: %d " % sum(map(lambda x: x.size(), mt))
        print mt

# some test parameters

flagSubParams = False

params = [('lo(1)', '1'), \
          ('lo(2)', '1'), \
          ('lo(3)', '1'), \
          ('hi(1)', '128'), \
          ('hi(2)', '128'), \
          ('hi(3)', '128'), \
          ('fine_lo(1)', '1'), \
          ('fine_lo(2)', '1'), \
          ('fine_lo(3)', '1'), \
          ('fine_hi(1)', '128'), \
          ('fine_hi(2)', '128'), \
          ('fine_hi(3)', '128'), \
          ('crse_lo(1)', '1'), \
          ('crse_lo(2)', '1'), \
          ('crse_lo(3)', '1'), \
          ('crse_hi(1)', '64'), \
          ('crse_hi(2)', '64'), \
          ('crse_hi(3)', '64'), \
          ('ilo1', '1'), \
          ('ilo2', '1'), \
          ('ilo3', '1'), \
          ('ihi1', '128'), \
          ('ihi2', '128'), \
          ('ihi3', '128'), \
          ('__BoxX__', '128'), \
          ('__BoxY__', '128'), \
          ('__BoxZ__', '128'), \
          ('ng' , '4'), \
          ('nspecies', '9'), \
         ]

blockParams = [('lo(1)', '1'), \
               ('lo(2)', '1'), \
               ('lo(3)', '1'), \
               ('hi(1)', '128'), \
               ('hi(2)', '128'), \
               ('hi(3)', '128'), \
               ('fine_lo(1)', '1'), \
               ('fine_lo(2)', '1'), \
               ('fine_lo(3)', '1'), \
               ('fine_hi(1)', '128'), \
               ('fine_hi(2)', '128'), \
               ('fine_hi(3)', '128'), \
               ('crse_lo(1)', '1'), \
               ('crse_lo(2)', '1'), \
               ('crse_lo(3)', '1'), \
               ('crse_hi(1)', '64'), \
               ('crse_hi(2)', '64'), \
               ('crse_hi(3)', '64'), \
               ('ilo1', '1'), \
               ('ilo2', '1'), \
               ('ilo3', '1'), \
               ('ihi1', '128'), \
               ('ihi2', '128'), \
               ('ihi3', '128'), \
               ('__BoxX__', '128'), \
               ('__BoxY__', '128'), \
               ('__BoxZ__', '128'), \
               ('ng' , '4'), \
               ('nspecies', '9'), \
              ]

conds_table = {'present(courno)': True}

bytesPerWord = 8
cache = 1*(2**15)

def main(args):
  if len(args) < 2:
    print "Usage: %s <xml-file>" % args[0]
    sys.exit(1)
  xmlfile = args[1]

  sa = StaticAnalysis(xmlfile)
  sa.dump(params, blockParams, cache, TableCondsChecker(conds_table))
 
if __name__ == '__main__':
  main(sys.argv)
