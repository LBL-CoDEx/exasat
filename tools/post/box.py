#!/usr/bin/env python

""" Helper class for box algebra. """
__author__ = "Cy Chan"
__copyright__ = "Copyright 2016, The Regents of the University of California, through Lawrence Berkeley National Laboratory"
__credits__ = ["Cy Chan"]
__license__ = "Modified BSD License (see LICENSE file)"
__version__ = "2.0"
__maintainer__ = "Cy Chan"
__email__ = "cychan@lbl.gov"
__status__ = "Production"

import operator
import math

inf = float("inf")

def prod(a):
  return reduce(operator.mul, a)

def collapse(a):
  return reduce(operator.add, a, BoxSet([]))

def diff(a):
  return a[1] - a[0]

def intervalIntersect(a, b):
  return [max(a[0], b[0]), min(a[1], b[1])]

def intervalIsEmpty(a):
  return a[1] <= a[0]


# Represents a box
class Box(object):

  __slots__ = ['dim', 'start', 'end', 'intervals', 'steps', 'boxid', 'locid', 'links', 'level', 'inLastLevel', 'node']

  def __init__(self, start=None, end=None, intervals=None, steps=None, boxid=None, level=None, locid=0):
    self.dim = 0
    if start and end:
      self.dim = len(start)
      if len(end) != self.dim:
        raise Exception('dimension mismatch')
      self.start = tuple(start)
      self.end = tuple(end)
      self.intervals = tuple(zip(start, end))
    elif intervals:
      self.dim = len(intervals)
      self.intervals = tuple(intervals)
      (self.start, self.end) = zip(*intervals)
    self.steps = steps
    if not steps:
      self.steps = (1,) * self.dim
    self.boxid = boxid
    self.level = level
    self.locid = locid

  def __str__(self):
    result = 'Box(%s, %s, NP=%d' % (self.start, self.end, self.numPoints())
    if self.boxid != None:
      result += ', ID=%d' % self.boxid
    if self.locid != None:
      result += ', LOC=%d' % self.locid
    if self.level != None:
      result += ', LVL=%d' % self.level
    result += ')'
    return result

  def __repr__(self):
    return self.__str__()

  def isEmpty(self):
    return any(map(intervalIsEmpty, self.intervals))

  def sizes(self):
    return map(diff, self.intervals)

  def size(self):
    return prod(self.sizes())

  def setSteps(self, steps):
    if len(steps) != self.dim:
      raise Exception('dimension mismatch')
    self.steps = steps

  def numPoints(self):
    return self.size() / prod(self.steps)

  def weight(self, refRatios):
    return self.numPoints() * prod(refRatios[:self.level+1])

  def numFlux(self):
    np = self.numPoints()
    sizes = map(lambda x,y: diff(x) / y, self.intervals, self.steps)
    return 2 * sum(map(lambda x: np / x, sizes))

  def __nonzero__(self):
    return not self.isEmpty()

  # ~A = spatial inverse of A
  def __invert__(self):
    result = []
    for i in xrange(0, self.dim):
      i1 = ((-inf, inf),) * (self.dim - (i+1)) + ((-inf, self.intervals[-(i+1)][0]),) + self.intervals[-i:]
      i2 = ((-inf, inf),) * (self.dim - (i+1)) + ((self.intervals[-(i+1)][1], inf),) + self.intervals[-i:]
      result += [Box(intervals=i1[:self.dim], steps=self.steps), Box(intervals=i2[:self.dim], steps=self.steps)]
    return BoxSet(result)

  def translate(self, shift):
    shiftedIntervals = tuple(map(lambda x, y: map(lambda z: z + y, x),
                                 self.intervals, shift))
    return Box(intervals = shiftedIntervals)

  def intersect(self, b):
    if isinstance(b, Box):
      return Box(intervals=map(intervalIntersect, self.intervals, b.intervals), steps=self.steps)
    elif isinstance(b, BoxSet):
      return BoxSet(map(self.intersect, b.contents)).nonempty()
    else:
      raise Exception('unexpected type')

  def union(self, b):
    return BoxSet(self).union(b)

  # A \ B = A \intersect ~B
  def __sub__(self, b):
    if isinstance(b, Box):
      i = self.intersect(b) 
      if i.isEmpty():
        result = BoxSet(self)
      else:
        result = self.intersect(~b)
      return result
    elif isinstance(b, BoxSet):
      result = BoxSet(self)
      for curBox in b.contents:
        result = collapse(map(lambda x: x - curBox, result.contents))
      return result
    else:
      raise Exception('unexpected type')

  # extend faces by d
  def extendFaces(self, d):
    result = BoxSet(self)
    for i in xrange(0, self.dim):
      x = list(self.intervals)
      x[i] = (self.intervals[i][0] - d[i], self.intervals[i][0])
      result += Box(intervals=list(x))
      x[i] = (self.intervals[i][1], self.intervals[i][1] + d[i])
      result += Box(intervals=list(x))
    return result

  # extend border by d
  def extend(self, d):
    return Box(intervals=map(lambda x, y: (x[0] - y, x[1] + y), self.intervals, d), steps=self.steps)

  # ghost region around box
  def ghost(self, d):
    return self.extend(d) - self

  def cover(self, steps):
    s = map(lambda x,y: math.floor(float(x) / y) * y, self.start, steps)
    e = map(lambda x,y: math.ceil(float(x) / y) * y, self.end, steps)
    return Box(s, e, steps = steps)

  # cut box with plane along dimension dim at val
  def cut(self, dim, val):
    if self.start[dim] < val and self.end[dim] > val:
      intervals = list(self.intervals)
      intervals[dim] = (self.start[dim], val)
      box1 = Box(intervals = intervals)
      intervals[dim] = (val, self.end[dim])
      box2 = Box(intervals = intervals)
      return BoxSet([box1, box2])
    else:
      return BoxSet(self)

  def shadow(self):
    return Shadow(self)

  def getLinks(self, tag):
    return sorted(self.links[tag].items(), key=lambda x: x[0].boxid)


# Represents a set of boxes
class BoxSet(object):

  __slots__ = ['contents', 'locid']

  def __init__(self, x):
    if type(x) == type([]):
      self.contents = x
    elif type(x) == type(Box()):
      self.contents = [Box(intervals=x.intervals, steps=x.steps, locid=x.locid)]
    else:
      raise Exception("BoxSet() called with unsupported argument: (%s) %s" %
                      str(type(x)), str(x))

  def __str__(self):
    result = ''
    result += 'Box Set (count=%d, NP=%d' % (len(self.contents), self.numPoints())
    if hasattr(self, 'locid') and self.locid != None:
      result += ', LOC=%d' % self.locid
    result += '):\n'
    for box in self.contents:
      result += '  ' + str(box) + '\n'
    return result

  def __add__(self, other):
    if isinstance(other, BoxSet):
      return BoxSet(self.contents + other.contents)
    elif isinstance(other, Box):
      return BoxSet(self.contents + [other])
    else:
      raise Exception('unexpected type')

  def __len__(self):
    return len(self.contents)

  def setSteps(self, steps):
    for box in self.contents:
      box.steps = steps

  def nonempty(self):
    self.contents = filter(lambda x: not x.isEmpty(), self.contents)
    return self

  def size(self):
    return sum(map(lambda x: x.size(), self.contents))

  def numPoints(self):
    return sum(map(lambda x: x.numPoints(), self.contents))

  def numFlux(self):
    if len(self.contents) == 1:
      return self.contents[0].numFlux()
    else:
      # flux for multiple boxes needs to subtract touching faces
      raise Exception('flux for multiple boxes not implemented yet')

  def isEmpty(self):
    return self.size() == 0

  def __nonzero__(self):
    return self.size() > 0

  def __sub__(self, b):
    return collapse(map(lambda x: (x - b), self.contents))

  def intersect(self, b):
    if isinstance(b, Box):
      b = BoxSet([b])
    return collapse(map(lambda x: x.intersect(b), self.contents))

  # A U B = A + B \ A
  def union(self, b):
    return self + (b - self)

  # extend borders by d (reduce by union so boxes are non-overlapping)
  def extend(self, d):
    temp = map(lambda x: x.extend(d), self.contents)
    return reduce(lambda x,y: x.union(y), temp, BoxSet([]))

  # extend faces by d (reduce by union so boxes are non-overlapping)
  def extendFaces(self, d):
    temp = map(lambda x: x.extendFaces(d), self.contents)
    return reduce(lambda x,y: x.union(y), temp, BoxSet([]))

  # ghost region around boxes
  def ghost(self, d):
    return self.extend(d) - self

  def cover(self, steps):
    temp = map(lambda x: x.cover(steps), self.contents)
    return reduce(lambda x,y: x.union(y), temp, BoxSet([]))

  # cut all boxes with plane along dimension dim at val
  def cut(self, dim, val):
    return collapse(map(lambda x:x.cut(dim, val), self.contents))

  def shadow(self):
    return Shadow(self)


class Shadow(object):

  __slots__ = ['np', 'locid']

  def __init__(self, b):
    self.np = b.numPoints()
    self.locid = b.locid

  def numPoints(self):
    return self.np

  def __nonzero__(self):
    return self.np > 0


def unitTests2D():

  a = Box([ 0, 0], [3, 3])
  tests = [Box([ 1,  1], [2, 2]), # interior
           Box([ 1,  1], [4, 2]), # one side
           Box([ 1, -1], [2, 4]), # two sides (opp)
           Box([ 2,  2], [4, 4]), # two sides (adj)
           Box([-1,  2], [4, 4]), # three sides
           Box([-1, -1], [4, 4])] # exterior

  print 'a: ' + str(a)
  for (i, t) in enumerate(tests):
    print 't%d: ' % i + str(t)
    print 'a ^ t%d: ' % i + str(a.intersect(t))
    print 'a \ t%d: ' % i + str(a - t)
    print 'a U t%d: ' % i + str(a.union(t))
  print

  print 'a \ (t2 U t1): ' + str(a - tests[2].union(tests[1]))

  # ghost region around ((0, 0), (10, 20))
  set1 = BoxSet([Box([-1, -1], [11,  0]), \
                 Box([-1, 20], [11, 21]), \
                 Box([-1,  0], [ 0, 20]), \
                 Box([10,  0], [11, 20])])

  # ghost region around ((10, 0), (20, 10))
  set2 = BoxSet([Box([ 9, -1], [21,  0]), \
                 Box([ 9, 10], [21, 11]), \
                 Box([ 9,  0], [10, 10]), \
                 Box([20,  0], [21, 10])])

  print 'set1 ^ set2: ' + str(set1.intersect(set2))
  print 'set1 \ set2: ' + str(set1 - set2)
  print 'set1 U set2: ' + str(set1.union(set2))

def unitTests3D():

  a = Box([0, 0, 0], [3, 3, 3])
  tests = [Box([ 1,  1,  1], [2, 2, 2]), # interior
           Box([ 1,  1,  1], [4, 2, 2]), # one face
           Box([ 1,  1, -1], [2, 2, 4]), # two faces (opp)
           Box([ 2,  2,  1], [4, 4, 2]), # two faces (adj)
           Box([ 1, -1,  1], [4, 4, 2]), # three faces (opp)
           Box([ 2,  2,  2], [4, 4, 4]), # three faces (adj)
           Box([-1,  2,  2], [4, 4, 4]), # four faces (~adj)
           Box([-1, -1,  1], [4, 4, 2]), # four faces (~opp)
           Box([-1, -1,  2], [4, 4, 4]), # five faces
           Box([-1, -1, -1], [4, 4, 4])] # exterior

  print 'a: ' + str(a)
  for (i, t) in enumerate(tests):
    print 't%d: ' % i + str(t)
    print 'a ^ t%d: ' % i + str(a.intersect(t))
    print 'a \ t%d: ' % i + str(a - t)
    print 'a U t%d: ' % i + str(a.union(t))
  print

  print 'a \ (t2 U t1): ' + str(a - tests[2].union(tests[1])) 


def unitTests():
  unitTests2D()
  unitTests3D()


if __name__ == "__main__":
  unitTests()
