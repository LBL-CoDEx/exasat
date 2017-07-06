import xml.dom.minidom
import ast
import operator
import copy
from box import Box

from common import options, collapse, numIters, rangeMerge, rangeDisjoint
from params import doSymSubs

def parseExpr(s):
  try:
    return int(s) # try to do the fast thing first
  except ValueError as e:
    return doSymSubs(s)

def dprint(s):
  if options.flag_debug:
    print s

def getChildren(node, tag):
  return filter(lambda x: x.nodeName == tag, node.childNodes)

def arrayName(name, component):
  if component != '':
    component = doSymSubs(component)
    return '%s.%s' % (name, component)
  else:
    return name

def makeName(node):
  return arrayName(str(node.getAttribute('name')), \
                   str(node.getAttribute('component')))

def parseTuple(s, f = lambda x: x):
  """Parses a string containing a tuple.
    
     Returns a tuple of elements with f applied to each element.
     Returns a tuple of strings by default."""
  s = s.strip('() ').split(',') # remove enclosing parens and whitespace and split on commas
  return tuple(map(lambda x: f(x.strip()), s)) # strip whitespace and cast to tuple of ints

def subToInt(x, params):
  if type(x) != int:
    return int(doSymSubs(x, params))
  return x


class Collection(object):
  """Generic collection class for use with visitor.
    
     Items must have a name member and implement __iadd__ (or __add__), and loop.
     Items need __str__ if collection is printed."""
  slots = ['d']
  def __init__(self, items = None, colls = None):
    '''Accepts list of items or list of collections of items.'''
    self.d = {}
    if items != None:
      map(self.consume, items)
    elif colls != None:
      map(self.merge, colls)
  def consume(self, item):
    if item.name in self.d:
      self.d[item.name] += item
    else:
      self.d[item.name] = item
  def merge(self, other):
    assert type(other) == Collection
    map(self.consume, other.d.itervalues())
  def loop(self, loop):
    return Collection(map(lambda x: x.loop(loop, self), self))
  def __str__(self):
    s = ""
    for item in sorted(self.d.values(), key=lambda x: x.name):
      s += str(item) + '\n'
    return s
  def __iter__(self):
    return self.d.itervalues()
  def __add__(self, other):
    result = Collection(self)
    result.merge(other)
    return result


class Flops(object):
  slots = ['adds', 'multiplies', 'divides', 'specials']
  def __init__(self, node):
    for x in self.slots:
      self.__setattr__(x, int(node.getAttribute(x)))
    dprint(self)
  def __str__(self):
    return "Flops(%d, %d, %d, %d)" % \
           tuple(map(lambda x: self.__getattribute__(x), self.slots))


class Scalar(object):
  slots = ['name', 'type', 'const', 'reads', 'writes']
  def __init__(self, node):
    self.name = str(node.getAttribute('name'))
    self.type = str(node.getAttribute('datatype'))
    self.const = str(node.getAttribute('isConstant'))
    self.reads = int(node.getAttribute('reads'))
    self.writes = int(node.getAttribute('writes'))
    dprint(self)
  def __str__(self):
    return "%s %s, R=%d, W=%d" % (self.type, self.name, self.reads, self.writes)


class ArrayAccess(object):
  slots = ['index', 'loopvars', 'reads', 'writes']
  def __init__(self, node=None, index=None, loopvars=None):
    if node:
      self.index = parseTuple(node.getAttribute('offset'), parseExpr)
      self.loopvars = parseTuple(node.getAttribute('dependentloopvar'), str)
      self.reads = int(node.getAttribute('reads' ))
      self.writes = int(node.getAttribute('writes'))
      dprint(self)
    else:
      self.index = index
      self.loopvars = loopvars
      self.reads = 0
      self.writes = 0
  def __str__(self):
    idxStr = '('+','.join(map(str, self.index))+')'
    lvStr = '('+','.join(self.loopvars)+')'
    return "  %s+%s, R=%d, W=%d" % \
           (idxStr, lvStr, self.reads, self.writes)
  def __sub__(self, other):
    assert self.loopvars == other.loopvars
    boxDiff = Box(intervals = self.index) - Box(intervals = other.index)
    return map(lambda x: x.intervals, boxDiff.contents)
  def isStateVar(self):
    return all(map(lambda x: x == '', self.loopvars))
  def subParams(self, params):
    result = ArrayAccess()
    result.index = tuple(map(lambda x: subToInt(x, params), self.index))
    result.loopvars = self.loopvars
    result.reads = self.reads
    result.writes = self.writes
    return result


class Array(object):
  slots = ['name', 'type', 'accesses']
  def __init__(self, node = None):
    if node:
      self.name = makeName(node)
      self.type = str(node.getAttribute('datatype'))
      dprint(self)
      self.accesses = map(ArrayAccess, getChildren(node, 'access'))
  def __str__(self):
    return "%s %s" % (self.type, self.name)
  def onlyStateVars(self):
    return all(map(ArrayAccess.isStateVar, self.accesses))
  def subParams(self, params):
    result = Array()
    result.name = self.name
    result.type = self.type
    result.accesses = map(lambda x: x.subParams(params), self.accesses)
    return result


class CodeBlock(object):
  """Branchless section of code with conditions for execution.

     Contains arithmetic, scalar and array accesses, and communication."""
  slots = ['flops', 'scalars', 'arrays', 'conds']
  def __init__(self, node = None, conds = []):
    if node:
      # XXX: gracefully handle code blocks missing attributes
      if node.getAttribute('adds'):
        self.flops = Flops(node)
      else:
        print "WARNING: No Flops attributes on code block!"
      self.conds = conds
      dprint(self)
      self.scalars = map(Scalar, getChildren(node, 'scalar'))
      self.arrays = map(Array, getChildren(node, 'array'))
  def __str__(self):
    return "Code Block (%s):" % str(map(str, self.conds))
  def visit(self, f):
    return Collection(colls = [f(self.flops, self.conds)] + \
                              map(lambda s: f(s, self.conds), self.scalars) + \
                              map(lambda a: f(a, self.conds), self.arrays))
  def subParams(self, params):
    result = CodeBlock()
    result.flops = self.flops
    result.scalars = self.scalars
    result.arrays = map(lambda x: x.subParams(params), self.arrays)
    result.conds = self.conds
    return result


class Conditional(object):
  """Encapsulates a conditional."""
  slots = ['linenum', 'condition', 'when']
  def __init__(self, node):
    assert node.nodeName == "if" or node.nodeName == "else"
    tag = 'linenum' if node.nodeName == 'if' else 'iflinenum'
    self.linenum = int(node.getAttribute(tag))
    self.condition = str(node.getAttribute('conditional'))
    self.when = (node.nodeName == "if")
  def __str__(self):
    return "Conditional: " + str((self.linenum, self.condition, self.when))


class Body(object):
  """A function or loop body: contains information on enclosed code blocks and loops."""
  slots = ['codeblocks', 'loops']
  def __init__(self, node = None, conds = []):
    """Does recursive traversal of conditional blocks within body."""
    def traverse(node, tag, conds):
      """Traverse nested if/else blocks in the body of a function or loop."""
      if tag:
        result = map(lambda x: (x, conds), getChildren(node, tag))
      else:
        result = [(node, conds)]
      for cnode in getChildren(node, 'if') + getChildren(node, 'else'):
        c2 = conds + [Conditional(cnode)]
        result.extend(traverse(cnode, tag, c2))
      return result

    if node:
      self.codeblocks = map(lambda x: CodeBlock(*x), traverse(node, None, conds))
      self.loops = map(lambda x: Loop(*x), traverse(node, 'loop', conds))

  def visit(self, f):
    return Collection(colls = map(lambda x: x.visit(f), self.codeblocks) + \
                              map(lambda x: x.visit(f), self.loops))
  def subParams(self, params):
    result = Body()
    result.codeblocks = map(lambda x: x.subParams(params), self.codeblocks)
    result.loops = map(lambda x: x.subParams(params), self.loops)
    return result


class Loop(object):
  """Contains information on loop variables, bounds, strides, and loop body."""
  slots = ['name', 'loopvar', 'linenum', 'range', 'stride', 'conds', 'body']
  def __init__(self, node = None, conds = []):
    if node:
      self.loopvar = str(node.getAttribute('loopvar'))
      self.linenum = int(node.getAttribute('linenum'))
      self.range = (parseExpr(str(node.getAttribute('lowerbound'))),
                    parseExpr(str(node.getAttribute('upperbound'))))
      self.stride = int(node.getAttribute('stride'))
      self.conds = conds
      dprint(self)
      self.body = Body(node, conds)
  def __str__(self):
    return "Loop %d: %s = [%s, %s] / %d, %s" % \
           (self.linenum, self.loopvar,
            self.range[0], self.range[1], self.stride, str(map(str, self.conds)))
  def iter_n(self):
    return numIters(self.range) / self.stride
  def visit(self, f):
    return self.body.visit(f).loop(self)
  def subParams(self, params, shallow = False):
    result = Loop()
    result.loopvar = self.loopvar
    result.linenum = self.linenum
    result.range = tuple(map(lambda x: subToInt(x, params), self.range))
    result.stride = subToInt(self.stride, params)
    result.conds = self.conds
    result.body = self.body
    if not shallow:
      result.body = self.body.subParams(params)
    return result


class Function(object):
  """Contains information on passed parameters and function body."""
  slots = ['name', 'params', 'local', 'body']
  def __init__(self, node = None):
    if node:
      self.name = str(node.getAttribute('name'))
      dprint("Parsing function %s ..." % self.name)
      self.params = map(makeName, getChildren(node, 'nonlocal'))
      self.local = map(makeName, getChildren(node, 'local'))
      self.body = Body(node)
  def visit(self, f):
    return self.body.visit(f)


class XMLParser(object):
  """Contains information on enclosed modules."""
  slots = ['doc', 'functions']
  def __init__(self, filename):
    assert type(filename) == type('')
    self.doc = xml.dom.minidom.parse(filename)
    program = getChildren(self.doc, 'program')[0]
    self.functions = map(Function, getChildren(program, 'function'))

class MachineXMLParser(object):
  __slots__  = ['doc', 'params']
  def __init__(self, filename):
    assert type(filename) == type('')
    self.doc = xml.dom.minidom.parse(filename)
    self.parse_params()
  def parse_params(self):
    self.params = {}
    machine = self.doc.firstChild
    for prop in getChildren(machine, 'prop'):
      key = str(prop.getAttribute('key'))
      value = float(prop.getAttribute('val'))
      self.params[key] = value
