import xml.dom.minidom
import ast
import operator
import copy
from box import Box

from common import options, collapse, numIters, rangeMerge, rangeDisjoint
from params import doSymSubs, doNameSubs
from collection import Collection

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
    component = doNameSubs(component)
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
  if type(x) == int:
    result = x
  else:
    try:
      result = int(doSymSubs(x, params))
    except Exception as e:
      print "Could not do integer parameter substitution for expression: ", x
      print "Parameters available: ", params
      raise e
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
  def collect(self, f):
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
    assert node.nodeName == 'if' or node.nodeName == 'else'
    tag = 'linenum' if node.nodeName == 'if' else 'iflinenum'
    self.linenum = int(node.getAttribute(tag))
    self.condition = str(node.getAttribute('conditional'))
    self.when = (node.nodeName == 'if')
    if options.flag_verbose_conditionals:
      print "Found", str(self)
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

    def sort_key(elt):
      (n,c) = elt
      # TODO: add linenum to function nodes in XML
      if n.nodeName == 'function':
        return 0
      # TODO: add linenum to else nodes in XML
      if n.getAttribute('nodeName') == 'else':
        attr_key = 'iflinenum'
      else:
        attr_key = 'linenum'
      return int(n.getAttribute(attr_key))

    if node:
      self.codeblocks = map(lambda x: CodeBlock(*x), sorted(traverse(node, None, conds), key=sort_key))
      self.loops = map(lambda x: Loop(*x), sorted(traverse(node, 'loop', conds), key=sort_key))

  def collect(self, f):
    return Collection(colls = map(lambda x: x.collect(f), self.codeblocks) + \
                              map(lambda x: x.collect(f), self.loops))
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
  def collect(self, f):
    return self.body.collect(f).loop(self)
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
  def collect(self, f):
    return self.body.collect(f)


class XMLParser(object):
  """Contains information on enclosed modules."""
  slots = ['doc', 'functions']
  def __init__(self, filename):
    assert type(filename) == type('')
    self.doc = xml.dom.minidom.parse(filename)
    program = getChildren(self.doc, 'program')[0]
    self.functions = map(Function, getChildren(program, 'function'))

class KeyValXMLParser(object):
  __slots__  = ['doc', 'items']
  def __init__(self, filename, val_type=str):
    assert type(filename) == type('')
    self.doc = xml.dom.minidom.parse(filename)
    self.parse_items(val_type)
  def parse_items(self, val_type):
    self.items = []
    machine = self.doc.firstChild
    for prop in getChildren(machine, 'prop'):
      key = str(prop.getAttribute('key'))
      value = val_type(prop.getAttribute('val'))
      self.items.append((key, value))

# parser for XML generated by Polly-ExaSAT Scop Analysis

class PollyArray(object):
  """Contains information about an array."""
  __slots__ = ['name', 'elt_type', 'elt_size', 'dim_n', 'dims']
  def __init__(self, node):
    self.name = str(node.getAttribute('name'))
    assert self.name[0:7] == 'MemRef_'
    self.name = self.name[7:]
    self.elt_type = str(node.getAttribute('elt_type'))
    self.elt_size = int(node.getAttribute('elt_size'))
    self.dim_n = int(node.getAttribute('dim_n'))
    self.dims = []
    for (i, dnode) in enumerate(getChildren(node, 'dim')):
      assert dnode.getAttribute('index') == str(i)
      dim = str(dnode.getAttribute('size')).lstrip('%')
      self.dims.append(dim)
    if self.dim_n == 0:
      dprint("  Parsing scalar (%s, %s, %s) ..." % \
        (self.name, self.elt_type, self.elt_size))
    else:
      dprint("  Parsing array (%s, %s, %s, %s, %s) ..." % \
        (self.name, self.elt_type, self.elt_size, self.dim_n, self.dims))

class PollyAccess(object):
  """Contains information about a statement."""
  __slots__ = ['type', 'reduction_type', 'is_scalar', 'relation']
  def __init__(self, node):
    self.type = str(node.getAttribute('type'))
    self.reduction_type = str(node.getAttribute('reduction_type'))
    self.is_scalar = str(node.getAttribute('is_scalar'))
    self.relation = str(node.getAttribute('relation'))
    dprint("    Parsing access (%s, %s, %s, %s) ..." % \
      (self.type, self.reduction_type, self.is_scalar, self.relation))

class PollyStatement(object):
  """Contains information about a statement."""
  __slots__ = ['name', 'domain', 'schedule', 'accesses']
  def __init__(self, node):
    self.name = str(node.getAttribute('name'))
    self.domain = str(node.getAttribute('domain'))
    self.schedule = str(node.getAttribute('schedule'))
    dprint("  Parsing statement (%s, %s, %s) ..." % \
      (self.name, self.domain, self.schedule))
    self.accesses = map(PollyAccess, getChildren(node, 'access'))

class PollyScop(object):
  """Contains information about a scop region."""
  __slots__ = ['filename', 'linenum', 'function', 'region', 'depth', 'arrays', 'statements']
  def __init__(self, node):
    self.filename = str(node.getAttribute('filename'))
    self.linenum = str(node.getAttribute('linenum'))
    self.function = str(node.getAttribute('function'))
    self.region = str(node.getAttribute('region'))
    self.depth = str(node.getAttribute('loop_depth'))
    dprint("Parsing scop (%s, %s, %s, %s, %s) ..." % \
      (self.filename, self.linenum, self.function, self.region, self.depth))
    self.arrays = map(PollyArray, getChildren(node, 'array'))
    self.statements = map(PollyStatement, getChildren(node, 'statement'))

class PollyXMLParser(object):
  """Contains information produced by Polly-ExaSAT Scop Analysis."""
  __slots__ = ['doc', 'scops']
  def __init__(self, filename):
    assert type(filename) == type('')
    self.doc = xml.dom.minidom.parse(filename)
    program = getChildren(self.doc, 'program')[0]
    self.scops = map(PollyScop, getChildren(program, 'scop'))
