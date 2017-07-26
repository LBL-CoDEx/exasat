import os
import string
import xml.dom.minidom
import ast

from distutils.util import strtobool
from common import options, int_types, float_types
from params import doSymSubs

# can override this if needed
flag_old = strtobool(os.getenv('old', 'True'))

def getChildren(node, tag):
  return filter(lambda x: x.nodeName == tag, node.childNodes)

def arrayName(name, component):
  if component != '':
    component = doSymSubs(component)
    return '%s.%s' % (name, component)
  else:
    return name

def processBlock(node):
  result = {}
  # process flops
  result['flops'] = {'adds'       : int(node.getAttribute('adds')), \
                     'multiplies' : int(node.getAttribute('multiplies')), \
                     'divides'    : int(node.getAttribute('divides')), \
                     'specials'   : int(node.getAttribute('specials'))}

  # process scalars
  result['scalars'] = []
  for scalarNode in getChildren(node, 'scalar'):
    scalar = {'name' : str(scalarNode.getAttribute('name')), \
              'datatype' : str(scalarNode.getAttribute('datatype')), \
              'const': str(scalarNode.getAttribute('isConstant')), \
              'reads': int(scalarNode.getAttribute('reads')), \
              'writes': int(scalarNode.getAttribute('writes')), \
             }
    accesstype = str(scalarNode.getAttribute('accesstype'))
    if accesstype == 'readonly' and scalar['writes'] > 0 or \
       accesstype == 'writeonly' and scalar['reads'] > 0 or \
       (accesstype == 'readwrite' or accesstype == 'writeread') and \
       (scalar['reads'] == 0 or scalar['writes'] == 0):
      raise Exception('scalar %s mismatched accesstype: %s, reads: %d, writes: %d' % \
                      (scalar['name'], accesstype, scalar['reads'], scalar['writes']))
    result['scalars'].append(scalar)

  # process arrays
  result['stateArrays'] = []
  result['arrays'] = []
  for arrayNode in getChildren(node, 'array'):
    name_raw = str(arrayNode.getAttribute('name'))
    type_raw = str(arrayNode.getAttribute('datatype'))
    # only analyze numerical data arrays
    if not (type_raw in float_types or
            type_raw in int_types):
      print 'Warning: ignoring array "%s" of type %s' % (name_raw, type_raw)
      continue
    access = {}
    for accessNode in getChildren(arrayNode, 'access'):
      index = parseStringTuple(accessNode.getAttribute('index' if flag_old else 'offset'))
      if not flag_old and len(index) == options.dim+1:
        index = index[:-1] # new format keeps component dimension, so drop it
      index = tuple(map(int, index))
      access[index] = {'reads' : int(accessNode.getAttribute('reads' )), \
                       'writes': int(accessNode.getAttribute('writes'))}
      # spatial (non-state) arrays have 3D relative indices
      if len(index) == options.dim and isRelative(accessNode):
        arraytype = 'array'
      # classify state arrays as relative or absolute
      elif isRelative(accessNode):
        arraytype = 'relStateArray'
      else:
        arraytype = 'absStateArray'
    entry = {'name' : arrayName(str(arrayNode.getAttribute('name')), \
                                str(arrayNode.getAttribute('component'))), \
             'copies' : 1, \
             'arraytype' : arraytype, \
             'datatype' : str(arrayNode.getAttribute('datatype')), \
             'accesstype' : str(arrayNode.getAttribute('accesstype')), \
             'access' : access}

    if arraytype == 'array':
      result['arrays'].append(entry)
    else:
      result['stateArrays'].append(entry)

  return result

def processLoop(lnode):
  loop = processBlock(lnode)
  loop['loopvar'] = str(lnode.getAttribute('loopvar'))
  loop['linenum'] = int(lnode.getAttribute('linenum'))
  loop['range'] = (str(lnode.getAttribute('lowerbound')), str(lnode.getAttribute('upperbound')))
  loop['stride'] = int(lnode.getAttribute('stride'))
  return loop

def processIf(inode):
  result = processBlock(inode)
  result['linenum'] = int(inode.getAttribute('linenum'))
  result['conditional'] = str(inode.getAttribute('conditional'))
  return result

def isRelative(anode):
  if flag_old:
    return anode.getAttribute('isRelative') == "true"
  else:
    tup = parseStringTuple(anode.getAttribute('dependentloopvar'))
    return any(map(lambda x: x != '', tup))

def traverseLoops(node):
  loops = []
  for lnode in getChildren(node, 'loop'):
    loop = processLoop(lnode)
    loop['loops'] = traverseLoops(lnode)
    loop['ifs'] = traverseIfs(lnode)
    loops.append(loop)
  return loops

def traverseIfs(node):
  ifs = []
  for inode in getChildren(node, 'if'):
    result = processIf(inode)
    result['loops'] = traverseLoops(inode)
    result['ifs'] = traverseIfs(inode)
    ifs.append(result)
  return ifs

def parseStringTuple(s):
  s = s.strip()
  if s[0] != '(' or s[-1] != ')':
    raise IOError('input string does not represent a tuple: %s' % s)
  return tuple(map(lambda x: string.strip(str(x)),
                   string.split(s[1:-1], ',')))

def parseIntTuple(s):
  return tuple(map(int, parseStringTuple(s)))

def processVars(fnode, varType):
  result = set()
  for node in getChildren(fnode, varType):
    aname = arrayName(str(node.getAttribute('name')), \
                      str(node.getAttribute('component')))
    result.add(aname)
  return result

def processComms(fnode):
  comms = []
  for cnode in getChildren(fnode, 'comm'):
    comm = {}
    comm['linenum'] = int(cnode.getAttribute('linenum'))
    comm['interface'] = str(cnode.getAttribute('interface'))
    comm['arrays'] = []
    for anode in getChildren(cnode, 'array'):
      array = {}
      array['name'] = str(anode.getAttribute('name'))
      array['numComps'] = doSymSubs(str(anode.getAttribute('numofcomponents')))
      array['type'] = str(anode.getAttribute('commtype'))
      array['ghost'] = ast.literal_eval(anode.getAttribute('ghost')) # ick
      comm['arrays'].append(array)
    comms.append(comm)
  return comms

# merge ifs into this node
def flattenIfs(x):
  def mergeFlops(a, b):
    for k in a.iterkeys():
      a[k] += b[k]
  for y in x['loops'] + x['ifs']:
    flattenIfs(y)
  for y in x['ifs']:
    if 'flops' in x: 
      mergeFlops(x['flops'], y['flops'])
      x['scalars'].extend(y['scalars'])
      x['arrays'].extend(y['arrays'])
      x['stateArrays'].extend(y['stateArrays'])
    x['loops'].extend(y['loops'])
  del x['ifs']
  x['loops'].sort(key=lambda x: x['linenum'])

# process functions in node and append info onto finfos
def processFunctions(node, finfos):
  functions = getChildren(node, 'function')
  for fnode in functions:
    fname = str(fnode.getAttribute('name'))
    print 'Processing %s ...' % fname
    finfo = {}
    # process local/nonlocal, and read/write/read-write vars
    for varType in ['local', 'nonlocal']:
      finfo[varType] = processVars(fnode, varType)
    # process comms
    finfo['comms'] = processComms(fnode)
    # process loops
    finfo['loops'] = traverseLoops(fnode)
    # process ifs
    finfo['ifs'] = traverseIfs(fnode)
    # TODO: add better support for conditionals, for now just flatten
    flattenIfs(finfo)
    finfos.append((fname, finfo))

class XMLParser(object):

  slots = ['doc', 'info']

  def __init__(self, filename):
    assert type(filename) == type('')
    self.doc = xml.dom.minidom.parse(filename)

  def info(self):
    # for functions not in modules
    finfos = []
    processFunctions(self.doc.firstChild, finfos)

    # for functions in modules
    modules = getChildren(self.doc.firstChild, 'module')
    for mnode in modules:
      processFunctions(mnode, finfos)

    return dict(finfos)

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
