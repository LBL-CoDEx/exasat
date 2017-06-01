#!/usr/bin/python

import sys
import os
import subprocess
import shlex
import re
import string

def call(cmd, outfile="-"):
  cmd_list = shlex.split(cmd)
  if outfile != "-":
    f = open(outfile, 'w')
    subprocess.call(cmd_list, stdout=f)
    f.close()
  else:
    subprocess.call(cmd_list)

# tried using this, but output doesn't get redirected, so can't capture
#   cmd_tem = "clang %s %s -O3 -mllvm -polly -c -mllvm -debug-only=polly-scops -mllvm -polly-process-unprofitable"

def gen_scops_file(filename, defines):
  src_file = filename + ".c"
  asm_file = filename + ".s"
  asm_opt_file = filename + ".s.opt.ll"
  scops_file = filename + ".scops"

  cmd_tem = "clang -S -emit-llvm %s %s -o %s"
  call(cmd_tem % (defines, src_file, asm_file))

  cmd_tem = "opt -S -polly-canonicalize %s"
  call(cmd_tem % (asm_file), asm_opt_file)

  cmd_tem = "opt -q -basicaa -polly-process-unprofitable -polly-scops -analyze %s"
  call(cmd_tem % (asm_opt_file), scops_file)

  return scops_file

def get_sections(pattern, in_buf):
  rx = re.compile(pattern, re.M)

  sec_start = None
  search_start = 0

  m = rx.search(in_buf, search_start)
  while m:
    if sec_start != None:
      sec_end = m.start()
      yield in_buf[sec_start:sec_end]
    sec_start = m.start()
    search_start = m.end()
    m = rx.search(in_buf, search_start)

  if sec_start != None:
    sec_end = len(in_buf)
    yield in_buf[sec_start:sec_end]

def get_scops(scops_file):
  f = open(scops_file)
  file_buffer = f.read()
  f.close()
  return get_sections("^\s*Function:\s*(\w+)\s*$", file_buffer)

class Variable(object):
  __slots__ = ["name", "type", "sizes"]
  def __init__(self, n, t, s):
    self.name = n
    self.type = t
    self.sizes = s
  def __str__(self):
    result = ""
    result += "Name: " + self.name + '\n'
    result += "Type: " + self.type + '\n'
    result += "Sizes: (" + string.join(self.sizes, ", ") + ")" + '\n'
    return result

def get_memrefs(scop_output):
  rx_str = "\s*Arrays\s*{\s*\n([^}]*)\s*}\s*\n"
  m = re.search(rx_str, scop_output)
  if not m:
    raise Exception("missing memrefs")
  memrefs_str = m.group(1).splitlines()
  memref_rx_str = "(?P<mtype>\w+) MemRef_(?P<name>\w+)(?P<sizes>.*);"
  result = []
  for memref_str in memrefs_str:
    m = re.search(memref_rx_str, memref_str)
    if m:
      mtype = m.group('mtype')
      name = m.group('name')
      if m.group('sizes') == '':
        sizes = []
      else:
        sizes = m.group('sizes').strip("[]").split("][")
        sizes = map(lambda x: x.lstrip('%'), sizes)
      result.append(Variable(name, mtype, sizes))
  return result
        

class Statement(object):
  __slots__ = ['name', 'parameters', 'ind_vars', 'bounds', 'schedule', 'accesses']
  def __init__(self, n, p, i, b, s, a):
    self.accesses = []
    self.name = n
    self.parameters = p
    self.ind_vars = i
    self.bounds = b
    self.schedule = s
    self.accesses = a
  def __str__(self):
    result = ''
    result += "Name: " + self.name + '\n'
    result += "Parameters: " + str(self.parameters) + '\n'
    result += "Ind Vars: " + str(self.ind_vars) + '\n'
    result += "Bounds: " + str(self.bounds) + '\n'
    result += "Schedule: " + self.schedule + '\n'
    result += "Accesses:\n"
    for access in self.accesses:
      result += "  " + str(access) + '\n'
    return result

class Access(object):
  __slots__ = ['name', 'scalar', 'idxs', 'conds', 'reads', 'writes']
  def __init__(self, n, s, i, c, r, w):
    self.name = n
    self.scalar = s
    self.idxs = i
    self.conds = c
    self.reads = r
    self.writes = w
  def __str__(self):
    result = string.join( \
               map(str, [self.name, \
                         self.scalar, \
                         self.idxs, \
                         self.conds, \
                         self.reads, \
                         self.writes]), \
               " : ")
    return result

def parse_domain(stmt_name, domain_str):
  domain_rx_str = "\[(?P<params>[^\]]*)\] -> { %s\[(?P<ind_vars>[^\]]*)\] : (?P<bounds>.*) }" % stmt_name
  m = re.match(domain_rx_str, domain_str)
  if not m:
    raise Exception('domain string mismatch')
  params = m.group('params').split(", ")
  ind_vars = m.group('ind_vars').split(", ")
  bounds = m.group('bounds').split(" and ")
  return (params, ind_vars, bounds)

def parse_access(label, stmt_name, ind_vars, access_str):
  label_to_kind = {"ReadAccess":      (1, 0),
                   "MayWriteAccess":  (0, 1),
                   "MustWriteAccess": (0, 1)}
  (is_read, is_write) = label_to_kind[label]
  access_rx_str = "\[Reduction Type: (?P<rdx>\w+)\] \[Scalar: (?P<sca>\d)\]\s+(?P<params>.*) -> { (?P<refs>.*) }"
  m = re.match(access_rx_str, access_str)
  if not m:
    raise Exception("missing accesses")
  (rdx, scalar, refs_str) = m.group('rdx', 'sca', 'refs')
  ref_rx_str = "\s*%s\[%s\] -> MemRef_(?P<aname>\w+)\[(?P<idxs>[^;]*)\](?: : (?P<conds>[^;]*))?" % (stmt_name, string.join(ind_vars, ", "))
  refs = re.findall(ref_rx_str, refs_str)
  (aname, idxs, conds) = (None, "", "")
  for ref in refs:
    (name, idxs, conds) = ref
    idxs = tuple(idxs.split(", "))
    conds = conds.split(" and ")
    if aname == None:
      aname = name
    assert name == aname # name should match in all refs
    # TODO: choose the correct idx/cond pair
    #       for now, last one seems correct
  return Access(aname, scalar, idxs, conds, is_read, is_write)

def group_accesses(accesses):
  seen = {}
  for a in accesses:
    key = (a.name, a.idxs)
    if key not in seen:
      seen[key] = a
    else:
      seen[key].reads += a.reads
      seen[key].writes += a.writes
  return sorted(seen.itervalues(), key=lambda x: x.name)

def parse_statement(stmt_str):
  m = re.match("\s*(?P<stmt_name>\w+)\s*\n(?P<attrs_str>.*)", stmt_str, re.S)
  stmt_name = m.group('stmt_name')
  attrs_str = m.group('attrs_str')
  attr_rx_str = "\s*(?P<label>\w+)\s*:=\s*(?P<val>.*?);\s*\n"
  attrs = re.findall(attr_rx_str, attrs_str, re.S)
  access_strs = ["ReadAccess", "MayWriteAccess", "MustWriteAccess"]
  accesses = []
  for attr in attrs:
    (label, val) = attr
    if label == "Domain":
      (params, ind_vars, bounds) = parse_domain(stmt_name, val)
    elif label == "Schedule":
      schedule = val
    elif label in access_strs:
      accesses.append(parse_access(label, stmt_name, ind_vars, val))
    else:
      raise Exception("unknown attribute type: " + label)
  accesses = group_accesses(accesses)
  return Statement(stmt_name, params, ind_vars, bounds, schedule, accesses)

def get_statements(scop_output):
  result = []
  stmt_rx_str = "\s*Statements\s*{\s*\n(.*)\s*}\s*\n"
  m = re.search(stmt_rx_str, scop_output, re.S)
  if not m:
    raise Exception("missing statements")
  # use indentation (ick) to split statement blocks
  stmts_str = get_sections("^    \t\w+\s*$", m.group(1))
  for stmt_str in stmts_str:
    result.append(parse_statement(stmt_str))
  return result

def analyze_scop(scop_output):
  print "Variables:\n"
  for memref in get_memrefs(scop_output):
    print memref
  print "Statements:\n"
  for statement in get_statements(scop_output):
    print statement

def main(args):
  filename = os.getenv("filename", "operators.fv2")
  scops_file = filename + ".scops"
  if (not os.path.exists(scops_file)):
    defines = os.getenv("defines", "-DUSE_GSRB -DGSRB_BRANCH -DUSE_HELMHOLTZ")
    scops_file = gen_scops_file(filename, defines)
  for scop_output in get_scops(scops_file):
    analyze_scop(scop_output)

if __name__ == "__main__":
  main(sys.argv)
