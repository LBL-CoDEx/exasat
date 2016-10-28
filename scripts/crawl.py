#!/usr/bin/python

import sys
import os
import re
import pprint
import shlex
import subprocess

home_dir = os.environ.get("HOME")
gcd_binary = os.getenv("GCD_BINARY",
                       os.path.join(home_dir, "projects/exasat/src/genCodeDescript"))

class MapEntry:
  """ Struct containing a filename and a list of names of modules used """
  __slots__ = ['filename', 'uses']
  def __init__(self, filename, uses):
    self.filename = filename
    self.uses = uses
  def __repr__(self):
    return pprint.pformat({'filename': self.filename, 'uses': self.uses})

def read_modules(filename, mapping):
  """ Reads the first (and only the first) module in a file, and adds
  (k, v) to the mapping dict, where k is the name of the module, and v
  is a list of modules used within. """
  f = open(filename, 'r')
  lines = f.readlines()
  f.close()
  module_name = None
  for line in lines:
    line = line.rstrip('\r\n')
    if not module_name:
      m = re.search("^\s*module ([^ ]+)", line)
      if m:
        module_name = m.group(1)
        if module_name == "procedure":
          continue
        if module_name in mapping:
          print "WARNING: duplicate module %s detected!" % m.group(1)
          return
        print "  " + module_name 
        mapping[module_name] = MapEntry(filename, [])
    else:
      m = re.search("^\s*use ([^ ,]+)", line)
      if m:
        use_name = m.group(1)
        if use_name not in mapping[module_name].uses:
          mapping[module_name].uses.append(use_name)
          print "    uses " + use_name
      if re.search("^\s*end module", line):
        return

def crawl_dir(dirname, mapping):
  """ Run read_modules over every .f90 file in dirname """
  for filename in os.listdir(dirname):
    if re.match(".*\.f90$", filename):
      filepath = os.path.join(dirname, filename)
      print filename
      read_modules(filepath, mapping)
      print
  print

def generate_files(source_file):
  """ Run GCD on the source file. """
  cmd = "%s %s" % (gcd_binary, source_file)
  print cmd
  f = open(os.path.split(source_file)[1] + ".xml", "w")
  subprocess.Popen(shlex.split(cmd), stdout=f).communicate()[0]
  f.close()

def traverse_module(module_name, mapping, indent=1):
  """ Recursively run GCD on the file containing module_name.
  Module dependencies have GCD run on them first. """
  if module_name not in mapping:
    print "ERROR: module %s not found in search path!" % module_name
    sys.exit(1)
  filename = mapping[module_name].filename
  print "  " * indent + "uses %s defined in %s" % (module_name, filename)
  for use_name in mapping[module_name].uses:
    traverse_module(use_name, mapping, indent+1)

  if not os.path.isfile(module_name + ".rmod"):
    generate_files(filename)

def traverse_file(filename, mapping):
  """ Recursively run GCD on the file and its dependencies. """
  uses = set()
  print "Traversing %s" % filename
  f = open(filename, 'r')
  lines = f.readlines()
  f.close()
  for line in lines:
    line = line.rstrip('\r\n')
    m = re.search("^\s*use ([^ ,]+)", line)
    if m:
      use_name = m.group(1)
      if use_name not in uses:
        uses.add(use_name)
        traverse_module(use_name, mapping)
  generate_files(filename)

if __name__ == "__main__":

  if len(sys.argv) < 3:
    print "Usage: %s <filename> [<include-dirs> ...]" % sys.argv[0]
    sys.exit(1)

  # gather mapping of modules to file names
  mapping = {}
  (target, include_dirs) = (sys.argv[1], sys.argv[2:])
  for include_dir in include_dirs:
    print include_dir
    crawl_dir(include_dir, mapping)

  # generate the XML code description for the target file
  traverse_file(target, mapping)
