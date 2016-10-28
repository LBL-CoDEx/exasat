#!/bin/env python

import sys
import os
import subprocess
import shlex

exe = '../../src/genCodeDescript'

def path_parts(path):
  (a, b) = os.path.split(path)
  if a == '':
    return [b]
  elif a == '/':
    return [a, b]
  else:
    return path_parts(a) + [b]

def rebase(src, src_base_dir, dst_base_dir):
  src_path = path_parts(src)
  src_base_path = path_parts(src_base_dir)
  assert(src_path[0:len(src_base_path)] == src_base_path)
  src_sub_path = src_path[len(src_base_path):]
  return os.path.join(dst_base_dir, *src_sub_path)

def main(src_base_dir, dst_base_dir):
  for (src_dir, sub_dirs, files) in os.walk(src_base_dir):
    dst_dir = rebase(src_dir, src_base_dir, dst_base_dir)
    if not os.path.exists(dst_dir):
      os.makedirs(dst_dir)
    for src_file in files:
      (src_root, src_ext) = os.path.splitext(src_file)
      if src_ext == '.f90':
        dst_file = src_root + '.xml'
        src_path = os.path.join(src_dir, src_file)
        dst_path = os.path.join(dst_dir, dst_file)
        if not os.path.exists(dst_path):
          cmd = '%s %s' % (exe, src_path)
          print cmd + ' > ' + dst_path
          dst_path_f = open(dst_path, 'w')
          env = dict(os.environ)
          env["NOFILTER"]="1"
          subprocess.call(shlex.split(cmd), stdout=dst_path_f, env=env)
          dst_path_f.close()

if __name__ == "__main__":

  src_base_dir = 'inputs'
  if len(sys.argv) > 1:
    src_base_dir = sys.argv[1]

  dst_base_dir = 'xml_new'
  if len(sys.argv) > 2:
    dst_base_dir = sys.argv[2]

  main(src_base_dir, dst_base_dir)
