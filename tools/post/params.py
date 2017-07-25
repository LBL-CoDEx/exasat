#!/usr/bin/env python

""" Parameter handling functions and defaults for code analysis. """
__author__ = "Cy Chan"
__copyright__ = "Copyright 2016, The Regents of the University of California, through Lawrence Berkeley National Laboratory"
__credits__ = ["Cy Chan"]
__license__ = "Modified BSD License (see LICENSE file)"
__version__ = "2.0"
__maintainer__ = "Cy Chan"
__email__ = "cychan@lbl.gov"
__status__ = "Production"

from sympy.parsing.sympy_parser import parse_expr

model_params = [('X Problem Size', 128), \
                ('Y Problem Size', 128), \
                ('Z Problem Size', 128), \
                ('X Block Size', 16), \
                ('Y Block Size', 16), \
                ('Z Block Size', 128), \
                ('Num Ghosts (NG)', 4), \
                ('Num Species', 9), \
                ('GB/s/thread', 1), \
                ('Gflop/s/thread', 10), \
                ('Int Regs/thread', 16), \
                ('FP Regs/thread', 16), \
                ('Threads', 1024), \
                ('$/thread group (kB)', 1024), \
                ('Word Size', 8), \
                ('Cache Line Size', 64), \
                ( 'R cost', 8), \
                ('RW cost', 16), \
                ( 'W cost', 16), \
                ('Division Cost', 39), \
                ('Special Cost', 125), \
                ('Streaming Writes', False), \
                ('NTA Hints', False), \
                ('NIC BW (GB/s)', 100), \
                ('NIC Latency (us)', 0.4), \
                ('CPU/DRAM Overlap', 1), \
                ('PWFlag', False)]

def defaultParams():
  params = {}
  for (name, value) in model_params:
    params[name] = value
  return params

def customizeParams(params, overrides):
  for (key, value) in overrides.iteritems():
    if key not in params:
      raise Exception('parameter: %s not valid' % key)
    if (value != params[key]):
      print 'Overriding parameter %s = %s' % (key, value)
      params[key] = value

def to_sym_dict(list_of_pairs):
  return dict(map(lambda (x,y): (parse_expr(x), parse_expr(y)), list_of_pairs))

# name replacements for better naming
namesubs = [('qx1+n-1', 'qxn'), \
            ('qy1+n-1', 'qyn'), \
            ('qh1+n-1', 'qhn'), \
            ('iry1+n-1', 'iryn'), \
           ]

# symbolic replacements and constants
symsubs = [('dlo(1)', 'lo(1) - ng'), \
           ('dhi(1)', 'hi(1) + ng'), \
           ('dlo(2)', 'lo(2) - ng'), \
           ('dhi(2)', 'hi(2) + ng'), \
           ('dlo(3)', 'lo(3) - ng'), \
           ('dhi(3)', 'hi(3) + ng'), \
           ('ncons', 'nspecies + 5'), \
           ('nprim', '3 * nspecies + 7'), \
           ('iryn', 'iry1+n-1'), \
           ('qyn', 'qy1+n-1'), \
           ('qxn', 'qx1+n-1'), \
           ('qhn', 'qh1+n-1'), \
           ('qx1', 'qy1+nspecies'), \
           ('qh1', 'qy1+2*nspecies'), \
           ('irho', '1'), \
           ('imx', '2'), \
           ('imy', '3'), \
           ('imz', '4'), \
           ('iene', '5'), \
           ('iry1', '6'), \
           ('qrho', '1'), \
           ('qu', '2'), \
           ('qv', '3'), \
           ('qw', '4'), \
           ('qpres', '5'), \
           ('qtemp', '6'), \
           ('qe', '7'), \
           ('qy1', '8'), \
           ('k3d', '0'), \
           ('kc' , '0'), \
           ('k'  , '0'), \
           ('offset', '0'), \
          ]

symrepl = to_sym_dict(symsubs)

# NOTE: this is really slow, try to use doSymRepl instead
def doSymSubs(expr, repl = symsubs):
  if type(expr) == type('') or type(expr) == unicode:
    expr = parse_expr(expr)
  for r in repl:
    expr = expr.subs(r[0], r[1])
  return expr

def doNameSubs(expr):
  return doSymSubs(expr, namesubs)

def doSymRepl(expr, repl = symrepl):
  if type(expr) == type('') or type(expr) == unicode:
    expr = parse_expr(expr)
  return expr.xreplace(repl)
