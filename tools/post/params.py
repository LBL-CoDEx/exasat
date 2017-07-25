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

def to_sym_dict(list_of_pairs):
  return dict(map(lambda (x,y): (parse_expr(x), parse_expr(y)), list_of_pairs))

def arrayName(name, component, namesubs):
  if component != '':
    component = parse_expr(component).subs(namesubs)
    return '%s.%s' % (name, component)
  else:
    return name

def doSymRepl(expr, repl):
  if type(expr) == type('') or type(expr) == unicode:
    expr = parse_expr(expr)
  return expr.xreplace(repl)

def parseExpr(s, symsubs):
  try:
    return int(s) # try to do the fast thing first
  except ValueError as e:
    return doSymRepl(s, symsubs)

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
      result = int(doSymRepl(x, params))
    except Exception as e:
      print "Could not do integer parameter substitution for expression: ", x
      print "Parameters available: ", params
      raise e
  return result
