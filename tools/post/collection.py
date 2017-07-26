#!/usr/bin/env python

""" Collection class for generic analysis traversal of code tree. """
__author__ = "Cy Chan"
__copyright__ = "Copyright 2016, The Regents of the University of California, through Lawrence Berkeley National Laboratory"
__credits__ = ["Cy Chan"]
__license__ = "Modified BSD License (see LICENSE file)"
__version__ = "2.0"
__maintainer__ = "Cy Chan"
__email__ = "cychan@lbl.gov"
__status__ = "Production"

class Collection(object):
  """Generic collection class for use with collect visitor.
    
     Items in collection are grouped by name, so must have a "name" field.
     When items with the same name are added, they are combined using the += operator.
     When the collection is looped over, it calls each item's loop() method.

     Items must implement __iadd__ (or __add__), and loop.
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
    # pass self into each item's loop() method in case each item needs to know
    # about other items in the collection (i.e. its siblings)
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
  def iterfirst(self):
    # return first item order of iteration
    for item in self.d.itervalues():
      return item
