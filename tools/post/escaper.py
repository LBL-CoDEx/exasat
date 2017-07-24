#!/usr/bin/env python

""" Helper program for escaping the special characters in an XML file. """
__author__ = "Cy Chan"
__copyright__ = "Copyright 2016, The Regents of the University of California, through Lawrence Berkeley National Laboratory"
__credits__ = ["Cy Chan"]
__license__ = "Modified BSD License (see LICENSE file)"
__version__ = "2.0"
__maintainer__ = "Cy Chan"
__email__ = "cychan@lbl.gov"
__status__ = "Production"

import sys
import re
from xml.sax.saxutils import escape

def main(args):
  infile = open(args[1], 'r')
  outfile = open(args[2], 'w')

  line = infile.readline()
  while line:

    substrings = line.split('"')
    values = substrings[1::2]
    substrings[1::2] = map(escape, values)
    line = '"'.join(substrings)
    outfile.write(line)

    line = infile.readline()

  infile.close()
  outfile.close()

if __name__ == "__main__":
  if len(sys.argv) < 3:
    print "Usage: <me> <input-xml> <output-xml>"
    sys.exit(1)
  main(sys.argv)
