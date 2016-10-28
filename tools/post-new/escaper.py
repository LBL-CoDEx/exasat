#!/usr/bin/env python

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
