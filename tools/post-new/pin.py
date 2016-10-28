DIM = 3
TEMPFILE = 'output.temp.csv'

class options:
  flag_debug = False

class PinAnalysis:

  def __init__(self, doc = None):
    self.dim = DIM
    if type(doc) == type(''):
      doc = xml.dom.minidom.parse(doc)
    if doc:
      self.info = self.docToInfo(doc)

  def docToInfo(self, doc):

    info = {}
    for fnode in getChildren(doc, 'function'):

      finfo = {}
      fname = str(fnode.getAttribute('name'))
      if not re.match('^advance', fname):
        continue
      filename = str(fnode.getAttribute('file'))

      loops = []
      for lnode in getChildren(fnode, 'loop'):
        loop = {}
        loop['fname'] = fname
        loop['file'] = filename
        loop['num'] = str(lnode.getAttribute('num'))
        if not re.match('^\d+\.\d+\.\d+$', loop['num']):
          continue
        loop['num'] = "'" + loop['num']
        loop['line'] = int(lnode.getAttribute('line'))
        loop['entries'] = int(lnode.getAttribute('entries'))
        loop['iters'] = int(lnode.getAttribute('iters'))
        loop['ops'] = int(lnode.getAttribute('ops'))
        loop['reads'] = int(lnode.getAttribute('reads'))
        loop['writes'] = int(lnode.getAttribute('writes'))
        loop['bytereads'] = int(lnode.getAttribute('bytereads'))
        loop['bytewrites'] = int(lnode.getAttribute('bytewrites'))
        loop['uniquereads'] = int(lnode.getAttribute('uniquereads'))
        loop['uniquewrites'] = int(lnode.getAttribute('uniquewrites'))
        loop['adds'] = int(lnode.getAttribute('adds'))
        loop['subs'] = int(lnode.getAttribute('subs'))
        loop['muls'] = int(lnode.getAttribute('muls'))
        loop['divs'] = int(lnode.getAttribute('divs'))
        loop['cachehits'] = int(lnode.getAttribute('cachehits'))
        loop['cachemisses'] = int(lnode.getAttribute('cachemisses'))
        loops.append(loop)

      finfo['loops'] = loops
      info[fname] = finfo

    return info

  def printInfo(self):
    pp.pprint(self.info)

  def exportSpreadsheet(self, outfile):

    def writeTable(f):
      f.write('PIN ANALYSIS\n')
      f.write('||||||Actual||Byte||Unique||Flops||||Cache||\n')
      f.write('Function|Num|Line|Entries|Iters|Ops|R|W|R|W|R|W|Adds|Subs|Muls|Divs|Hits|Misses|\n')

      # write summary info
      for fname in sorted(self.info.keys()):
        function = self.info[fname]
        for loop in function['loops']:

          def lookup(loop, *args):
            return tuple(map(lambda x:loop[x], args))

          f.write(('%s|'*2+'%d|'*16) % lookup(loop, 'fname', 'num', 'line', 'entries', \
                  'iters', 'ops', 'reads', 'writes', 'bytereads', 'bytewrites', 'uniquereads', \
                  'uniquewrites', 'adds', 'subs', 'muls', 'divs', 'cachehits', 'cachemisses'))
          f.write('\n')
        f.write('\n')

    def writePerIterTable(f):
      f.write('Per Iteration:\n')
      f.write('||||||Actual||Byte||Bytes/Actual||Flops||||\n')
      f.write('Function|Num|Line|Entries|Iters|Ops|R|W|R|W|R|W|Adds|Subs|Muls|Divs|\n')

      # write summary info
      for fname in sorted(self.info.keys()):
        function = self.info[fname]
        for loop in function['loops']:

          def lookup1(loop, *args):
            return tuple(map(lambda x:loop[x], args))
          def lookup2(loop, *args):
            return tuple(map(lambda x:loop[x]/loop['iters'], args))

          f.write(('%s|'*2+'%d|'*3) % lookup1(loop, 'fname', 'num', 'line', 'entries', 'iters'))
          f.write(('%d|'*5) % lookup2(loop, 'ops', 'reads', 'writes', 'bytereads', 'bytewrites'))
          f.write(('=%s/%s|' % getRefs(-2, -4)) * 2)
          f.write(('%d|'*4) % lookup2(loop, 'adds', 'subs', 'muls', 'divs'))
          f.write('\n')
        f.write('\n')

    f = open(TEMPFILE, 'wt')
    writeTable(f)
    writePerIterTable(f)
    f.close()

    spreadsheet.processRefs(outfile, TEMPFILE)
    os.remove(TEMPFILE)

def main(args):
  if len(args) < 3:
    print "Usage: <me> <pinxml-input> <pintsv-output>"
    sys.exit(1)

  (pinxml, pintsv) = args[1:3]

  pa = PinAnalysis(pinxml)
  if options.flag_debug:
    pa.printInfo()
  pa.exportSpreadsheet(pintsv)
 
if __name__ == '__main__':
  main(sys.argv)
