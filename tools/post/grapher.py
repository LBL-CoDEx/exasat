import re

from common import options
import analyze

noGraphs = True
noSlideShow = True
hasGV = False
if not noGraphs:
  try:
    import pygraphviz as pgv
    hasGV = True
  except Exception as e:
    print e

class Grapher(object):

  def __init__(self, sa):
    self.info = sa.info

  # draw graph showing dependencies between loops only
  def generateLoopDepGraphs(self):
    if not hasGV:
      print "Skipping graph generation ..."
      return
    for fname in self.info.keys():
      loops = self.info[fname]['loops']
      G = pgv.AGraph(name = fname, directed = True)
      for (i, loop1) in enumerate(loops):
        l1Arrays = loop1['arrays']
        l1 = 'L%d' % (i+1)
        G.add_node(l1)
        for (j, loop2) in enumerate(loops[i+1:]):
          l2Arrays = loop2['arrays']
          j += i+1
          l2 = 'L%d' % (j+1)
          if set.intersection(analyze.getRNames(l1Arrays), analyze.getRWNames(l2Arrays)) or \
             set.intersection(analyze.getRNames(l1Arrays), analyze.getWNames(l2Arrays)) or \
             set.intersection(analyze.getRWNames(l1Arrays), analyze.getRNames(l2Arrays)) or \
             set.intersection(analyze.getRWNames(l1Arrays), analyze.getRWNames(l2Arrays)) or \
             set.intersection(analyze.getRWNames(l1Arrays), analyze.getWNames(l2Arrays)) or \
             set.intersection(analyze.getWNames(l1Arrays), analyze.getRNames(l2Arrays)) or \
             set.intersection(analyze.getWNames(l1Arrays), analyze.getRWNames(l2Arrays)) or \
             set.intersection(analyze.getWNames(l1Arrays), analyze.getWNames(l2Arrays)):
            G.add_edge(l1, l2)
      G.tred() # compute transitive reduction
      G.layout(prog = 'dot')
      G.draw(fname + '.png')


  # draw dependency graphs showing both loops and data
  def generateDepGraphs(self):

    loopNodeProps = {'shape': 'oval', 'style': 'filled', 'fillcolor': 'cornsilk'}
    dataNodeProps = {'shape': 'box', 'style': 'filled', 'fillcolor': 'cadetblue1'}
    edgeProps = {}
    ghostEdgeProps = {'style': 'bold', 'color': 'red'}
    subGraphProps = {'shape': 'box', 'name': 'sub', 'style': 'filled', 'color': 'lightgrey'}
    graphIteration = 0

    def nodeType(n):
      if re.match('L\d+', n):
        return 'loop'
      else:
        return 'data'

    # TODO: this is really inefficient
    def reachable(a, b):
      if a == b:
        return True
      for (t, data) in G.out_edges(a):
        for (t, c) in G.out_edges(data):
          if reachable(c, b):
            return True
      return False

    def iterpairs(x):
      for i in xrange(0, len(x)):
        for j in xrange(i+1, len(x)):
          yield (x[i], x[j])

    def countArrays(l):
      return sum(map(lambda x: len(x.split(',')), l))

    def findMerge(G):
      loopNodes = filter(lambda x: nodeType(x) == 'loop', G.nodes())
      m = 0
      toMerge = None
      for (n1, n2) in iterpairs(sorted(loopNodes)):
        if not (reachable(n1, n2) or reachable(n2, n1)):
          s1 = set(map(lambda x: x[0], G.in_edges(n1)))
          s2 = set(map(lambda x: x[0], G.in_edges(n2)))
          common = s1.intersection(s2)
          if countArrays(common) > m:
            m = countArrays(common)
            toMerge = (n1, n2)
      return toMerge

    def mergeNodes(G, n1, n2, nodeProps = {}, edgeProps = {}):
        if options.flag_debug:
          print "Merging: %s, %s" % (n1, n2)
        newNode = "%s,%s" % (n1, n2)
        G.add_node(newNode, **nodeProps)
        map(lambda x:G.add_edge(x[0], newNode, **edgeProps), G.in_edges(n1))
        map(lambda x:G.add_edge(x[0], newNode, **edgeProps), G.in_edges(n2))
        map(lambda x:G.add_edge(newNode, x[1], **edgeProps), G.out_edges(n1))
        map(lambda x:G.add_edge(newNode, x[1], **edgeProps), G.out_edges(n2))
        G.remove_node(n1)
        G.remove_node(n2)

    def printGraph(G):
      if noSlideShow:
        return
      if not os.path.exists('png.out'):
        os.mkdir('png.out')
      if not hasattr(printGraph, 'counter'):
        printGraph.counter = 0
      printGraph.counter += 1
      G.layout(prog = 'dot')
      G.draw('png.out/printGraph%02d.png' % printGraph.counter)

    def horizontalMergeLoops(G):
      printGraph(G)
      while True:
        temp = findMerge(G)
        if not temp:
          break
        (n1, n2) = temp
        mergeNodes(G, n1, n2, nodeProps = loopNodeProps)
        printGraph(G)

    def mergeDataNodes(G):
      printGraph(G)
      success = True
      while success:
        success = False
        dataNodes = filter(lambda x: nodeType(x) == 'data', G.nodes())
        for (n1, n2) in iterpairs(sorted(dataNodes)):
          s1 = set(map(lambda x: x[0], G.in_edges(n1)))
          s2 = set(map(lambda x: x[0], G.in_edges(n2)))
          s3 = set(map(lambda x: x[1], G.out_edges(n1)))
          s4 = set(map(lambda x: x[1], G.out_edges(n2)))
          if s1 == s2 and s3 == s4:
            mergeNodes(G, n1, n2, nodeProps = dataNodeProps)
            printGraph(G)
            success = True
            break

    # construct raw graph
    def makeGraph(finfo):
      G = pgv.AGraph(name = fname, directed = True)
      for (i, loop) in enumerate(finfo['loops']):
        lname = 'L%d' % (i+1)
        G.add_node(lname, **loopNodeProps)
        for aname in analyze.getRNames(loop['arrays']):
          G.add_node(aname, **dataNodeProps)
          G.add_edge(aname, lname, **edgeProps)
        for aname in analyze.getRWNames(loop['arrays']):
          G.add_node(aname, **dataNodeProps)
          G.add_edge(lname, aname, dir='both', **edgeProps)
        for aname in analyze.getWNames(loop['arrays']):
          G.add_node(aname, **dataNodeProps)
          G.add_edge(lname, aname, **edgeProps)
      return G

    # clear labels
    def clearLabels(G):
      for n in G.iternodes():
        n.attr['label'] = ''

    if not hasGV:
      print "Skipping graph generation ..."
      return

    for (fname, finfo) in self.info.iteritems():
      # make and draw raw graph
      G = makeGraph(finfo)
      G.layout(prog = 'dot')
      G.draw(fname + '_0.png')

      # optimize and draw graph
      mergeDataNodes(G)
      G.layout(prog = 'dot')
      G.draw(fname + '_1.png')
      horizontalMergeLoops(G)
      mergeDataNodes(G)
      G.layout(prog = 'dot')
      G.draw(fname + '_2.png')
