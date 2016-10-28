#!/usr/bin/python

import pygraphviz as pgv

# note: subgraphs must have a name that starts with 'cluster' to be
# spatially clustered in the output 
graphs = [ \
  {'deps': {'ctoprim': {'in': ['*U(1:5)'], 'out': ['Q(2:6)'] }, \
            'diff_1' : {'in': [], 'out': ['D(1)']}, \
            'diff_2' : {'in': ['*Q(2:6)'], 'out': ['Q\'(1:9)']}, \
            'diff_3' : {'in': ['*Q(2:6)'], 'out': ['Q\'(1:9)']}, \
            'diff_4' : {'in': ['*Q(2:6)'], 'out': ['Q\'(1:9)']}, \
            'diff_5' : {'in': ['*Q(2:6)', '*Q\'(1:9)'], 'out': ['D(2:4)']}, \
            'diff_6' : {'in': ['*Q(2:6)', '*Q\'(1:9)'], 'out': ['D(2:4)']}, \
            'diff_7' : {'in': ['*Q(2:6)', '*Q\'(1:9)'], 'out': ['D(2:4)']}, \
            'diff_8' : {'in': ['*Q(2:6)', 'Q\'(1:9)', 'D(2:4)'], 'out': ['D(5)']}, \
            'hyp_1'  : {'in': ['*U(1:5)', '*Q(2:6)'], 'out': ['F(1:5)']}, \
            'hyp_2'  : {'in': ['*U(1:5)', '*Q(2:6)'], 'out': ['F(1:5)']}, \
            'hyp_3'  : {'in': ['*U(1:5)', '*Q(2:6)'], 'out': ['F(1:5)']}, \
            'advance': {'in': ['U(1:5)', 'F(1:5)', 'D(1)', 'D(2:4)', 'D(5)'], 'out': ['Unew(1:5)']}, \
#            'fill_boundary': {'in': ['Unew(1:5)'], 'out': ['Unew4(1:5)']}, \
           }, \
    'subgraphs': [{'name': 'cluster_diff_234', 'nbunch': ['diff_2', 'diff_3', 'diff_4']}, \
                  {'name': 'cluster_diff_567', 'nbunch': ['diff_5', 'diff_6', 'diff_7']}, \
                  {'name': 'cluster_hyp_123', 'nbunch': ['hyp_1', 'hyp_2', 'hyp_3']}, \
                  {'name': 'cluster_updateU_2', 'nbunch': ['diff_1', 'diff_8'], 'color': 'white'}], \
  }, \
  {'deps': {'ctoprim': {'in': ['*U(1:5)'], 'out': ['Q(2:6)'] }, \
            'diff_1' : {'in': [], 'out': ['D(1)']}, \
            'diff_2' : {'in': ['*Q(2:6)'], 'out': ['Q\'(1:9)']}, \
            'diff_3' : {'in': ['*Q(2:6)'], 'out': ['Q\'(1:9)']}, \
            'diff_4' : {'in': ['*Q(2:6)'], 'out': ['Q\'(1:9)']}, \
            'diff_5' : {'in': ['*Q(2:6)', '*Q\'(1:9)'], 'out': ['D(2:4)']}, \
            'diff_6' : {'in': ['*Q(2:6)', '*Q\'(1:9)'], 'out': ['D(2:4)']}, \
            'diff_7' : {'in': ['*Q(2:6)', '*Q\'(1:9)'], 'out': ['D(2:4)']}, \
            'diff_8' : {'in': ['*Q(2:6)', 'Q\'(1:9)', 'D(2:4)'], 'out': ['D(5)']}, \
            'hyp_1'  : {'in': ['*U(1:5)', '*Q(2:6)'], 'out': ['F(1:5)']}, \
            'hyp_2'  : {'in': ['*U(1:5)', '*Q(2:6)'], 'out': ['F(1:5)']}, \
            'hyp_3'  : {'in': ['*U(1:5)', '*Q(2:6)'], 'out': ['F(1:5)']}, \
            'advance': {'in': ['U(1:5)', 'F(1:5)', 'D(1)', 'D(2:4)', 'D(5)'], 'out': ['Unew(1:5)']}, \
#            'fill_boundary': {'in': ['Unew(1:5)'], 'out': ['Unew4(1:5)']}, \
           }, \
    'subgraphs': [{'name': 'cluster_updateU_1', 'nbunch': ['diff_2', 'diff_3', 'diff_4']}, \
                  {'name': 'cluster_updateU_2', \
                   'nbunch': ['diff_1', 'diff_5', 'diff_6', 'diff_7', 'diff_8', \
                              'hyp_1', 'hyp_2', 'hyp_3', 'advance', \
                              'D(1)', 'D(2:4)', 'D(5)', 'F(1:5)']}], \
  }, \
  {'deps': {'ctoprim'  : {'in': ['*U(1:5)'], 'out': ['Q(2:6)'] }, \
            'updateU_1': {'in': ['*Q(2:6)'], 'out': ['Q\'(1:9)']}, \
            'updateU_2': {'in': ['*U(1:5)', '*Q(2:6)', '*Q\'(1:9)'], 'out': ['Unew(1:5)']}, \
           }, \
   'subgraphs': [], \
  }, \
  {'deps': {'updateU': {'in': ['*U(1:5)'], 'out': ['Unew(1:5)']}, \
           }, \
   'subgraphs': [], \
  }]

funcNodeProps = {'shape': 'oval', 'style': 'filled', 'fillcolor': 'cornsilk'}
dataNodeProps = {'shape': 'box', 'style': 'filled', 'fillcolor': 'cadetblue1'}
edgeProps = {}
# ghostEdgeProps = {'style': 'bold', 'color': 'red'}
ghostEdgeProps = {'style': 'dashed', 'color': 'black'}
subGraphProps = {'shape': 'box', 'style': 'filled', 'color': 'grey'}

# duplicate keys in d2 overwrite those in d1
def mergeDicts(d1, d2):
  d3 = {}
  for (k,v) in d1.iteritems():
    d3[k] = v
  for (k,v) in d2.iteritems():
    d3[k] = v
  return d3

suppressLabels = True
suppressSubgraphs = True

for (i, graph) in enumerate(graphs, start = 1):
  G = pgv.AGraph(name = '', directed = True, nodesep=.1)
  deps = graph['deps']
  for (fname, fval) in deps.iteritems():
    print fname
    if suppressLabels:
      fval['label'] = ''
    if 'label' in fval:
      G.add_node(fname, label=fval['label'], **funcNodeProps)
    else:
      G.add_node(fname, **funcNodeProps)
    for data in fval['in']:
      p = edgeProps
      if data[0] == "*":
        data = data[1:]
        p = ghostEdgeProps
      dataLabel = data
      if suppressLabels:
        dataLabel = ''
      G.add_node(data, label=dataLabel, **dataNodeProps)
      G.add_edge(data, fname, **p)
    for data in fval['out']:
      dataLabel = data
      if suppressLabels:
        dataLabel = ''
      G.add_node(data, label=dataLabel, **dataNodeProps)
      G.add_edge(fname, data, **edgeProps)
  if not suppressSubgraphs:
    for subgraph in graph['subgraphs']:
      Gsub = G.add_subgraph(**mergeDicts(subGraphProps, subgraph))
  G.layout(prog = 'dot')
  G.draw('depgraph%d.png' % i)
