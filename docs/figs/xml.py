#!/usr/bin/python

import pygraphviz as pgv

graph = {'program'  : ['function'], \
         'function' : ['local', 'nonlocal', 'loop', 'comm', 'functioncall'], \
         'loop'     : ['scalar', 'array', 'pininfo', 'loop'], \
         'array'    : ['access'], \
         'comm'     : ['commscalar', 'commarray'], \
        }

nodeProps = {'shape': 'box', 'style': 'filled', 'fillcolor': 'cornsilk'}
edgeProps = {}

G = pgv.AGraph(name = '', directed = True)
for element in graph.keys():
  G.add_node(element, **nodeProps)
  for child in graph[element]:
    G.add_node(child, **nodeProps)
    G.add_edge(element, child, **edgeProps)
G.layout(prog = 'dot')
G.draw('xml.png')
