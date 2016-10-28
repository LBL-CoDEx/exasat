/*
  by Didem Unat 
  */

#ifndef _NODEHANDLER_H
#define _NODEHANDLER_H

#include <rose.h>
#include <iostream>
#include <string>
#include <sstream>


using namespace std;


class NodeHandler
{
 protected:

  NodeHandler* parent;
  vector<NodeHandler*> children;
  SgNode* myNode;
  bool visited; 

 public:

  //gathering data
  virtual void process() = 0;

  //Output related stuff
  virtual void dump() = 0;
  virtual void closeXMLTag () = 0; 
  
  virtual void removeDuplicates() {}

  //parent - child relationships
  void setParent(NodeHandler* p )     { parent = p;}
  NodeHandler* getParent()            { return parent;}
  void addChild(NodeHandler* c)       { children.push_back(c) ;}
  vector<NodeHandler*>* getChildren() { return &children;}

  //node traversal
  SgNode* getNode()            { return myNode; }
  void setNode(SgNode* n)      { myNode = n; } 
  void setVisited(bool v=true) { visited = v;}
  bool getVisited()            { return visited;}

};

#endif //_NODEHANDLER_H
