/*
  by Didem Unat 
  */

#ifndef _FUNCCALLHANDLER_H
#define _FUNCCALLHANDLER_H

#include <rose.h>
#include <iostream>
#include <string>
#include <sstream>

#include "../handler/NodeHandler.h"
#include "../datatypes.h"
#include "ASTtools.hh"


using namespace std;
using namespace SageInterface;

class FunctionCallHandler : public NodeHandler
{
  
 private:
  //line number and function information
  FuncCallInfo_t fCallInfo; 

  bool matchArgumentsWithParameters(SgExpressionPtrList callArgList, 
				    SgInitializedNamePtrList funcParamList);

 public:

  void process();
  void dump();
  void closeXMLTag();
  
  FunctionCallHandler(SgNode* node);
  ~FunctionCallHandler();

  static bool isFortranMathFunction(SgNode* node);
  static bool isIgnoredCFunction(SgNode* node);
};

#endif //_FUNCCALLHANDLER_H
