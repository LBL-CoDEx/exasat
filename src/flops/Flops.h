/*
  by Didem Unat 
  */

#ifndef _NOTAFLOPS_H
#define _NOTAFLOPS_H

//#include <rose.h>
#include <iostream>
#include <string>
#include <sstream>

using namespace std;
using namespace SageInterface;

class Flops
{
public:

  static bool isDoubleOrFloat(SgExpression* exp);
  static int computeAddFlops(SgNode* node);
  static int computeMulFlops(SgNode* node);
  static int computeDivFlops(SgNode* node);
  static bool isValidFlop(SgBinaryOp* binOp);
  static int getNumberOfMathFunctions(SgNode* node);
};

#endif //_NOTAFLOPS_H
