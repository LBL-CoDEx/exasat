/*
   by Didem Unat 
  */

#ifndef _WHOLEARRAYEXP_H
#define _WHOLEARRAYEXP_H

//#include <rose.h>
#include <iostream>
#include <string>
#include <sstream>
#include "../datatypes.h"
#include "VarSym.hh" 
#include "ASTtools.hh"  

using namespace std;
using namespace SageInterface;

class WholeArrayExp
{
 public:

  void convertWholeArraysToFortranDos(SgNode* node);

  //selected statements by the isWholeArrayStatement become a nested do loop
  void convertWholeArrayStatementToFortranDoStatement(SgStatement* stmt);

  //checks if a statement contains any whole array references, call isWholeArrayExpression
  static bool isWholeArrayStatement(SgStatement* s, int* depth);

 private:
  //check if the expression is a whole array expression A or A(:), computes the depth of the loop
  static bool isWholeArrayExpression(SgExpression* ref, int* depth);
  
  
  //gets the lower and upper bounds (range of the array that is accessed in the whole array statement
  //these bounds are calculated per array 
  void getLoopBounds(SgExpression* array, vector<SgExpression*>& lower, 
		     vector<SgExpression*>& upper, vector<SgExpression*>& stride);

  //get the subscript expression in the alloc statement
  SgExpression* getAllocExpression(SgExpression* array);

  //tries to get the bounds from the variable declaration and alloc statement
  void parseSubscript(SgExpression* sub, vector<SgExpression*>& lower, 
		      vector<SgExpression*>& upper, vector<SgExpression*>& stride, int index);
  //for whole array operations
  void buildFortranDoLoop(SgBasicBlock* body, RangeInfo_t* lInfo);

  SgVariableDeclaration* buildIndexVariable(SgBasicBlock* func_body);

};

#endif //_WHOLEARRAYEXP_H
