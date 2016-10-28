/*
   by Didem Unat 
  */

#ifndef _ARRAYUTILS_H
#define _ARRAYUTILS_H

//#include <rose.h>
#include <iostream>
#include <string>
#include <sstream>
#include "../datatypes.h"
#include "VarSym.hh" 
#include "ASTtools.hh"  

using namespace std;
using namespace SageInterface;

namespace ArrayUtils
{
  //returns the type of an array if the array is part of a named type such as class, union or struct
  //type info for an array includes its rank and element type
  SgType* getTypeInNamedType(SgInitializedName* name, SgExpression* exp);

  //get the subscript expression in the alloc statement
  SgExpression* getAllocExpression(SgExpression* array);

  int get_array_rank(SgExpression* array);

  void convertSubscriptExpIntoString(vector<SgExpression*> subs, 
				     vector<string>& subsInStr);

  //int getLoopDepth(SgStatement* );

  int getArrayDimension(SgPntrArrRefExp* exp);

  bool isFullArrayReference(SgPntrArrRefExp* exp);

  bool isArrayOrPointerType(SgType* type);

  bool isPntrType(SgType* type);

  bool isArrayReference(SgExpression* ref,
			SgExpression** arrayName/*=NULL*/,
			vector<SgExpression*>* subscripts/*=NULL*/);

  void collectLocalVarSyms (SgStatement* s, ASTtools::VarSymSet_t& syms);
  
  void getOnlyArraySymbols(ASTtools::VarSymSet_t& V, 
			   ASTtools::VarSymSet_t& newSym);
};

#endif //_ARRAYUTILS_H
