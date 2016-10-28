#include "AnalysisInterface.h"

#include "./FuncSideEffect.h"
#include <iostream>
#include "../constants.h"

using namespace std;

bool FunctionSideEffectAnalysis::get_modify(AstInterface& fa, const AstNodePtr& fc,
					   CollectObject<AstNodePtr>* collect ) 
{
  AstInterface::AstNodeList args;
 
  std::string func_name;
  AstNodePtr f;
  if (fa.IsFunctionCall(fc, &f,&args) && fa.IsVarRef(f,0,&func_name)) {

    if(math_func_list.find(func_name) != math_func_list.end())
      return true; //not modified, read-only
    if(math_func_low_lat_list.find(func_name) != math_func_low_lat_list.end())
      {
	return true; //read-only 
      }
    if(fortran_intrinsic_func_list.find(func_name) != fortran_intrinsic_func_list.end())
      return true; //read-only
  }
  return false;
}

bool FunctionSideEffectAnalysis::get_read(AstInterface& fa, const AstNodePtr& fc,
					  CollectObject<AstNodePtr>* collect ) 
{
  //arguments to a function are assumed to be read
  return false; 
  /* 
     AstInterface::AstNodeList args; 
     std::string sig; 
     AstNodePtr f; 
     
     if (fa.IsFunctionCall(fc, &f,&args) && fa.IsVarRef(f,0,&sig)) {
     if (collect != 0)  {
     AstInterface::AstNodeList::const_iterator argp = args.begin();
     for ( ++argp; argp != args.end(); ++argp) {
     (*collect)(*argp);
     }
     }
     return false;
     }
  return true;
  */
}
