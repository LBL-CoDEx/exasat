#include "AnalysisInterface.h"
#include <iostream>

using namespace std;

class FunctionSideEffectAnalysis: public FunctionSideEffectInterface 
{
public: 

  FunctionSideEffectAnalysis() {};

  bool get_modify(AstInterface& fa, const AstNodePtr& fc,
		  CollectObject<AstNodePtr>* collect = 0);
  bool get_read(AstInterface& fa, const AstNodePtr& fc,
		CollectObject<AstNodePtr>* collect = 0) ;
  
};
/*
bool OperatorSideEffectAnnotation::
get_modify( AstInterface& fa, const AstNodePtr& fc,
	    CollectObject< AstNodePtr >* collect)
{
  AstInterface::AstNodeList args;
  OperatorSideEffectDescriptor mod;
  if (modInfo.known_operator( fa, fc, &args, &mod)) {
    if (collect != 0)
      mod.get_side_effect(fa, args, *collect);
    return true;
  }
  return false;
}

bool OperatorSideEffectAnnotation::
get_read( AstInterface& fa, const AstNodePtr& fc,
	  CollectObject< AstNodePtr >* collect)
{
  AstInterface::AstNodeList args;
  OperatorSideEffectDescriptor read;
  if  (readInfo.known_operator( fa, fc, &args, &read)) {
    if (collect != 0)
      read.get_side_effect(fa, args, *collect);
    return true;
  }
  return false;
}


bool ArrayUseAccessFunction::get_modify(AstInterface& fa, const AstNodePtr& fc,
					CollectObject<AstNodePtr>* collect)
{
  AstInterface::AstNodeList args;
  if (prev1 != 0 && prev1->get_modify(fa, fc, collect))
    return true;
  std::string sig;
  AstNodePtr f;
  if (fa.IsFunctionCall(fc, &f,&args) && fa.IsVarRef(f,0,&sig) && sig == funcname) {
    return true;
  }
  return false;
}

bool ArrayUseAccessFunction::get_read(AstInterface& fa, const AstNodePtr& fc,
				      CollectObject<AstNodePtr>* collect)
{
  AstInterface::AstNodeList args;
  if (prev1 != 0 && prev1->get_read(fa, fc, collect))
    return true;
  std::string sig;
  AstNodePtr f;
  if (fa.IsFunctionCall(fc, &f,&args) && fa.IsVarRef(f,0,&sig) && sig == funcname) {
    if (collect != 0)  {
      AstInterface::AstNodeList::const_iterator argp = args.begin();
      for ( ++argp; argp != args.end(); ++argp) {
	(*collect)(*argp);
      }
    }
    return true;
  }
  return false;
}

bool ArrayAnnotation :: get_modify(AstInterface& _fa, const AstNodePtr& fc,
	   CollectObject<AstNodePtr>* collect)
{
  CPPAstInterface& fa = static_cast<CPPAstInterface&>(_fa);
  if ( is_access_array_elem(fa, fc) || is_access_array_length(fa, fc))
    return true;
  AstNodePtr array;
  if (is_reshape_array( fa,fc, &array)) {
    if (collect != 0)
      (*collect)(array);
    return true;
  }
  return OperatorSideEffectAnnotation::get_inst()->get_modify(fa, fc, collect);
}

bool ArrayAnnotation :: get_read(AstInterface& _fa, const AstNodePtr& fc, CollectObject<AstNodePtr>* collect)
{
  CPPAstInterface& fa = static_cast<CPPAstInterface&>(_fa);
  AstNodePtr dim;
  if (is_access_array_length( fa, fc, 0, &dim)) {
    if (collect != 0)
      (*collect)(dim);
    return true;
  }
  CPPAstInterface::AstNodeList args;
  if (is_access_array_elem( fa, fc, 0, &args)) {
    if (collect != 0) {
      for (CPPAstInterface::AstNodeList::iterator p = args.begin();
	   p != args.end(); ++p)
	(*collect)(*p);
    }
    return true;
  }
  return OperatorSideEffectAnnotation::get_inst()->get_read(fa, fc, collect);
}
*/
