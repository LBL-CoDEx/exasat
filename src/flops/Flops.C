/*
  Didem Unat
  
 */

#include <rose.h>


#include "Flops.h"
#include "../constants.h"
#include "../utils/ExaSatOptions.h"
#include "../utils/OutputFormatting.h"

using namespace std;
using namespace SageInterface;


bool Flops::isDoubleOrFloat(SgExpression *exp)
{
  Rose_STL_Container<SgNode*> vars = NodeQuery::querySubTree(exp, V_SgVarRefExp); 
  Rose_STL_Container<SgNode*>::iterator v = vars.begin();
  
  for( ; v != vars.end(); v++)
    {
      SgVarRefExp* var = isSgVarRefExp(*v);

      SgInitializedName* iname = convertRefToInitializedName(var);
      ROSE_ASSERT(iname);

      SgType* type = iname->get_type();
      ROSE_ASSERT(type);

      if(isSgTypeDouble(type) || isSgTypeFloat(type))
	return true;
      
      SgType* element_type = type->findBaseType();

      if(element_type != NULL)
	if(isSgTypeDouble(element_type) || isSgTypeFloat(element_type))
	  return true;
    }
  return false;
}
bool Flops::isValidFlop(SgBinaryOp* binOp)
{
  //TODO:
  //Returns false if operands are integers for ex. index exp. 
  //Returning true doesn't mean it is a valid flop. Need to add more tests

  //Get the lhs and rhs of the bindary op
  //if lhs or rhs is a VarRefExp and it's type is integer, then return false
  SgExpression* lhs = binOp->get_lhs_operand();
  SgExpression* rhs = binOp->get_rhs_operand();

  ROSE_ASSERT(lhs);
  ROSE_ASSERT(rhs);

  if(isDoubleOrFloat(lhs))
    return true;

  if(isDoubleOrFloat(rhs))
    return true;
  /*
  if(isSgVarRefExp(lhs))
    {
      SgVarRefExp* varExp = isSgVarRefExp(lhs);
      SgType* type = varExp->get_type();
      if(isSgTypeInt(type))
	return false;
    }
  if(isSgVarRefExp(rhs))
    {
      SgVarRefExp* varExp = isSgVarRefExp(rhs);
      SgType* type = varExp->get_type();
      if(isSgTypeInt(type))
	return false;
    }  
  if(isSgIntVal(rhs) && isSgIntVal(lhs))
    return false;

  if(isSgIntVal(rhs) || isSgIntVal(lhs))
    {
      
      return false;
    }

  return true;
  */
  return false;
}


int Flops::getNumberOfMathFunctions(SgNode* node)
{
  Rose_STL_Container<SgNode*> funcCalls = NodeQuery::querySubTree(node, V_SgFunctionRefExp); 
  Rose_STL_Container<SgNode*>::iterator funcItr = funcCalls.begin();
  
  int count=0; 
  for( ; funcItr != funcCalls.end(); funcItr++)
    {
      SgFunctionRefExp* exp = isSgFunctionRefExp(*funcItr);
      SgFunctionSymbol* sym = exp->get_symbol();
      ROSE_ASSERT(sym);
      
      string func_name = sym->get_name().str();

      if(math_func_list.find(func_name) != math_func_list.end())
	{
	  count++;
	}
    }
  return count; 
}

int Flops::computeMulFlops(SgNode* node)
{

  Rose_STL_Container<SgNode*> binOps = NodeQuery::querySubTree(node, V_SgBinaryOp); 
  Rose_STL_Container<SgNode*>::iterator opItr = binOps.begin();

  int count_mul = 0;
  for( ; opItr != binOps.end(); opItr++)
    {
      SgBinaryOp* binOp = isSgBinaryOp(*opItr);
      
      if(isValidFlop(binOp))
	{
	  switch(binOp -> variantT())
	    {
	    case V_SgExponentiationOp:
	      {
		count_mul++;  //check if the exponent is 2, then count that as multiply 
		break;
	      }
	    case V_SgMultiplyOp:
	      {
		count_mul++;
		break;
	      }
	    default:
	      break;
	    }
	}
    }
  return count_mul; 
}

int Flops::computeDivFlops(SgNode* node)
{

  Rose_STL_Container<SgNode*> binOps = NodeQuery::querySubTree(node, V_SgBinaryOp); 
  Rose_STL_Container<SgNode*>::iterator opItr = binOps.begin();

  int count_div = 0;
  for( ; opItr != binOps.end(); opItr++)
    {
      SgBinaryOp* binOp = isSgBinaryOp(*opItr);
      
      if(isValidFlop(binOp))
	{
	  switch(binOp -> variantT())
	    {
	    case V_SgDivideOp:
	      {
		count_div++;
		break;
	      }
	    default:
	      break;
	    }
	}
    }
  return count_div; 
}


int Flops::computeAddFlops(SgNode* node)
{

  Rose_STL_Container<SgNode*> binOps = NodeQuery::querySubTree(node, V_SgBinaryOp); 
  Rose_STL_Container<SgNode*>::iterator opItr = binOps.begin();

  int count_add = 0;
  for( ; opItr != binOps.end(); opItr++)
    {
      SgBinaryOp* binOp = isSgBinaryOp(*opItr);
      
      if(isValidFlop(binOp))
	{
	  switch(binOp -> variantT())
	    {
	    case V_SgAddOp:
	      {
		count_add++;
		break;
	      }
	    case V_SgSubtractOp:
	      {
		count_add++;
		break;
	      }
	    default:
	      break;
	    }
	}
    }
  return count_add; 
}



// eof


