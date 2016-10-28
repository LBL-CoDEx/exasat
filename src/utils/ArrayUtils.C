/*
  Didem Unat
 */

#include <rose.h>

#include "ArrayUtils.h"
#include "Utils.h"
#include "stringify.h"
#include "../datatypes.h"

using namespace std;
using namespace SageInterface;
using namespace Utils; 
using namespace SageBuilder;



int ArrayUtils::getArrayDimension(SgPntrArrRefExp* exp)
{
  //ROSE cannot parse the Fortran arr refs 
  //(i,j,k) and (i) has the same dim, we count the commas which is not ellegant 
  string str = exp->unparseToString();

  vector<string> vec_str; 

  Utils::stringSplit(str, ",", vec_str);

  return vec_str.size();
}


int ArrayUtils::get_array_rank(SgExpression* exp )
{
  SgInitializedName* name = convertRefToInitializedName(exp);
  ROSE_ASSERT(name); 
  
  SgType* type = name->get_type();
  
  if ( !isSgArrayType(type) ) // check if the array is part of a named type struct, class or union
    type = getTypeInNamedType(name, exp); //including class, union or struct
  
  if(type != NULL &&  isSgArrayType(type)) 
    {
      SgExprListExp* expList = isSgArrayType(type)->get_dim_info();
      SgExpressionPtrList& list = expList->get_expressions();	
      
      return list.size();
    }
  return 0;
}

SgType* ArrayUtils::getTypeInNamedType(SgInitializedName* name, SgExpression* exp)
{
  SgType* type = name->get_type();
  ROSE_ASSERT(type);
  SgType* var_type = NULL;
   
  if(isSgNamedType(type))
    {
      SgNamedType* namedType= isSgNamedType(type);
      
      SgDeclarationStatement* named_decl = namedType->get_declaration();

      if(named_decl != NULL) {
	// the case of class declaration, including struct, union
	SgClassDeclaration* classDeclaration = isSgClassDeclaration(named_decl);
	
	if (classDeclaration != NULL)
	  {
	    SgClassDeclaration* defining_declaration = isSgClassDeclaration(classDeclaration -> get_definingDeclaration());
	    ROSE_ASSERT(defining_declaration); 
	    
	    SgClassDefinition* def = isSgClassDefinition(defining_declaration -> get_definition());
	    ROSE_ASSERT(def) ; 
	    
	    vector<SgDeclarationStatement *> decls = def -> get_members();
	    for (vector<SgDeclarationStatement*>::iterator d = decls.begin(); d != decls.end(); d++)
	      
	      if(isSgVariableDeclaration(*d))
		{
		  vector<SgInitializedName *> vars = isSgVariableDeclaration(*d) -> get_variables();
		  
		  for (int k = 0; k < vars.size(); k++) {
		    
		    SgInitializedName *var = vars[k];
		    ROSE_ASSERT(var); 
		    
		    if (exp-> unparseToString () == var-> get_name().str())
		      {
			SgType* var_type = var -> get_type();
			
			return var_type; 
		      }
		  }
		}
	  }
      }//end of class declaration   
    }// end of Named Type 
  return var_type; 
}

bool ArrayUtils::isFullArrayReference(SgPntrArrRefExp* exp)
{
  return !isSgPntrArrRefExp( exp -> get_parent());
}
void ArrayUtils::getOnlyArraySymbols(ASTtools::VarSymSet_t& V,
				     ASTtools::VarSymSet_t& newSym)
{
  for(ASTtools::VarSymSet_t::iterator i= V.begin() ; i!= V.end() ; ++i)
    {
      const SgVariableSymbol* sym = (*i);
      SgType* type = sym ->get_type();

      if(isArrayOrPointerType(type))
	newSym.insert(sym);
    }
}

void ArrayUtils::collectLocalVarSyms (SgStatement* stmt, 
				      ASTtools::VarSymSet_t& locals)
{

  SgFunctionDeclaration* func_decl = getEnclosingFunctionDeclaration(isSgNode(stmt));
  ROSE_ASSERT(func_decl);

  const SgInitializedNamePtrList arg_list = func_decl->get_args();
 
  ASTtools::VarSymSet_t syms;
  ASTtools::collectDefdVarSyms (stmt, syms);
  
  //not all the defined variables are local. 
  //in Fortran, variables with intent modifier are not local 

  for(ASTtools::VarSymSet_t::iterator s = syms.begin() ; s != syms.end(); s++)
    {
      const SgVariableSymbol* var_sym = (*s);
      ROSE_ASSERT(var_sym);
      SgInitializedName* iname = var_sym->get_declaration();
      ROSE_ASSERT(iname);

      SgDeclarationStatement* var_decl = iname ->get_declaration();;  
      ROSE_ASSERT(var_decl);

      //Fortran variables can have intent modifier, we skip those because they are not local
      if(var_decl -> get_declarationModifier().get_typeModifier().isIntent_in())
	{
	  continue; 
	}
      if(var_decl -> get_declarationModifier().get_typeModifier().isIntent_out())
	{
	  continue; 
	}
      if(var_decl -> get_declarationModifier().get_typeModifier().isIntent_inout())
	{
	  continue; 
	}

      //in Fortran, we need to exclude the function parameters 
      //because they are declared again in the function body      
      if(find(arg_list.begin(), arg_list.end(), iname) != arg_list.end())
	{
	  continue;
	}
      locals.insert((var_sym));
    }
}

bool ArrayUtils::isArrayOrPointerType(SgType* type)
{
  //we are interested in only array types
  if(isSgArrayType(type) || isSgPointerType(type))
    {
      return true;
    }
  else if(isSgTypedefType(type))
    {
      
      SgType* base_type = isSgTypedefType(type)->get_base_type();

      return isArrayOrPointerType(base_type);
    }
  return false;
}


bool ArrayUtils::isPntrType(SgType* type)
{
  //we are interested in only array types
  if(isSgPointerType(type))
    {
      return true;
    }
  else if(isSgTypedefType(type))
    {
      
      SgType* base_type = isSgTypedefType(type)->get_base_type();

      return ArrayUtils::isPntrType(base_type);
    }
  return false;
  }

bool ArrayUtils::isArrayReference(SgExpression* ref,
				   SgExpression** arrayName/*=NULL*/,
				   vector<SgExpression*>* subscripts/*=NULL*/)
{
  SgExpression* arrayRef=NULL;

  if(ref->variantT() == V_SgPntrArrRefExp)
    {
      SgPntrArrRefExp* arr = isSgPntrArrRefExp(ref);

      *arrayName = arr->get_lhs_operand();
      
      SgExpression* subs = arr->get_rhs_operand(); 

      if(isSgExprListExp(subs))
	{
	  SgExprListExp* expList = isSgExprListExp(subs);  
	  //For fortran, we need to split the expression into ExprList 
	  SgExpressionPtrList& list = expList->get_expressions();
	  SgExpressionPtrList::iterator i;
	  for (i = list.begin(); i != list.end(); ++i) {
	    SgExpression* expr = *i;
	    subscripts->push_back(expr);
	  }
	}
      else 
	subscripts->push_back(subs);

      return true;
    }
  return false; 
}


// eof


