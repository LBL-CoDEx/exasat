
#include "MintArrayInterface.h"


using namespace std;
using namespace SageInterface;
using namespace SageBuilder;
//! Check if an expression is an array access. If so, return its name and subscripts if requested.
// Based on AstInterface::IsArrayAccess()                                                                                                 

bool MintArrayInterface::isArrayReference(SgExpression* ref, 
					  SgExpression** arrayName/*=NULL*/, 
					  vector<SgExpression*>* subscripts/*=NULL*/)
{
  SgExpression* arrayRef=NULL;

  if (ref->variantT() == V_SgPntrArrRefExp) {
    if (subscripts != 0 || arrayName != 0) {
      SgExpression* n = ref;
      while (true) {
	SgPntrArrRefExp *arr = isSgPntrArrRefExp(n);
	if (arr == 0)
	  break;
	n = arr->get_lhs_operand();
	// store left hand for possible reference exp to array variable                                                      
	if (arrayName!= 0)
	  arrayRef = n;
	// right hand stores subscripts                                                                                      
	if (subscripts != 0){
	  subscripts->push_back(arr->get_rhs_operand());
	  //cout << "sub: " << (arr->get_rhs_operand())->unparseToString() << endl;
	}
      } // end while                                       
      if  (arrayName !=NULL)
        {
          *arrayName = arrayRef;
        }
    }
    return true;
  }
  return false;
}

void MintArrayInterface::linearizeArrays(SgSourceFile* file)
{
 
  Rose_STL_Container<SgNode*> kernelList = NodeQuery::querySubTree(file, V_SgFunctionDefinition);
  Rose_STL_Container<SgNode*>::iterator kernel ;

  for (kernel = kernelList.begin() ; kernel != kernelList.end();  kernel++)
    {      
      SgFunctionDefinition* kernel_def = isSgFunctionDefinition(*kernel);

      SgFunctionDeclaration* kernel_decl = kernel_def -> get_declaration();

      string func_name = kernel_decl->get_name().getString() ;
    
      if(kernel_decl ->get_functionModifier().isCudaKernel())
	{
	  linearizeArrays(kernel_decl);
	}
    }

}


void MintArrayInterface::linearizeArrays(SgFunctionDeclaration* kernel)
{

  //replaces all the array references in the basic block into 1-dim array references
  //do not change the array name?

  //Step 1: Find all the pointer array reference expressions like E[i][j] 
  ROSE_ASSERT(kernel);

  SgBasicBlock* kernel_body = kernel->get_definition()->get_body();

  Rose_STL_Container<SgNode*> arrList = NodeQuery::querySubTree(kernel_body, V_SgPntrArrRefExp);//V_SgExpression);

  Rose_STL_Container<SgNode*>::iterator arr;

  std::map<SgVariableSymbol*, string> sizeList; 

  for(arr = arrList.begin(); arr != arrList.end(); arr++)
    {
      SgExpression* exp = isSgExpression(*arr);
	
	if(isSgPntrArrRefExp(exp))
	  {
	    //cout << exp->unparseToString() << endl;

	    //Step2 : Find array reference name and index expressions
	    SgExpression* arrayName; 
	    vector<SgExpression*>  subscripts; //first index is i if E[j][i]
	    
	    //index list are from right to left 
	    bool yes = MintArrayInterface::isArrayReference(isSgExpression(*arr), &arrayName, &subscripts);
	    assert(yes == true);
	    
	    //do not linearize the shared memory variables 
	    //TODO:med use the shared memory type to differentiate it 
	    //rather than its name. if the type is __shared__ then exclude that ref
	    
	    SgVarRefExp* arrRef = isSgVarRefExp(deepCopyNode(arrayName));

	    SgVariableSymbol* sym= isSgVarRefExp(arrRef)->get_symbol();

	    ROSE_ASSERT(sym);

	    string arr_str = sym->get_name().str() ; 
	    
	    //step 3: if the array dimension is greater than 1, then convert the reference to 1-dim array
	    if(subscripts.size() > 1  && arr_str.find("_sh_block") == -1 ) {
	      
	      //step 4: compute the new index expression from the original index expressions
	      //For example: E[j][i] becomes j*n + i, E[j][i][k] becomes j*n*m+ i*m + k
	      //might be better if I use reverse iterator 
	      
	      string sizeStr = arr_str ; 
	      sizeStr = "size" + sizeStr ;

	      if(sizeList.find(sym) == sizeList.end()) // this is the first time we encounter this arr
		  sizeList[sym] = sizeStr ; 

	      SgVarRefExp* sizeExp = buildVarRefExp(sizeStr, kernel_body); 

	      vector<SgExpression*>::iterator index= subscripts.begin();
	      SgExpression* indexExp = deepCopy(*index);
	      index++;
	      
	      for(; index != subscripts.end(); index++)
		{
		  indexExp = buildAddOp(indexExp, buildMultiplyOp( deepCopy(*index), sizeExp)); //buildMultiplyOp(*index, size));
		}
	      
	      SgExpression* newExp = buildPntrArrRefExp(arrRef, indexExp);

	      replaceExpression(exp, newExp);
	    }
	    //skip the sub  array references like E[j] in E[j][i]
	    
	    int arrayDim = subscripts.size() ;
	    
	    arr = arr + arrayDim - 1;
	  }
    }


  std::map<SgVariableSymbol*, string>::iterator it ;

  for(it = sizeList. begin(); it != sizeList.end(); it++)
    {
      SgVariableSymbol* sym = it->first;

      string arrStr = sym->get_name();
      string sizeStr = it->second;

      SgVarRefExp* arrRef = buildVarRefExp(arrStr, kernel_body);

      SgExpression* glmem1 = buildAddressOfOp(buildPntrArrRefExp( arrRef, buildVarRefExp("0", kernel_body)));
      SgExpression* glmem2 = buildAddressOfOp(buildPntrArrRefExp( arrRef, buildVarRefExp("1", kernel_body)));
      
      SgAssignInitializer* init=buildAssignInitializer(buildSubtractOp(glmem2, glmem1));
      
      SgVariableDeclaration* size_decl = buildVariableDeclaration(sizeStr , buildIntType(), init, kernel_body); 
      
      SgStatement* stmt = SageInterface::getFirstStatement(kernel_body);
      
      insertStatementAfter(stmt,size_decl);      
    }
}
