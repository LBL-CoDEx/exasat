/*
  Didem Unat
  
 */

#include <rose.h>


#include "./Scalars.h"
#include "../utils/ArrayUtils.h"
#include "../memory/AccessAnalysis.h"

using namespace std;
using namespace SageInterface;
using namespace ArrayUtils;


void Scalars::getActiveScalarVars(SgNode* node, 
				  std::set<SgVarRefExp*>& scalars)
{
  Rose_STL_Container<SgNode*> varList= NodeQuery::querySubTree(node, V_SgVarRefExp);
  Rose_STL_Container<SgNode*>::iterator v = varList.begin();

  for(; v!= varList.end() ; v++)
    {
      SgVarRefExp* varExp = isSgVarRefExp(*v);
if(!varExp)exit(0);
      ROSE_ASSERT(varExp);
      
      //this function returns the struct not the field if the field is a scalar 
      //it returns the array if the field is an array 
      //SgInitializedName* iname = convertRefToInitializedName(varExp);
      //ROSE_ASSERT(iname);

      //instead get the declaration from the symbol 
      //this is safe to call this because we queried SgVarRefExp 
      SgVariableSymbol* symbol = varExp-> get_symbol(); 

      if(symbol != NULL){
	SgInitializedName* name = symbol->get_declaration();
	ROSE_ASSERT(name);

	SgType* type = name->get_type();
	ROSE_ASSERT(type);

	//arrays are not scalar
	if(isArrayOrPointerType(type))
	  continue; 
	
	//might be struct, need to see if the field is a scalar
	if(isSgNamedType(type)) 
	  {
	    //this function returns the type of the field in a struct (or union or class) 
	    SgType* typeInStruct = getTypeInNamedType(name, varExp);  
	    
	    if(typeInStruct == NULL) //accessing the struct not to its field, can be argument to a function
	      continue; 
	    else { //this is an access to the field of a struct
	      
	      //it looks like if the VarRefExp is an array and part of a struct, its iname is the field not the struct
	      //it is not true for scalars, iname of a sclara which is part of a struct returns its struct not the field 
	      //thus the following check is unnecessary 
	      if(isArrayOrPointerType(typeInStruct)) //if it is array, skip it 
		continue;
	      if(isSgNamedType(typeInStruct)) //a struct nested in other struct, skip this too A%b%x
		continue;             //actual field should be caught in another VarRefExp
	    }
	  }// end of SgNamedType
	scalars.insert(varExp);

      }//symbol is not NULL
    }
}

void Scalars::getScalarReadWriteCounts(ScalarRWMap_t& scalarRW, 
				       FirstAccessMap_t& firstAccess,
				       std::set<SgVarRefExp*>& diffScalars)
{
  AccessAnalysis* access = new AccessAnalysis(); 
  for(std::set<SgVarRefExp*>::iterator var= diffScalars.begin(); var!= diffScalars.end(); var++)
    {
      SgVarRefExp* varExp= *var;
      ROSE_ASSERT(varExp);
      
      bool read = access->isRead(isSgExpression(varExp));
      
      SgVariableSymbol* symbol = varExp-> get_symbol(); 
      if(symbol != NULL){

	SgInitializedName* iname = symbol->get_declaration();
	ROSE_ASSERT(iname);

	//if varExp is a field of a struct, when we convert it to initialized name, it will return 
	// the struct name not the field name, instead get the declaration from the symbol 
	//SgInitializedName* iname = convertRefToInitializedName(varExp);
	//ROSE_ASSERT(iname);

	//scalars do not have a component, so create an empty component field for scalars
	access->firstWrittenOrReadAnalysis(firstAccess, isSgExpression(varExp), make_pair(iname, ""));

	ScalarRWMap_t::iterator rw = scalarRW.find(iname);
	if(rw != scalarRW.end()) //found
	  {
	    AccessPattern_t* rwacc = &(rw->second); 
	    if(read) rwacc->readCount++;
	    else rwacc->writeCount++;
	  }
	else //first appearance 
	  {
	    AccessPattern_t rwacc;
	    rwacc.readCount = (read)? 1 : 0;
	    rwacc.writeCount = (read)? 0 : 1; 
	    scalarRW.insert(make_pair(iname, rwacc));
	  }
      } // end of symbol
    }
  delete access; 
}

void Scalars::getScalarInames(std::set<SgInitializedName*>& scalarIname, 
			      std::set<SgVarRefExp*>& diffScalars)
{

  for(std::set<SgVarRefExp*>::iterator var= diffScalars.begin(); var!= diffScalars.end(); var++)
    {
      SgVarRefExp* varExp= *var;
      ROSE_ASSERT(varExp);
            
      //if varExp is a field of a struct, when we convert it to initialized name, it will return 
      // the struct name not the field name, instead get the declaration from the symbol 
      //SgInitializedName* iname = convertRefToInitializedName(varExp);
      //ROSE_ASSERT(iname);

      SgVariableSymbol* symbol = varExp-> get_symbol(); 
      if(symbol != NULL){
	SgInitializedName* iname = symbol->get_declaration();
	ROSE_ASSERT(iname);
      
	scalarIname.insert(iname);
      }
    }
}

// eof


