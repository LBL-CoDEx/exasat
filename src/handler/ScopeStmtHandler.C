
#include <rose.h>
#include "./ScopeStmtHandler.h"
#include "./NodeHandler.h"

#include "../utils/OutputFormatting.h"
#include "../utils/Utils.h"

#include "../flops/Flops.h"
#include "../memory/ReadWriteProp.h"
#include "../memory/AccessAnalysis.h"


using namespace std;
using namespace SageInterface;


void ScopeStmtHandler::computeFlops(SgNode* node)
{
  //count the flops in this loop
  int add = Flops::computeAddFlops(node);
  int multiply = Flops::computeMulFlops(node);
  int division = Flops::computeDivFlops(node);
  //count the math intrinsic functions 
  int special = Flops::getNumberOfMathFunctions(node);
  
  flops.addOps = add; 
  flops.mulOps = multiply; 
  flops.divOps = division; 
  flops.specialOps = special;
}


void ScopeStmtHandler::computeArrayAccesses(SgNode* node)
{
  //read-write props of the loop body 
  ReadWriteProp* rwprop= new ReadWriteProp();

  //we don't really need this class, I can get them from accessanalysis, it is there for validation
  rwprop->run(node);
  rwprop->getReadOnlyVars(readOnlyVars);
  rwprop->getWriteOnlyVars(writeOnlyVars);
  rwprop->getBothReadWrittenVars(bothReadWrittenVars);
  
  //get access and ghost cell info 
  AccessAnalysis* accesses= new AccessAnalysis();
  
  std::vector<SgInitializedName*> loopIndices;
  getLoopIndicesOfParentNest(this, loopIndices);

  //last loop in the list is the cur loop
  accesses->setLoopNo(loopIndices.size()- 1); 
  accesses->setLoopIndices(&loopIndices);

  vector< set< SgVarRefExp*> *> allScalarVars;
  getScalarsOfParentNest(this, allScalarVars);
  accesses->setScalarVars(&allScalarVars);

  ROSE_ASSERT(allScalarVars.size() == loopIndices.size());
  //assumes the size of allScalarVars and loopIndices are the same
  accesses->arrayAccessAnalysis(node, accessInfo, firstAccess);  

  loopIndices.clear();
  allScalarVars.clear();

  delete rwprop; 
  delete accesses;
}

void ScopeStmtHandler::removeDuplicates()
{
  vector<NodeHandler*>* children = this->getChildren();
 
  vector<NodeHandler*>::iterator ch ;

  for(ch= children->begin(); ch!= children->end(); ch++)
    {
      NodeHandler* childNode = *ch;      
      ScopeStmtHandler* child = dynamic_cast<ScopeStmtHandler*>(childNode);

      if(child)
	{
	  removeDuplicatesInScalars(this, child);
	  removeDuplicatesInArrayAccesses(this, child);
	  removeDuplicatesInFlops(this, child);
	}
    } // end of children
}

void ScopeStmtHandler::removeDuplicatesInScalars(ScopeStmtHandler* parent, 
						 ScopeStmtHandler* child)
{
  set<SgVarRefExp*> diffScalars ;
  
  set<SgVarRefExp*>* pScalar = parent->getScalars();
  set<SgVarRefExp*>* cScalar = child ->getScalars();

  set_difference(pScalar->begin(), pScalar->end(), 
  		 cScalar->begin(), cScalar->end(), 
  		 std::inserter(diffScalars, diffScalars.begin()));
  
  pScalar->clear();
  set<SgVarRefExp*>::iterator it; 
  for(it = diffScalars.begin(); it!= diffScalars.end(); it++)
    pScalar-> insert(*it);

}


void ScopeStmtHandler::removeDuplicatesInArrayAccesses(ScopeStmtHandler* parent,
						       ScopeStmtHandler* child)
{ 

  ArraysAccessInfo_t* arraysAccessInfo_ParentLoop = parent ->getAccessInfo();
  ArraysAccessInfo_t* arraysAccessInfo_ChildLoop  = child ->getAccessInfo() ;

  vector<NameComp_t> arrayDeleteList;
  
  ArraysAccessInfo_t::iterator lacc = arraysAccessInfo_ParentLoop->begin();
  for(; lacc != arraysAccessInfo_ParentLoop->end(); lacc++)
    {
      NameComp_t arr_comp= lacc->first;
      AccessListPerArray_t* aListParent = &(lacc->second);
      
      ArraysAccessInfo_t::const_iterator innerlacc = arraysAccessInfo_ChildLoop->find(arr_comp);
      
      vector<string> deleteList;
      
      //same array appears in the inner most loop as well
      //need to look at if the access appears as well
      if(innerlacc != arraysAccessInfo_ChildLoop->end())
	{								
	  const AccessListPerArray_t innerListChild = innerlacc->second;
	  
	  for(AccessListPerArray_t::iterator ac = aListParent->begin(); ac!= aListParent->end(); ac++)
	    {
	      string access_str = ac->first;
	      
	      AccessListPerArray_t::const_iterator innerac = innerListChild.find(access_str);
	      
	      if(innerac != innerListChild.end())
		{
		  //same access appears in the inner loop as well
		  AccessPattern_t inneraccess = innerac ->second;
		  AccessPattern_t* access = &(ac -> second); 
		  
		  access->readCount = access->readCount - inneraccess.readCount;
		  access->writeCount = access->writeCount - inneraccess.writeCount;
		  
		  if(access->readCount <= 0 && access->writeCount <= 0 )
		    {
		      //remove this access because, all the references happen in the inner loop
		      deleteList.push_back(access_str);
		    }
		}
	    }
	  
	  for(vector<string>::iterator it = deleteList.begin(); it!= deleteList.end(); it++)
	    {
	      string access_str = *it;
	      AccessListPerArray_t::iterator ac = aListParent->find(access_str);
	      
	      if(ac != aListParent->end())
		aListParent->erase(ac);
	    }
	  
	  if(aListParent->size() == 0)
	    {
	      //remove the array-comp pair as well
	      aListParent->clear();
	      arrayDeleteList.push_back(arr_comp);
	    }
	} // end of if access appears in the inner loop
      
    } //end of var-ghost pairs
  
  for(vector<NameComp_t> ::iterator it = arrayDeleteList.begin(); it!= arrayDeleteList.end(); it++)
    {
      NameComp_t array_comp = *it;
      ArraysAccessInfo_t::iterator item = arraysAccessInfo_ParentLoop->find(array_comp);
      if(item != arraysAccessInfo_ParentLoop->end())
	arraysAccessInfo_ParentLoop->erase(item);
    }
}


void ScopeStmtHandler::removeDuplicatesInFlops(ScopeStmtHandler* parent, 
					       ScopeStmtHandler* child)
{  
  Flops_t* flops_parent = parent -> getFlops();
  Flops_t* flops_child = child -> getFlops(); 
    
  if(flops_parent->addOps > 0){
    flops_parent -> addOps -= flops_child-> addOps; 
    ROSE_ASSERT(flops_parent-> addOps >=0 );
  }
  if(flops_parent->mulOps > 0){
    flops_parent -> mulOps -= flops_child-> mulOps; 
    ROSE_ASSERT(flops_parent-> mulOps >=0 );
  }
  if(flops_parent->divOps > 0){
    flops_parent -> divOps -= flops_child-> divOps; 
    ROSE_ASSERT(flops_parent-> divOps >=0 );
  }
  if(flops_parent->specialOps > 0){
    flops_parent -> specialOps -= flops_child-> specialOps; 
    ROSE_ASSERT(flops_parent-> specialOps >=0 );
  }
}

void ScopeStmtHandler::getLoopIndicesOfParentNest(ScopeStmtHandler* cur_handler,
						  std::vector<SgInitializedName*>& loopIndices)
{
  //This function gets all the loop indices from parent's loop
  NodeHandler* parentNode = cur_handler -> getParent() ; 
  ScopeStmtHandler* parent = dynamic_cast<ScopeStmtHandler*>(parentNode);

  if(parent)
    getLoopIndicesOfParentNest(parent, loopIndices);

  SgInitializedName* index =  cur_handler -> getLoopIndex();
  loopIndices.push_back(index);
}

void ScopeStmtHandler::getScalarsOfParentNest(ScopeStmtHandler* cur_handler,
					      vector< std::set<SgVarRefExp*> * >& allScalarVars)
{
  NodeHandler* parentNode =  cur_handler -> getParent() ; 
  ScopeStmtHandler* parent = dynamic_cast<ScopeStmtHandler*>(parentNode);
  
  if(parent)
    getScalarsOfParentNest(parent, allScalarVars);

  std::set<SgVarRefExp*>* scalars = cur_handler -> getScalars();
  allScalarVars.push_back(scalars);
}

void ScopeStmtHandler::dumpScalarInfo(set<SgInitializedName*>& V,
				      ScalarRWMap_t& scalarRW,
				      FirstAccessMap_t& firstAccess,
				      const string xml_tag)
{

  for(set<SgInitializedName*>::iterator it= V.begin(); it != V.end() ; it++)
    {
      SgInitializedName* name = (*it);
      ROSE_ASSERT(name);

      SgDeclarationStatement* decl_stmt = name-> get_declaration(); 
      ROSE_ASSERT(decl_stmt);

      //For the structs, we would like to report variable names togather with the struct name
      string struct_name = "";

      SgClassDefinition* class_def = getEnclosingClassDefinition(decl_stmt);

      if(class_def != NULL){
	SgClassDeclaration* class_decl = class_def -> get_declaration();

	if(class_decl != NULL && isSgDerivedTypeStatement(class_decl))
	  struct_name = class_def-> get_qualified_name(); 
      }
      //end of getting struct name 
      
      ScalarRWMap_t::const_iterator rw = scalarRW.find(name);
      int reads = 0;
      int writes = 0;

      if(rw != scalarRW.end()) //found                                                                  
	{
	  AccessPattern_t rwacc = rw->second;
	  reads = rwacc.readCount; 
	  writes = rwacc.writeCount++;
	}

      string accesstype = "readwrite";
      if(reads > 0 && writes ==0)
	accesstype = "readonly";
      else if(reads == 0 && writes >0)
	accesstype = "writeonly";
      else 
	{
	  //scalars do not have a component, so create an empty component
	  FirstAccessMap_t::iterator fa = firstAccess.find(make_pair(name, ""));
	  if(fa != firstAccess.end())
	    {
	      FirstAccess_t acc = fa->second;
	      if(!acc.firstRead)
		accesstype = "writeread";
	    }
	}

      ostringstream read_str, write_str ;
      read_str << reads; 
      write_str << writes;

      string const_name = "false";
      
      if((decl_stmt -> get_declarationModifier().get_typeModifier().get_constVolatileModifier().isConst()))
	const_name = "true";
      
      SgType* type = name->get_type();
      if(isConstType(type)) const_name = "true";
      ROSE_ASSERT(type);
      
      string s_name; 
      s_name = name ->get_name().str();      

      string t_name = type->unparseToString();
      
      std::vector < pair < string, string> >attr;
      
      attr.push_back(make_pair("name", s_name));
      if(struct_name != "")
	attr.push_back(make_pair("struct", struct_name));
      attr.push_back(make_pair("datatype", t_name));
      attr.push_back(make_pair("isConstant", const_name));
      attr.push_back(make_pair("accesstype", accesstype));
      attr.push_back(make_pair("reads", read_str.str()));
      attr.push_back(make_pair("writes", write_str.str()));
      OutputFormat::outputBeginEnd(xml_tag, attr);
    }
}

void ScopeStmtHandler::dumpVectorInfo(set<NameComp_t> vars,
				      string rose_accesstype)
{  
  for(set<NameComp_t>::const_iterator rIt = vars.begin(); rIt != vars.end() ; rIt++)
    {
      const NameComp_t arr_comp = *rIt;

      std::vector< pair<string, string> > attributes;       
      const SgInitializedName* array_name = arr_comp.first; 
      const string component = arr_comp.second; 
      //const SgType* type = array_name->get_type()->findBaseType();
      const SgType* type;
      string type_str="";
      if(array_name->get_type()){
        if(isSgPointerType(array_name->get_type())){
          type = isSgPointerType(array_name->get_type())->get_base_type();
          type_str+= type->unparseToString();
        }
      }

      ArraysAccessInfo_t::const_iterator arrAccInfo_itr = accessInfo.find(arr_comp);

      if(arrAccInfo_itr != accessInfo.end()) //not found 
	{ 

	  AccessListPerArray_t accessList = arrAccInfo_itr->second; 

	  //Didem: Dec 19, 2013
	  //correct readonly, writeonly analysis based on the access analysis
	  //ROSE's collectReadWriteRefs is not reliable for structs
	  //This step also fixes readwrite cases where it should be really be writeread
	  string accesstype = rose_accesstype; 

	  FirstAccessMap_t::iterator fa_itr = firstAccess.find(arr_comp);
	  if(fa_itr != firstAccess.end()){
	    FirstAccess_t facc = fa_itr -> second; 
	    accesstype = fixAccessType(accessList, facc, rose_accesstype);
	  }
	  //end of fix access type

	  //For the structs, we would like to report variable names together with the struct name
	  SgDeclarationStatement* decl_stmt = array_name-> get_declaration(); 
	  ROSE_ASSERT(decl_stmt);

	  string struct_name = "";	  

	  SgClassDefinition* class_def = getEnclosingClassDefinition(decl_stmt);
	  if(class_def != NULL){
	  
	    SgClassDeclaration* class_decl = class_def -> get_declaration();

	    if(class_decl != NULL && isSgDerivedTypeStatement(class_decl))
	      struct_name = class_def-> get_qualified_name(); 
	  }
	  //end of getting struct name 
	  
	  //report array attributes
	  attributes.push_back(make_pair("name", array_name->get_name().str()));
	  if(struct_name != "")
	    attributes.push_back(make_pair("struct", struct_name));
	  attributes.push_back(make_pair("component", component));
	  attributes.push_back(make_pair("datatype", type_str));
	  attributes.push_back(make_pair("accesstype", accesstype));
	  
	  //output the array attributes if there is at least one access
	  OutputFormat::outputBegin("array", attributes);

	  //starting accesses 
	  for(AccessListPerArray_t::iterator ac = accessList.begin(); ac!= accessList.end(); ac++)
	    {
	      const string access_str = ac->first; 	        
	      const AccessPattern_t access = ac->second;
	      const string offset = Utils::convertStringVecToString(access.offset);
	      const string dependentLoopVar = Utils::convertStringVecToString(access.dependentLoopVar);
	      //const string relative = (access.isRelative)? "true" : "false"; 
	      
	      ostringstream readstr, writestr ;
	      int reads = access.readCount; 
	      int writes = access.writeCount; 

	      readstr << reads;
	      writestr << writes;
	      
	      std::vector< pair<string, string> > access_attr;
	      //access_attr.push_back(make_pair("index", access_str));
	      //access_attr.push_back(make_pair("isRelative", relative ));
	      access_attr.push_back(make_pair("offset", offset));
	      access_attr.push_back(make_pair("dependentloopvar", dependentLoopVar ));
	      access_attr.push_back(make_pair("reads", readstr.str()));
	      access_attr.push_back(make_pair("writes", writestr.str()));
	      
	      OutputFormat::outputBeginEnd("access", access_attr);
	    }
	}
      else
	{
	  continue;
	  //attributes.push_back(make_pair("ghost", ghostVal));
	  //outputBegin("array", attributes);
	}
      OutputFormat::outputEnd("array");

    }//end of for 
}

string ScopeStmtHandler::fixAccessType(AccessListPerArray_t accessList,
				       FirstAccess_t fa,
				       string rose_accesstype)
{
  //correct readonly, writeonly analysis based on the access analysis
  //ROSE's collectReadWriteRefs is not reliable for structs
	
  string accesstype = rose_accesstype; 
  int readcounts =0;
  int writecounts = 0;

  for(AccessListPerArray_t::iterator ac = accessList.begin(); ac!= accessList.end(); ac++)
    {
      const AccessPattern_t access = ac -> second;
      readcounts = access.readCount; 
      writecounts = access.writeCount; 
    }

  if(writecounts != 0 && readcounts != 0 )
    {
      if(!fa.firstRead)
	accesstype = "writeread";
      else 
	accesstype = "readwrite";
    }
  if(writecounts == 0 && readcounts != 0 ){
    accesstype = "readonly";
    if(!fa.firstRead){
      cerr << "WARNING: should be first read because it is readonly" << endl; 
    }
    //if(rose_accesstype != accesstype)
    //  cerr << "WARNING: rose's accesstype doesn't match with mine" << endl; 
  }

  if(writecounts != 0 && readcounts == 0 ){
    accesstype = "writeonly";

    if(fa.firstRead){
      cerr << "WARNING: should be first written because it is writeonly" << endl; 
    }
    //if(rose_accesstype != accesstype)
    //  cerr << "WARNING: rose's accesstype doesn't match with mine" << endl; 
  }

  return accesstype; 
}

//eof








