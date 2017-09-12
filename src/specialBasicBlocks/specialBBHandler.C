#include <rose.h>
#include "specialBBHandler.h"

#include "../scalars/Scalars.h"
#include "../utils/OutputFormatting.h"


using namespace std;
using namespace SageInterface;
using namespace OutputFormat; 

SpecialBBHandler::SpecialBBHandler(SgNode* node)
{
  myNode = node;
  visited = false;
  parent = NULL;
}

SpecialBBHandler::~SpecialBBHandler()
{
  //clean up 
  accessInfo.clear();
  readOnlyVars.clear();
  writeOnlyVars.clear();
  bothReadWrittenVars.clear();
  scalarVars.clear();
  firstAccess.clear();
  children.clear();
}

void SpecialBBHandler::process()
{
  SgNode* node = this-> myNode; 

  ROSE_ASSERT(isSgBasicBlock(node));
  
  SgBasicBlock* body= isSgBasicBlock(node);

  computeFlops(isSgNode(body));

  computeArrayAccesses(isSgNode(body));

  Scalars::getActiveScalarVars(isSgNode(body), scalarVars);
}

void SpecialBBHandler::closeXMLTag()
{
  //do nothing
}

void SpecialBBHandler::dump()
{
  SgName func_name = isSgFunctionDefinition(isSgBasicBlock(this-> myNode)->get_parent()) -> get_declaration()->get_name();
  string funcStr= "function name=\"";
  funcStr = funcStr + func_name.str();
  funcStr = funcStr + "\"";

  //if attributes 
  std::vector< pair<string, string> > attributes;
  
  ostringstream count_add_str, count_mul_str, count_special_str, count_div_str; 
  //flops 
  count_add_str << flops.addOps;
  count_mul_str << flops.mulOps;
  count_div_str << flops.divOps;
  count_special_str << flops.specialOps;
  
  attributes.push_back(make_pair("adds", count_add_str.str()));
  attributes.push_back(make_pair("multiplies", count_mul_str.str() ));
  attributes.push_back(make_pair("divides", count_div_str.str() ));
  attributes.push_back(make_pair("specials", count_special_str.str()));
  outputBegin(funcStr, attributes);

  //we need to remove the dublicates first for scalars because those are varref
  //when we convert them into iname, they are unique
  ScalarRWMap_t scalarRW;
  FirstAccessMap_t scalarFirstAccess;
  set<SgInitializedName*> scalarInames; 
  
  Scalars::getScalarReadWriteCounts(scalarRW, scalarFirstAccess, scalarVars);
  Scalars::getScalarInames(scalarInames, scalarVars);

  //dump scalar variables
  dumpScalarInfo(scalarInames, scalarRW, scalarFirstAccess, "scalar");
  
  dumpVectorInfo(readOnlyVars, "readonly");
  dumpVectorInfo(writeOnlyVars, "writeonly");
  dumpVectorInfo(bothReadWrittenVars, "readwrite");
}


//eof








