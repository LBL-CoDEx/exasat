#include <rose.h>
#include "IfHandler.h"

#include "../scalars/Scalars.h"
#include "../utils/OutputFormatting.h"


using namespace std;
using namespace SageInterface;
using namespace OutputFormat; 

IfHandler::IfHandler(SgNode* node, bool iforelse)
{
  myNode = node;
  visited = false;
  parent = NULL;
  isTrueBody = iforelse;

  ifAttr.conditional = NULL;
}

IfHandler::~IfHandler()
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


void IfHandler::process()
{

  SgNode* node = this-> myNode; 

  SgIfStmt* ifstmt = isSgIfStmt(node);
  SgBasicBlock* body; 

  //Mar 2015 Tan: I modified the code that gets the true and false bodies of an if statement so it works with C code as well
  //in particular, an if statement in C does not necessarily go with a basic block. 
  //in that case, we create a basic block to have a unified IR for both C and Fortran 
  if(isTrueBody){ 
    SgStatement* true_stmt = ifstmt -> get_true_body();
    SgBasicBlock* true_body=NULL;
    if(!isSgBasicBlock(true_stmt)){
      true_body = buildBasicBlock(true_body);
      true_body->set_parent(ifstmt);
      ifstmt->set_true_body(true_body);
    }else true_body= isSgBasicBlock(true_stmt);
    ROSE_ASSERT(true_body);

    int lineno = node ->get_file_info()->get_line();
    ifAttr.lineno = lineno;

    SgStatement* condition = ifstmt-> get_conditional();
    ifAttr.conditional = condition;

    body = isSgBasicBlock(true_body);
  }
  else 
    {
      SgStatement* false_stmt = ifstmt -> get_false_body();
      SgBasicBlock* false_body=NULL;
      if(!isSgBasicBlock(false_stmt)){
        false_body= buildBasicBlock(false_stmt);
        false_body->set_parent(ifstmt);
        ifstmt->set_false_body(false_body);
      }else false_body= isSgBasicBlock(false_stmt);
      ROSE_ASSERT(false_body);
//May 04: Tan fixed this bug on behalf of Cy. Exasat now reports the line number of the if statement for both branches
      int lineno = /*false_body*/node->get_file_info()->get_line();
      ifAttr.lineno = lineno;

      body = false_body;
    }

  computeFlops(isSgNode(body));

  //scalar variables in the loop
  //because of the else statments, we do not include the scalars of the condition 
  //in the current if but include it in the parent
  Scalars::getActiveScalarVars(isSgNode(body), scalarVars);

  computeArrayAccesses(isSgNode(body));

}

void IfHandler::closeXMLTag()
{
  if(isTrueBody)
    outputEnd("if");
  else 
    outputEnd("else");
}


inline string find_replace(string &originalStr, string replacedSubStr, string replacementStr){
    size_t position= originalStr.find(replacedSubStr);
    while(position!= std::string::npos){
      originalStr.replace(position, replacedSubStr.length(), replacementStr);
      position += replacementStr.length();
      position= originalStr.find(replacedSubStr, position);
    }
    return originalStr;
}

void IfHandler::dump()
{
  //if attributes 
  std::vector< pair<string, string> > attributes;
  
  ostringstream linenum_str ;
  linenum_str << ifAttr.lineno;

  if(isTrueBody)
    attributes.push_back(make_pair("linenum", linenum_str.str()));
  else 
    attributes.push_back(make_pair("iflinenum", linenum_str.str()));

  if(ifAttr.conditional){
    string conditionalStr(ifAttr.conditional->unparseToString());
    find_replace(conditionalStr, "&", "&amp;");
    find_replace(conditionalStr, "<", "&lt;");
    find_replace(conditionalStr, ">", "&gt;");
    find_replace(conditionalStr, "\"", "&quot;");
    attributes.push_back(make_pair("conditional", conditionalStr));
    //attributes.push_back(make_pair("conditional", ifAttr.conditional->unparseToString()));
  }
  
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

  //Output if attributes 
  if(isTrueBody)
    outputBegin("if", attributes);
  else 
    outputBegin("else", attributes);

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








