#include "rose.h"
#include "DefUseAnalysis.h"
#include <string>
#include <iostream>
using namespace std;

int main( int argc, char * argv[] ) 
{
  vector<string> argvList(argv, argv + argc);
  SgProject* project = frontend(argvList);

  // Call the Def-Use Analysis
  DFAnalysis* defuse = new DefUseAnalysis(project);
  bool debug = false;
  defuse->run(debug);
  // Output def use analysis results into a dot file
  // defuse->dfaToDOT();

  // Find all variable references
  NodeQuerySynthesizedAttributeType doLoops = NodeQuery::querySubTree(project, V_SgFortranDo); 
  NodeQuerySynthesizedAttributeType::const_iterator loop= doLoops.begin();

  for(; loop != doLoops.end() ; loop++)
    {
      NodeQuerySynthesizedAttributeType vars = NodeQuery::querySubTree(*loop, V_SgVarRefExp); 
      NodeQuerySynthesizedAttributeType::const_iterator i = vars.begin();
      
      for (; i!=vars.end();++i) 
	{
	  SgVarRefExp * varRef = isSgVarRefExp(*i);
	  SgInitializedName* initName = isSgInitializedName(varRef->get_symbol()->get_declaration());
	  std::string name = initName->get_qualified_name().str();
	  SgType* type = initName->get_type();
	
	  //	  if(isSgPointerType(type) || name == "ux" ){
	    if(name == "cons" || name == "flux" || name == "q" || name == "unp1"){
	    // Find reaching definition of initName at the control flow node varRef
	    vector<SgNode* > vec = defuse->getDefFor(varRef, initName);
	    //ROSE_ASSERT (vec.size() >0 ); // each variable reference must have a definition somewhere
	    
	    // Output each reaching definition node and the corresponding statement.
	    std::cout<<"---------------------------------------------"<<std::endl;
	    std::cout << vec.size() << " definition entry/entries for " << varRef->unparseToString() <<  
	      " @ line " << varRef->get_file_info()->get_line()<<":"<<varRef->get_file_info()->get_col() 
		      << "????these are"<< std::endl;
	    for (size_t j =0; j<vec.size(); j++)
	      {
		cout<<vec[j]->class_name()<<" "<<vec[j]<<endl;
		SgStatement * def_stmt = SageInterface::getEnclosingStatement(vec[j]);
		ROSE_ASSERT(def_stmt);
		cout<<def_stmt->unparseToString()<<"  @ line "<<def_stmt->get_file_info()->get_line()
		    <<":"<<def_stmt->get_file_info()->get_col() <<endl;
	      }
	  }
	}
    }
      return 0;
}
