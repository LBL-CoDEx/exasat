#include "rose.h"
#include "DefUseAnalysis.h"
#include "LivenessAnalysis.h"

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
  LivenessAnalysis* liv = new LivenessAnalysis(debug, (DefUseAnalysis*)defuse);

  NodeQuerySynthesizedAttributeType funs = NodeQuery::querySubTree(project, V_SgFunctionDefinition);
  NodeQuerySynthesizedAttributeType::const_iterator func_itr= funs.begin();
  
  for(; func_itr != funs.end() ; func_itr++){
    SgFunctionDefinition* func = isSgFunctionDefinition(*func_itr);
    ROSE_ASSERT(func);
    
    bool abortme = false; 
    FilteredCFGNode< IsDFAFilter> rem_source = liv->run(func, abortme); 
    
    //if(!abortme) 
    //   cerr << "aborting " << endl; 
    // Find all variable references
    NodeQuerySynthesizedAttributeType doLoops = NodeQuery::querySubTree(func, V_SgFortranDo); 
    NodeQuerySynthesizedAttributeType::const_iterator loop= doLoops.begin();
    
    for(; loop!= doLoops.end(); loop++)
      {
	std::vector<SgInitializedName*> liveIns = liv->getIn(*loop);
	std::vector<SgInitializedName*> liveOuts = liv->getOut(*loop);
        
	cout << "live outs " << endl; 
	for (std::vector<SgInitializedName*>::const_iterator iter = liveOuts.begin();
	     iter!=liveOuts.end(); iter++)
	  cout<<" "<<(*iter)->get_name().getString()<<" ";
	cout<<endl;
	
	cout << "live ins " << endl ; 
	for (std::vector<SgInitializedName*>::const_iterator iter = liveIns.begin();
	     iter!=liveIns.end(); iter++)
	  cout<<" "<<(*iter)->get_name().getString()<<" ";
	cout<<endl;
      }
  }
  delete liv;
  delete defuse;
      
  return 0;
}
