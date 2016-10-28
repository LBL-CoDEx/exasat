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
  //DFAnalysis* defuse = new DefUseAnalysis(project);
  //bool debug = false;
  //defuse->run(debug);
  // Output def use analysis results into a dot file
  // defuse->dfaToDOT();
  //LivenessAnalysis* liv = new LivenessAnalysis(debug, (DefUseAnalysis*)defuse);
      
  return 0;
}
