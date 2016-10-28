
#ifndef OUTPUTFORMAT_H
#define OUTPUTFORMAT_H

#include <rose.h>

#include <vector>
#include <string>
#include "ASTtools.hh"

using namespace std;
using namespace SageInterface; 

namespace OutputFormat
{
  void dump (const ASTtools::VarSymSet_t& V, const std::string& tag, const string xml_tag);
  void dump (const set<SgInitializedName*>& V, const std::string& tag, const string xml_tag);

  void outputBegin(string tag, string name="", 
		   string component="", string value="");

  void outputBegin(string tag, std::vector< pair <string, string> > attributes);
  void outputBeginEnd(string tag, std::vector< pair <string, string> > attributes);

  
  void outputEnd(string tag="");
  
};

#endif 
