
#include <rose.h>
#include "./OutputFormatting.h"

#include "VarSym.hh"
#include "./ExaSatOptions.h"


using namespace std;


void OutputFormat::dump (const ASTtools::VarSymSet_t& V, const std::string& tag, const string xml_tag)
{
  if(!ExaSatOptions::xmlOutput)
    cout << tag << " " << ASTtools::toString (V)  << endl ;
  else
    {
      for(ASTtools::VarSymSet_t::const_iterator it= V.begin(); it != V.end() ; it++)
        {
          const SgVariableSymbol* sym = (*it);
          if(sym){
            const SgInitializedName* name = sym->get_declaration();
            ROSE_ASSERT(name);
            string s_name = name ->get_name().str();
            
	    std::vector < pair < string, string> >attr;
	    
	    attr.push_back(make_pair("name", s_name));
	    
	    outputBeginEnd(xml_tag, attr);

          }
        }
    }
}
void OutputFormat::dump (const set<SgInitializedName*>& V, 
			 const std::string& tag, const string xml_tag)
{
  if(!ExaSatOptions::xmlOutput){
    //cout << tag << " " << ASTtools::toString (V)  << endl ;
  }
  else
    {
      for(set<SgInitializedName*>::const_iterator it= V.begin(); it != V.end() ; it++)
        {
	  const SgInitializedName* name = (*it);
	  ROSE_ASSERT(name);
	  SgDeclarationStatement* var_decl = name -> get_declaration();
	  string const_name = "false";

	  if((var_decl -> get_declarationModifier().get_typeModifier().get_constVolatileModifier().isConst()))          
	    const_name = "true";

	  SgType* type = name->get_type();
	  ROSE_ASSERT(type);

	  string s_name = name ->get_name().str();

	  string t_name = type->unparseToString(); 

	  std::vector < pair < string, string> >attr;

	  attr.push_back(make_pair("name", s_name));
	  attr.push_back(make_pair("elementtype", t_name));
	  attr.push_back(make_pair("isConstant", const_name));

	  outputBeginEnd(xml_tag, attr);

        }
    }
}

void OutputFormat::outputBegin(string tag, 
			       std::vector< pair<string, string> > attributes)
{
  if(!ExaSatOptions::xmlOutput)
    return ; 
  
  cout << "<" << tag ;
  
  std::vector< pair <string, string> > ::const_iterator at = attributes.begin();

  for( ; at != attributes.end() ; at++)
    {
      pair < string, string> attribute = (*at);
      
      cout << " "<<attribute.first << "=\"" << attribute.second << "\""; 
    }
  cout << ">" << endl; 
}

void OutputFormat::outputBeginEnd(string tag, 
			       std::vector< pair<string, string> > attributes)
{
  if(!ExaSatOptions::xmlOutput)
    return ; 
  
  cout << "<" << tag ;
  
  std::vector< pair <string, string> > ::const_iterator at = attributes.begin();

  for( ; at != attributes.end() ; at++)
    {
      pair < string, string> attribute = (*at);
      
      cout << " "<<attribute.first << "=\"" << attribute.second << "\""; 
    }
  cout << "/>" << endl; 
}

void OutputFormat::outputBegin(string tag, string name, 
			       string component, string value)
{
  if(!ExaSatOptions::xmlOutput)
    return ; 
  
  cout << "<" << tag ;
  
  if(!name.empty())
    cout << " name=\"" << name << "\""; 
  
  if(!component.empty())
    cout << " component=\"" << component << "\""; 
  
  if(!value.empty())
    cout << " value=\"" << value << "\""; 
    
  cout << ">" << endl; 
}


void OutputFormat::outputEnd(string tag)
{
    if(!ExaSatOptions::xmlOutput)
      return ; 
    
    //ex: </add>
    if(!tag.empty())
      cout << "</" << tag << ">" << endl; 
    else
      cout << "/>" << endl; 
}
