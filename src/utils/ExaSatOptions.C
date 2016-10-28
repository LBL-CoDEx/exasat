#include "ExaSatOptions.h"

#include <CommandOptions.h>

#include <cassert>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>

using namespace std;

//init static members 
ExaSatOptions* ExaSatOptions::inst=0;
bool ExaSatOptions::xmlOutput= true; 
bool ExaSatOptions::generateAST= false;

ExaSatOptions* ExaSatOptions::GetInstance()
{
  if(inst == 0)
    {
      inst = new ExaSatOptions();
    }
  return inst;
}

bool ExaSatOptions::HasOption(const std::string& opt)
{
  for(std::vector<std::string>::const_iterator i = opts.begin(); i!= opts.end(); i++)
    {
      if((*i).substr(0, opt.size()) == opt)
	return true;
    }
  return false; 
}

void ExaSatOptions::SetExaSatOptions  (const std::vector<std::string>& opts)
{
  CmdOptions::GetInstance()->SetOptions(opts);
  assert (!opts.empty());
  this->opts = opts;
  this->opts.erase(this->opts.begin());

  isXMLOutput();
  isGenerateAST();

}

bool ExaSatOptions::isXMLOutput()
{
  static int r = 0;
  if(r == 0)
    {
      if (HasOption("-xml")){
	xmlOutput = true; 
	r=1;
      }
      else{ 
	r=-1; 
      }
    }
  return r==1;
}

bool ExaSatOptions::isGenerateAST()
{
  static int r = 0;
  if(r == 0)
    {
      if (HasOption("-ast")){
	generateAST = true; 
	r=1;
      }
      else{ 
	r=-1; 
      }
    }
  return r==1;
}

/*
what ExaSatOptpions::getBoxSizes()
{
  getBoxSize("-bx");
  getBoxSize("-by");
  getBoxSize("-bz");
}

int ExaSatOptions::getBoxSize(string bxyz)
{
  vector<string>::const_iterator config = CmdOptions::GetInstance()->GetOptionPosition(bxyz);

  size_t size = 0;

  if (config == CmdOptions::GetInstance()->opts.end())
    return size;

  assert (config != CmdOptions::GetInstance()->opts.end());

  string sizeStr = (*config).c_str() ;

  int cutAt = sizeStr.find_first_of("=");

  if(cutAt > 0){ //box size
    //get the substring                                                                                                                                               
    string tmp = sizeStr.substr(cutAt+1, sizeStr.length());
    sscanf(tmp.c_str(), "%zu", &size);
  }

  return size;
}
*/
