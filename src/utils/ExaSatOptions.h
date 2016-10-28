#ifndef OPTIONS_H
#define OPTIONS_H

#include <iostream>
#include <string>
#include <vector>

#include <AstInterface.h>


class ExaSatOptions
{ 
 private:
  static ExaSatOptions *inst;
  
 public:
  
  ExaSatOptions() : opts() {}

  static bool xmlOutput; 
  static bool generateAST;

  std::vector<std::string> opts; // So its .end() method is accessible
  
  static ExaSatOptions* GetInstance () ;
  
  void SetExaSatOptions( const std::vector<std::string>& argvList);

  bool isXMLOutput();
  bool isGenerateAST();

  bool HasOption( const std::string& opt);
  void PrintUsage(std::ostream& stream) const ;

};


#endif 
