#!/usr/bin/env python

""" Analyze the XML generated by the Compiler Analysis pass.

    This module contains classes that do various types of analysis over
    the code, including flop count, state and streaming variable access,
    working set, and memory traffic estimates.
"""
__author__ = "Cy Chan"
__copyright__ = "Copyright 2016, The Regents of the University of California, through Lawrence Berkeley National Laboratory"
__credits__ = ["Cy Chan"]
__license__ = "Modified BSD License (see LICENSE file)"
__version__ = "2.0"
__maintainer__ = "Cy Chan"
__email__ = "cychan@lbl.gov"
__status__ = "Production"

import sys
import os
from distutils.util import strtobool

from analyze import StaticAnalysis, TableCondsChecker
from parser import KeyValXMLParser, to_sym_dict

# default args
def default_args():
  return {
    "xml"          : None,
    "fns"          : None,
    "polly"        : None,
    "symsubs"      : None,
    "namesubs"     : None,
    "params"       : None,
    "block_params" : None,
    "conds"        : None,
    "machine"      : "../../examples/machine.xml",
    "subparams"    : "False",
  }

# setup for CNS code
def cns_args():
  return {
    "xml"          : "../../examples/cns-smc/xml/advance-flat.xml",
    "fns"          : None,
#    "polly"        : "../../examples/cns-smc/xml/advance-flat.polly.xml",
    "polly"        : "",
    "symsubs"      : "../../examples/cns-smc/symsubs.xml",
    "namesubs"     : "../../examples/cns-smc/namesubs.xml",
    "params"       : "../../examples/cns-smc/params.xml",
    "block_params" : "../../examples/cns-smc/block_params.xml",
    "conds"        : "../../examples/cns-smc/conds.xml",
    "machine"      : "../../examples/machine.xml",
    "subparams"    : os.getenv("subparams", "False"),
  }

# setup for SMC code
def smc_args():
  # use most of the CNS setup
  result = cns_args()
  # change the source input files
  result["xml"] = "../../examples/cns-smc/xml/smc/advance-smc-modified.xml"
#  result["polly"] = "../../examples/cns-smc/xml/smc/advance-smc-modified.polly.xml"
  result["polly"] = ""
  return result

# setup for HPGMG code
def hpgmg_args():
  return {
    "xml"          : "../../examples/hpgmg/xml/hpgmg.27pt.s2.xml",
    "fns"          : "smooth",
    "polly"        : "",
    "symsubs"      : "../../examples/hpgmg/symsubs.xml",
    "namesubs"     : "../../examples/hpgmg/namesubs.xml",
    "params"       : "../../examples/hpgmg/params.xml",
    "block_params" : "../../examples/hpgmg/block_params.xml",
    "conds"        : "../../examples/hpgmg/conds.xml",
    "machine"      : "../../examples/machine.xml",
    "subparams"    : os.getenv("subparams", "False"),
  }

def get_env_args():
  result = default_args()
  for tag in ["xml", "polly", # which XML files to analyze
              "fns", # which functions to analyze
              "symsubs", "namesubs", # substitutions made upon parsing the XML
              # problem and machine parameters used to evaluate performance
              "params", "block_params", "conds", "machine",
              # bool: T: substitute parameters first for numeric results (faster)
              #       F: generate symbolic results in terms of the parameters (slower)
              "subparams", 
             ]:
    val = os.getenv(tag, None)
    if val:
      result[tag] = val
  return result

def load_args(args):
  return (
    { "xml"             : args["xml"],
      "polly_xml"       : args["polly"],
      "symsubs"         : to_sym_dict(KeyValXMLParser(args["symsubs"]).items),
      "namesubs"        : to_sym_dict(KeyValXMLParser(args["namesubs"]).items),
    },
    { "fns"             : set(args["fns"].split(',')),
      "params"          : to_sym_dict(KeyValXMLParser(args["params"]).items),
      "block_params"    : to_sym_dict(KeyValXMLParser(args["block_params"]).items),
      "machine"         : dict(KeyValXMLParser(args["machine"], float).items),
      "conds_chk"       : TableCondsChecker(KeyValXMLParser(args["conds"], float).items),
      "flag_sub_params" : strtobool(args["subparams"]),
    })

def main(cl_args):

  if len(cl_args) > 1:
    if cl_args[1] == "cns":
      args = cns_args()
    elif cl_args[1] == "smc":
      args = smc_args()
    elif cl_args[1] == "hpgmg":
      args = hpgmg_args()
  else:
    args = get_env_args()

  (sa_kw_args, dump_kw_args) = load_args(args)

  # parse the XML files and do substitutions
  sa = StaticAnalysis(**sa_kw_args)

  # do performance analysis at granularity of first-level loops in all functions
  sa.dump(**dump_kw_args)
 
if __name__ == '__main__':
  main(sys.argv)
