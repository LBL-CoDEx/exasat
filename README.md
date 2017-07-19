ExaSAT: Exascale Static Analysis Tool
======
---
### Prerequisites: ###
- ROSE and its prerequisites (check website):
  - gcc-4.4.7, Boost, Java
  - python
- OMEGA and CHILL (optional)

### Compiler Analysis Component Usage: ###
- To build, run make from src/
-  Set environment variables:
    - `FILTER=<fname>`: fname is function that we want to analyze, or
    - `NOFILTER=1`: analyze all functions in input file
  - Add Boost and Java lib dirs to `LD_LIBRARY_PATH`
  - Run: `genCodeDescript <input_file>`
    - Produces XML on stdout
  - Example: run `./genxml.py` from `examples/cns-smc/`
    - Crawls the `inputs/` sub-directory and tries to run the compiler analysis on all files
      - May produce errors for some because it can’t find some dependent module files (it’s safe to ignore these)
    - Stores the produced XMLs in the `xml_new/` sub-directory

### Performance Model Component (old version) Usage: ###

##### Set environment variables: #####
  - `problem_xml=<problem-XML-file>`:
    The problem XML specifies problem parameters and software
    optimizations to be evaluated, such as problem size, number of
    ghost cells, how many chemical species to evaluate, cache blocking
    parameters, whether to utilize streaming writes and non-temporal
    access hints.  Defaults will be chosen if no file is specified.
  - `machine_xml=<machine-XML-file>`:
    The machine XML specifies machine parameters, such as memory
    bandwidth, FP compute throughput, number registers, cache sizes, word
    size, cache line size, and relative costs of arithmetic operations
    (+,-,*,/,specials) and memory operations (R,W,RW).  Defaults will
    be chosen if no file is specified.
  - `old=0` or `old=1`:
    What format XML to expect

##### Generate analysis spreadsheet (in tools/post-old/ directory): #####
- `./analyze.py <xml-input> <tsv-output>`
- Example: problem_xml=../../examples/problem.xml ./analyze.py ../../examples/cns-smc/xml/advance-nomod.xml advance-nomod.tsv
  - Generates an analysis spreadsheet advance-nomod.tsv
- Example script: `./gentsvs.sh`
  - Analyze advance routines from CNS and SMC codes as well as chemistry routines from SMC code
- Can generate dependency graphs by setting `options.flag_graph=True` in common.py

##### Or run in script mode: #####
  - Call analysis from another Python script through analyze and runner modules
    - Step 1: analyze an input XML file to recognize loop structures and calculate working sets, etc. (analyze module)
    - Step 2: use static analysis info to run a performance model for a given problem and hardware configuration (runner module)
  - Multiple example analyses given by running `./runModel.py <mode>`
    - `mode = {micro, cns, smc, smc-regs, smc-summary, diffterm}`

---
## Copyright ##
"Exascale Static Analysis Tool (ExaSAT)" Copyright (c) 2016, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy).  All rights reserved.

If you have questions about your rights to use or distribute this software, please contact Berkeley Lab's Innovation & Partnerships Office at  IPO@lbl.gov.

NOTICE.  This Software was developed under funding from the U.S. Department of Energy and the U.S. Government consequently retains certain rights. As such, the U.S. Government has been granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, distribute copies to the public, prepare derivative works, and perform publicly and display publicly, and to permit other to do so.
