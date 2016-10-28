#!/bin/bash

# need to point these to the right locations
export GCD_BINARY=$HOME/projects/exasat/src/genCodeDescript
SCRIPT=$HOME/projects/exasat/scripts/crawl.py 
CASTRO=$HOME/projects/Castro
BOXLIB=$HOME/projects/BoxLib

# default include directories
DIRS="$CASTRO/Source $CASTRO/Source/Src_3d $CASTRO/Networks/ignition_simple $CASTRO/EOS/gamma_law_general $CASTRO/EOS $CASTRO/constants $CASTRO/Exec/HCBubble $BOXLIB/Src/F_BaseLib $BOXLIB/Src/LinearSolvers/F_MG"

if (( $# < 1 )); then
  echo "Usage: $1 <target> <include-dirs>..."
else
  TARGET=$1
  shift
  $SCRIPT $TARGET $DIRS $@
fi
