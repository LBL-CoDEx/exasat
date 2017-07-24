#!/bin/bash

ORIG_DIR=$(pwd)
: ${EXASAT_TOP="${ORIG_DIR}/../.."}
EXE="${EXASAT_TOP}/src/genCodeDescript"

: ${HPGMG_TOP="${ORIG_DIR}/../../../lbl/hpgmg"}
cd ${HPGMG_TOP}/finite-volume/source

for op_type in 7pt 27pt fv2 fv4; do
  ${EXE} -I./ -DUSE_GSRB -DGSRB_BRANCH  -DUSE_HELMHOLTZ hpgmg-fv.c mg.c operators.${op_type}.c | tee ${ORIG_DIR}/hpgmg.${op_type}.br.xml
  ${EXE} -I./ -DUSE_GSRB -DGSRB_FP      -DUSE_HELMHOLTZ hpgmg-fv.c mg.c operators.${op_type}.c | tee ${ORIG_DIR}/hpgmg.${op_type}.fp.xml
  ${EXE} -I./ -DUSE_GSRB -DGSRB_STRIDE2 -DUSE_HELMHOLTZ hpgmg-fv.c mg.c operators.${op_type}.c | tee ${ORIG_DIR}/hpgmg.${op_type}.s2.xml
done
