#!/bin/bash

ANALYZE=./analyze.py
XML_DIR=../../examples/cns-smc/xml
OUT_DIR=output
inputs="advance-nomod advance-fused-nomod smc/advance-smc-modified smc/advance-smc-fused chemistry/drm19 chemistry/grimech30 chemistry/lidryer"

mkdir -p $OUT_DIR/smc
mkdir -p $OUT_DIR/chemistry
for input in $inputs; do
  $ANALYZE $XML_DIR/$input.xml $OUT_DIR/$input.tsv
done
