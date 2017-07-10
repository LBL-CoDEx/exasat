#!/bin/bash

ANALYZE=./analyze.py
XML_DIR=../../examples/cns-smc/xml_old
OUT_DIR=output
inputs="advance-flat advance-smc-modified drm19 grimech30 hai lidryer prf_ethanol"

export old=1
mkdir -p $OUT_DIR
for input in $inputs; do
  $ANALYZE $XML_DIR/$input.xml $OUT_DIR/$input.tsv
done
