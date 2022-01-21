#!/bin/bash

#yield: why hadd != sum one by one
for f in $(ls tree_dataSB_OO_*root)
do
  root -q -l yield_dataBin.cxx\(\"${f}\"\)
done
