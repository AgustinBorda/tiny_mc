#!/bin/bash
touch out.csv
for i in {0..$1} do
 ./tiny_mc | grep "photons" | grep -o [0-9]*\\.[0-9]* >> out.csv
done
