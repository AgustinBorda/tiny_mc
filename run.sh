#!/bin/bash
touch out.csv
for i in $(seq 1 $1) ; do
	echo seq 1 $1
       	./tiny_mc | grep "photons" | grep -o [0-9]*\\.[0-9]* >> out.csv
done
