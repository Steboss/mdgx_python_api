#!/bin/bash


for f in * ; do 
cd $f
sbatch structure_$f.sh
cd ../
done

