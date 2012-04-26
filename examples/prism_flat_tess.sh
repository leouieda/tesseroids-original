#!/bin/bash

echo "Calculate the potential, gravity, and gradient tensor of a prism model generated from a flattened tesseroid model"
echo '# Test tesseroid model file
-5 5 -10 10 0 -50000 1000
-12 -16 -12 -16 0 -30000 -1000' | tess2prism --flatten > model.txt
echo "The model file:"
cat model.txt
echo "Calculating..."
time tessgrd -r-20/20/-20/20 -b50/50 -z250e03 | \
prismpot model.txt | \
prismgx model.txt | prismgy model.txt | prismgz model.txt | \
prismgxx model.txt  | prismgxy model.txt  | prismgxz model.txt  | \
prismgyy model.txt  | prismgyz model.txt  | prismgzz model.txt  -v > output.txt
python plot.py output.txt 50 50
