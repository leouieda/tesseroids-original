#!/bin/bash

echo "Calculate the potential, gravity, and gradient tensor of a prism model"
echo '# Test prism model file
2000 5000 2000 15000 0 5000 1000
10000 18000 10000 18000 0 5000 -1000' > model.txt
echo "The model file:"
cat model.txt
echo "Calculating..."
time tessgrd -r0/20000/0/20000 -b50/50 -z1000 | \
prismpot model.txt | \
prismgx model.txt | prismgy model.txt | prismgz model.txt | \
prismgxx model.txt  | prismgxy model.txt  | prismgxz model.txt  | \
prismgyy model.txt  | prismgyz model.txt  | prismgzz model.txt > output.txt
python plot.py output.txt 50 50
