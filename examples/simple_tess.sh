#!/bin/bash

echo "Calculate the potential, gravity, and gradient tensor of a tesseroid model"
echo '# Test tesseroid model file
-5 5 -10 10 0 -50000 1000
-12 -16 -12 -16 0 -30000 -1000' > model.txt
echo "The model file:"
cat model.txt
echo "Calculating..."
time tessgrd -r-20/20/-20/20 -b50/50 -z250e03 | \
tesspot model.txt | \
tessgx model.txt | tessgy model.txt | tessgz model.txt | \
tessgxx model.txt  | tessgxy model.txt  | tessgxz model.txt  | \
tessgyy model.txt  | tessgyz model.txt  | tessgzz model.txt  -v > output.txt
python plot.py output.txt 50 50
