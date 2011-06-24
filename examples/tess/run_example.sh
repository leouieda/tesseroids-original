#!/bin/bash

tessgrd -r-5/5/-5/5 -b100/100 -z250e03 | tessgz model.txt  | \
tessgxx model.txt  | tessgxy model.txt  | tessgxz model.txt  | \
tessgyy model.txt  | tessgyz model.txt  | tessgzz model.txt  -v -lgzz.log > output.txt
