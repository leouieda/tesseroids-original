#!/bin/bash

echo "" > time_normal.txt

for i in {1..3}
do
    (time tessgzz model.txt -v $1 < grid.txt > normal_data.txt) 2>> time_normal.txt
done