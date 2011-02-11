#!/bin/bash

for i in {1..1}
do
    time 2> time_normal.txt echo "0 0 250000" | tessgzz model.txt -a  > normal_data.txt
done