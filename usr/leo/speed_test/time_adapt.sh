#!/bin/bash

echo "" > time_adapt.txt

for i in {1..3}
do
    (time tessgzz model.txt -a -v < grid.txt > adapt_data.txt) 2>> time_adapt.txt
done