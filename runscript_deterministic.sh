#!/bin/bash

for i in {30..150..1}
do
for j in `seq 0.1 0.005 0.2`
    do 

   ./main3.exe 1000 $i 1000-$i 5000 1.1 1.0 0.5 0.1 $j 1.0 1000
done
done