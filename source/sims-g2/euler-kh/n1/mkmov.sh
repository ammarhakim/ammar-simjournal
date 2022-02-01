#!/bin/bash

for i in {0..200}
do
    fName=`printf "n1-euler-kh_fluid_%05d.png" $i`
    cmd="pgkyl -f n1-euler-kh_fluid_$i.bp sel -c0 pl --fix-aspect --no-show --saveas $fName"
    echo $cmd
    $cmd
done
