#!/bin/bash

for i in {158..200}
do
    fName=`printf "s6-euler-kh_fluid_%05d.png" $i`
    cmd="pgkyl -f s6-euler-kh_fluid_$i.bp sel -c0 pl --fix-aspect --no-show --saveas $fName"
    echo $cmd
    $cmd
done
