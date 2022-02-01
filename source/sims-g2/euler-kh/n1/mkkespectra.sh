#!/bin/bash

for i in {190..200}
do
    fName=`printf "n1-euler-kh_ke-spect_%05d.bp" $i`
    cmd="pgkyl -f n1-euler-kh_fluid_$i.bp euler -v ke -g 1.4 integrate 0 fft -p write -f $fName"
    echo $cmd
    $cmd
done
