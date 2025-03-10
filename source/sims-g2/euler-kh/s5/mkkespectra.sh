#!/bin/bash

for i in {190..200}
do
    fName=`printf "s5-euler-kh_ke-spect_%05d.bp" $i`
    cmd="pgkyl -f s5-euler-kh_fluid_$i.bp euler -v ke -g 1.4 integrate 0 fft -p write -f $fName"
    echo $cmd
    $cmd
done
