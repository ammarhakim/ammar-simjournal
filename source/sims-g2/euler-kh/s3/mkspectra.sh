#!/bin/bash

for i in {190..200}
do
    fName=`printf "s3-euler-kh_spect_%05d.bp" $i`
    cmd="pgkyl -f s3-euler-kh_fluid_$i.bp sel -c0 integrate 0 fft -p write -f $fName"
    echo $cmd
    $cmd
done
