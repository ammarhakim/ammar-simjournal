#!/bin/bash
for i in {190..190}
do
    fName=`printf "s5-euler-kh_chi-spect_%05d.bp" $i`
    cmd="pgkyl -f s5-euler-kh_fluid_$i.bp euler -v vel -g 1.4 sel -c0:2 ev \"f0 curl sq\" integrate 0 fft -p write -f $fName"
    echo $cmd
    $cmd
done
