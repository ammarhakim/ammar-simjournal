# density
for i in {0..15}
do
    fname=`printf s448-blast-density_%05d.png $i`
    cmd="pgkyl -f s448-blast-2d_q_$i.h5 comp 0 plot --no-show --fixed-axis --saveas $fname"
    echo $cmd
    $cmd
done

# pressure
for i in {0..15}
do
    fname=`printf s448-blast-pressure_%05d.png $i`
    cmd="pgkyl -f s448-blast-2d_q_$i.h5 euler -v pressure -g 1.666666667 plot --no-show --fixed-axis --saveas $fname"
    echo $cmd
    $cmd
done
	  
