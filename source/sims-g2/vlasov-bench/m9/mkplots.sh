cmd="pgkyl -f m9-max-3d_field_0.bp -f m9-max-3d_field_2.bp  -c interpolate -b ms -p 2 comp 2 fix --c1 0.5 --c2 0.5 hold plot"
echo $cmd
$cmd
