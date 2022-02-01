cmd="pgkyl -f m8-max-3d_field_0.bp -f m8-max-3d_field_2.bp  -c interpolate -b mo -p 1 comp 2 fix --c1 0.5 --c2 0.5 hold plot"
echo $cmd
$cmd
