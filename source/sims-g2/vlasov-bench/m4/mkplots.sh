cmd="pgkyl -f m4-max-2d_field_0.bp -f m4-max-2d_field_2.bp interpolate -b mo -p 1 comp 2 fix --c1 0.5 hold plot"
echo $cmd
$cmd
