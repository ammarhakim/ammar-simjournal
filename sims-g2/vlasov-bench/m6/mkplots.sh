cmd="pgkyl -f m6-max-2d_field_0.bp -f m6-max-2d_field_2.bp interpolate -b mo -p 2 comp 2 fix --c1 0.5 hold plot"
echo $cmd
$cmd
