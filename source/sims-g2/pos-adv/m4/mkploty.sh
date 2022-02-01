cmd="pgkyl -f m4-2d-adv-dg_distf_0.bp -f m4-2d-adv-dg_distf_1.bp interpolate -p 1 -b ms fix --c1 0.5 hold plot --save"
echo $cmd
$cmd
