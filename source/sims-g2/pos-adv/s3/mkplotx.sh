cmd="pgkyl -f s3-2d-adv-dg_distf_0.bp -f s3-2d-adv-dg_distf_1.bp interpolate -p 1 -b ms fix --c0 0.5 hold plot --save"
echo $cmd
$cmd 
