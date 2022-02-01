cd ../s3
pgkyl -f s3-2d-adv-dg_absDist write -m txt
cd ../m3
pgkyl -f m3-2d-adv-dg_absDist write -m txt
python plot-sol-cmp.py
