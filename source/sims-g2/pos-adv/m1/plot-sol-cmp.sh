cd ../s1
pgkyl -f s1-2d-adv-dg_absDist write -m txt
cd ../m1
pgkyl -f m1-2d-adv-dg_absDist write -m txt
python plot-sol-cmp.py
