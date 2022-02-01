cd ../s5
pgkyl -f s5-2d-adv-dg_absDist write -m txt
cd ../m5
pgkyl -f m5-2d-adv-dg_absDist write -m txt
python plot-sol-cmp.py
