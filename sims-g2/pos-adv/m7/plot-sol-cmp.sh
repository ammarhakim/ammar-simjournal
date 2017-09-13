cd ../s6
pgkyl -f s6-2d-adv-dg_deltaChange write -m txt
cd ../m7
pgkyl -f m7-2d-adv-dg_density write -m txt
pgkyl -f m7-2d-adv-dg_deltaChange write -m txt
python plot-sol-cmp.py
