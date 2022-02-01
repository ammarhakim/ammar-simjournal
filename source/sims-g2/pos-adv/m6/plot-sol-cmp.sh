cd ../s6
pgkyl -f s6-2d-adv-dg_deltaChange write -m txt
pgkyl -f s6-2d-adv-dg_rescaledCells write -m txt
cd ../m6
pgkyl -f m6-2d-adv-dg_density write -m txt
pgkyl -f m6-2d-adv-dg_deltaChange write -m txt
pgkyl -f m6-2d-adv-dg_rescaledCells write -m txt
python plot-sol-cmp.py
