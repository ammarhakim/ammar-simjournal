cd ../s4
pgkyl -f s4-2d-adv-dg_deltaChange write -m txt
pgkyl -f s4-2d-adv-dg_rescaledCells write -m txt
cd ../m4
pgkyl -f m4-2d-adv-dg_density write -m txt
pgkyl -f m4-2d-adv-dg_deltaChange write -m txt
pgkyl -f m4-2d-adv-dg_rescaledCells write -m txt
python plot-sol-cmp.py
