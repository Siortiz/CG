python Recortar.py
chmod 777 dopsfex_mask_SPLUS-n12s30.sh
./dopsfex_mask_SPLUS-n12s30.sh
chmod 777 galfit_1650_550_SPLUS-n12s30.input
chmod 777 galfit_1650_1650_SPLUS-n12s30.input
chmod 777 galfit_1650_2750_SPLUS-n12s30.input
chmod 777 galfit_1650_4950_SPLUS-n12s30.input
chmod 777 galfit_2750_1650_SPLUS-n12s30.input
chmod 777 galfit_2750_3850_SPLUS-n12s30.input
chmod 777 dopsfex_mask_gal.sh
./dopsfex_mask_gal.sh
python mascara.py
./galfitm-1.4.4-linux-x86_64 galfit_1650_550_SPLUS-n12s30.input
./galfitm-1.4.4-linux-x86_64 galfit_1650_1650_SPLUS-n12s30.input
./galfitm-1.4.4-linux-x86_64 galfit_1650_2750_SPLUS-n12s30.input
./galfitm-1.4.4-linux-x86_64 galfit_1650_4950_SPLUS-n12s30.input
./galfitm-1.4.4-linux-x86_64 galfit_2750_1650_SPLUS-n12s30.input
./galfitm-1.4.4-linux-x86_64 galfit_2750_3850_SPLUS-n12s30.input
chmod 777 ejecutable_gal.sh
./ejecutable_gal.sh
python leer_header_output.py
