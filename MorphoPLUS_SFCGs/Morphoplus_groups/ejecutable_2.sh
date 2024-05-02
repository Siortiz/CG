#python Recortar_groups.py
chmod 777 dopsfex_mask.sh
./dopsfex_mask.sh
python inputs_galfitm.py
chmod 777 inputs/galfit_267_MC0137.input
./galfitm-1.4.4-linux-x86_64 inputs/galfit_267_MC0137.input
chmod 777 inputs/galfit_230_SPLUS-s22s19.input
./galfitm-1.4.4-linux-x86_64 inputs/galfit_230_SPLUS-s22s19.input
chmod 777 inputs/galfit_236_SPLUS-s36s36.input
./galfitm-1.4.4-linux-x86_64 inputs/galfit_236_SPLUS-s36s36.input
python leer_outputs.py
