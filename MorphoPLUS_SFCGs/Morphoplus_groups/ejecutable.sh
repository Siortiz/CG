python Recortar_groups.py
chmod 777 dopsfex_mask.sh
./dopsfex_mask.sh
python inputs_galfitm.py
chmod 777 galfit_234_MC0011.input
./galfitm-1.4.4-linux-x86_64 inputs/galfit_234_MC0011.input
python leer_outputs.py
