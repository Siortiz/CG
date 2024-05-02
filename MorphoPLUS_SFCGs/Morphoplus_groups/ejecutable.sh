python Recortar_groups.py
chmod 777 dopsfex_mask.sh
./dopsfex_mask.sh
python inputs_galfitm.py
chmod 777 inputs/galfit_270_SPLUS-s29s37.input
./galfitm-1.4.4-linux-x86_64 inputs/galfit_270_SPLUS-s29s37.input
python leer_outputs.py
