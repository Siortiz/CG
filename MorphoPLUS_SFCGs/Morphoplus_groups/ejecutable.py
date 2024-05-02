
from astropy.table import *
from astropy.io import ascii
import numpy as np
import csv

S=Table.read('/home/seba/Documents/CG/Observaciones_SPLUS/Galaxies/Galaxies_SPLUS_270.csv')

size=500
#archivo_csv = "Catalogos/g_S.csv"

Datos_S= S.group_by('Field')
Fields=Datos_S.groups.keys 
Data=['python Recortar_groups.py']
Data.append('chmod 777 dopsfex_mask.sh')
Data.append('./dopsfex_mask.sh')
Data.append('python inputs_galfitm.py')

def ejecutable(Fields):
	for f in Fields['Field']:	
		SF=S[S['Field']==f]
		Datos_SF= SF.group_by('Group')
		GS=Datos_SF.groups.keys #grupos
		for g in GS['Group']:
			Data.append('chmod 777 inputs/galfit_%s_%s.input'%(g,f))
			Data.append('./galfitm-1.4.4-linux-x86_64 inputs/galfit_%s_%s.input'%(g,f))
	Data.append('python leer_outputs.py')
	fic = open('ejecutable.sh','w')
	for line in Data:
		print(line, file=fic)
	fic.close()
	return


ejecutable(Fields)
