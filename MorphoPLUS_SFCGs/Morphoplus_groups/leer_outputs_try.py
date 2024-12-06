from astropy.table import *
from astropy.table import Table
from astropy.io import ascii
import numpy as np
import pandas as pd
from astropy.io import fits
import matplotlib.pyplot as plt
import glob
import os
import re
from ejecutable import S
from astropy.table import Table

Filtros=np.array(['R','F378','F395','F410','F430','F515','F660','F861','G','I','Z','U']) #filtros de splus
Fil_name=np.array(['R','J0378','J0395','J0410','J0430','J0515','J0660','J0861','G','I','Z','U']) #filtros de splus
def leer_header(NGAL,HEADER):
    data = {}
    data['CHI2NU'] = HEADER['CHI2NU']
    for i in range(12):
        xc=HEADER['%i_XC_%s'%(NGAL,Fil_name[i])].split( )
        yc=HEADER['%i_YC_%s'%(NGAL,Fil_name[i])].split( )
        R=HEADER['%i_RE_%s'%(NGAL,Fil_name[i])].split( )
        m=HEADER['%i_MAG_%s'%(NGAL,Fil_name[i])].split( )
        N=HEADER['%i_n_%s'%(NGAL,Fil_name[i])].split( )
        ar=HEADER['%i_AR_%s'%(NGAL,Fil_name[i])].split( )
        PA=HEADER['%i_PA_%s'%(NGAL,Fil_name[i])].split( )
		#print(HEADER['%i_XC_%s'%(NGAL,Fil_name[i])])
        data['XC_%s' % Fil_name[i]] = xc[0]
        data['e_XC_%s' % Fil_name[i]] = xc[2]
        data['YC_%s' % Fil_name[i]] = yc[0]
        data['e_YC_%s' % Fil_name[i]] = yc[2]
        data['MAG_%s' % Fil_name[i]] = m[0]
        data['e_MAG_%s' % Fil_name[i]] = m[2]
        data['RE_%s' % Fil_name[i]] = R[0]
        data['e_RE_%s' % Fil_name[i]] = R[2]
        data['n_%s' % Fil_name[i]] = N[0]
        data['e_n_%s' % Fil_name[i]] = N[2]
        data['AR_%s' % Fil_name[i]] = ar[0]
        data['e_AR_%s' % Fil_name[i]] = ar[2]
        data['PA_%s' % Fil_name[i]] = PA[0]
        data['e_PA_%s' % Fil_name[i]] = PA[2]
        '''
        data.append(xc[0]) #crea una fila con los paremetros
        data.append(xc[2])
        data.append(yc[0])
        data.append(yc[2])
        data.append(m[0])
        data.append(m[2])
        data.append(R[0])
        data.append(R[2])
        data.append(N[0])
        data.append(N[2])
        data.append(ar[0])
        data.append(ar[2])
        data.append(PA[0])
        data.append(PA[2])'''
    return data

def grafico(f,grupo,field):
	fig, axs = plt.subplots(3, 12,figsize=(20,8),  sharex=True, sharey=True)
	vmax=[0.03,0.03,0.03,.1,.1,0.7,0.5,0.7,0.5,0.7,0.7,0.7]
	for j  in range(12):
		axs[0, j].imshow(f[j].data,cmap='gray', vmin=-vmax[j],vmax=vmax[j])
		axs[0, j].set_title(f'Filter {Bands[j]}')
		axs[1, j].imshow(f[12+j].data,cmap='gray',vmin=-vmax[j],vmax=vmax[j])
		axs[2, j].imshow(f[24+j].data,cmap='gray',vmin=-vmax[j],vmax=vmax[j])
		axs[0,0].text( 60, 140, '5.5\"', fontsize=10,color='blue') 
		axs[0,0].plot([60, 80], [150, 150], 'b-', lw=3)
	plt.savefig('Out_Img/%s_%s.svg'%(grupo,field),format='svg', dpi=1200)
	#plt.show()
	plt.close()
	return
Tabla=[]
header_names = ['Group', 'Gal', 'ID', 'CHI2NU']


for i in range(12):
    header_names.append('XC_%s'%(Fil_name[i]))
    header_names.append('e_XC_%s'%(Fil_name[i]))
    header_names.append('YC_%s'%(Fil_name[i]))
    header_names.append('e_YC_%s'%(Fil_name[i]))
    header_names.append('MAG_%s'%(Fil_name[i]))
    header_names.append('e_MAG_%s'%(Fil_name[i]))
    header_names.append('RE_%s'%(Fil_name[i]))
    header_names.append('e_RE_%s'%(Fil_name[i]))
    header_names.append('n_%s'%(Fil_name[i]))
    header_names.append('e_n_%s'%(Fil_name[i]))
    header_names.append('AR_%s'%(Fil_name[i]))
    header_names.append('e_AR_%s'%(Fil_name[i]))
    header_names.append('PA_%s'%(Fil_name[i]))
    header_names.append('e_PA_%s'%(Fil_name[i]))

#$S=Table.read('/home/seba/Documents/CG/Observaciones_SPLUS/Input_MorphoPlus.csv')
Datos_S= S.group_by('Field')


Fields=Datos_S.groups.keys 
Bands=np.array(['U','F378','F395','F410','F430','G','F515','R','F660','I','F861','Z'])
no_ajustados = []
for f in Fields['Field']:
    SF=S[S['Field']==f]
    Datos_SF= SF.group_by('Group')
    GS=Datos_SF.groups.keys #grupos
    for g in GS['Group']:
        try:
            mask=Datos_SF.groups.keys['Group'] ==  g #'cCGs-4007' #Gt['Groups'][i]
            Tablef=Datos_SF.groups[mask]
            fi=fits.open("galfitm_%s_%s.fits"%(g,f))
        
            for i in range(len(Tablef)):
                #fi=fits.open("galfitm_%s_%s.fits"%(g,f))
                gal = Tablef[i]['Gal']
                ID = Tablef[i]['ID']
                header_data = leer_header(i+1,fi[13].header)
                #print(type(g))
                header_data['Group'] = g
                header_data['Gal'] = gal
                header_data['ID'] = ID
                Tabla.append(header_data)
		        #grafico(fi,g,f)
        except:
            no_ajustados.append(g)
print(f"Los grupos no ajustados son: {no_ajustados}")


#header_names = ['Group', 'Gal', 'ID', 'CHI2NU']


'''for i in range(12):
	header_names.append('XC_%s'%(Fil_name[i]))
	header_names.append('e_XC_%s'%(Fil_name[i]))
	header_names.append('YC_%s'%(Fil_name[i]))
	header_names.append('e_YC_%s'%(Fil_name[i]))
	header_names.append('MAG_%s'%(Fil_name[i]))
	header_names.append('e_MAG_%s'%(Fil_name[i]))
	header_names.append('RE_%s'%(Fil_name[i]))
	header_names.append('e_RE_%s'%(Fil_name[i]))
	header_names.append('n_%s'%(Fil_name[i]))
	header_names.append('e_n_%s'%(Fil_name[i]))
	header_names.append('AR_%s'%(Fil_name[i]))
	header_names.append('e_AR_%s'%(Fil_name[i]))
	header_names.append('PA_%s'%(Fil_name[i]))
	header_names.append('e_PA_%s'%(Fil_name[i]))'''
	 

				
table = Table(rows=Tabla,
              names=header_names)

						
#Datos_S= S.group_by('Group')
#GS=Datos_S.groups.keys
#g = int(GS['Group'])
ascii.write(table, f'Catalogos_try/GalfitM_output_final_ahorasi.csv', format='csv', fast_writer=False, overwrite=True) #guarda la tabla completa con los 
		
		
		
	
		
