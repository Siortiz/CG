import re 
from astropy.table import Table

# Función para ajustar minutos y segundos con un solo dígito y dos espacios
def adjust_single_digit_minutes_seconds(line):
    # Ajustar minutos con un solo dígito (dos espacios antes del minuto)
    line = re.sub(r'(?<!\S)\b(\d{1,2})\s\s(\d)\s(\d{1,2}\.\d+)', r'\1 0\2 \3', line)
    # Ajustar segundos con un solo dígito (dos espacios antes del segundo)
    line = re.sub(r'(?<!\S)\b(\d{1,2})\s(\d{1,2})\s\s(\d{1,2})\.(\d+)', r'\1 \2 0\3.\4', line)
    # Ajustar casos donde ambos minutos y segundos tienen un solo dígito (dos espacios antes de minutos y segundos)
    line = re.sub(r'(?<!\S)\b(\d{1,2})\s\s(\d)\s\s(\d{1,2})\.(\d+)', r'\1 0\2 0\3.\4', line) 
    return line
def replace_sexagesimal_spaces(line):
    #Expresión regular para valores sexagesimales
    pattern_ra = re.compile(r'(?<![A-Z\s])\b(\d{1,2})\s(\d{1,2})\s(\d{1,2}\.\d+)')
    pattern_dec = re.compile(r'(?<![A-Z\s])\b([+-]?\d{1,3})\s(\d{1,2})\s(\d{1,2}\.\d+)')
    
    # Reemplazar espacios con dos puntos
    line = pattern_ra.sub(r'\1:\2:\3', line)
    line = pattern_dec.sub(r'\1:\2:\3', line)

    return line

#Nombre del archivo original y el archivo temporal
original_filename = '/home/seba/Downloads/best.observations.idz'
temporal_filename = 'temp_best.observations.idz'

#Leer el archivo original y escribir en el archivo temporal

with open(original_filename, 'r') as original_file, open(temporal_filename, 'w') as temp_file:
    for line in original_file:
        line = adjust_single_digit_minutes_seconds(line)
        
        line = replace_sexagesimal_spaces(line)
        temp_file.write(line)

# Leer la tabla desde el archivo temporal especificando el delimitador y desactivando el adivinador
try:
    table = Table.read(temporal_filename, format='ascii', delimiter=' ', guess=False, fast_reader=False)
    # Obtener el número de columnas
    num_columns = len(table.colnames)

    # Mostrar el número de columnas
    print(f'Número de columnas: {num_columns}')

    # Opcional: mostrar los nombres de las columnas
    print('Nombres de las columnas:')
    print(table.colnames)
except Exception as e:
    print(f'Error al leer la tabla: {e}')
