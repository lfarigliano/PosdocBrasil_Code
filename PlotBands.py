#import pandas as pd
#import numpy as np
import math
#import matplotlib as mpl

import matplotlib.pyplot as plt

#import seaborn as sns
#import csv
#from scipy.spatial.transform import Rotation


message = """
This plot is version 1 of the code. Things to note:
- It can only be used for semiconductors or insulators.
- It must be used with the 'bands_1.sh' and 'bands_2.sh' protocols.
- It is only suitable for scanR2 calculations.
- The path in the Brillouin zone must be continuous.


FROM: Santos Morto, Grande Pedro, Rey de Copas 
"""

print(message)

print("Please provide the energy limits for the band plot:")

# Solicitar y capturar los valores de ymin y ymax desde el usuario
ymin = float(input("Enter ymin: "))
ymax = float(input("Enter ymax: "))


# parametros para las graficas
plt.rcParams['font.size'] = 26
plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams['xtick.major.size'] = 26
plt.rcParams['ytick.major.size'] = 26
plt.rcParams['xtick.major.width'] = 2.5
plt.rcParams['ytick.major.width'] = 2.5
plt.rcParams['axes.labelsize'] = 26  # Tamaño de la etiqueta de los ejes
plt.rcParams['xtick.labelsize'] = 26  # Tamaño de las etiquetas del eje x
plt.rcParams['ytick.labelsize'] = 26  # Tamaño de las etiquetas del eje y
plt.rcParams['legend.fontsize'] = 26  # Tamaño de letra de la leyenda


with open("EIGENVAL", "r") as archivo:
    lines_eigenval = archivo.readlines()
    electrons, nkpoints, nbands = lines_eigenval[5].split()

electrons=int(electrons)
nkpoints= int(nkpoints)
nbands=int(nbands)

with open("OUTCAR", "r") as file_outcar:
    lines_outcar = file_outcar.readlines()
    for i, line in enumerate(lines_outcar):
        if 'Following cartesian coordinates' in line:
            line_kpoints_outcar=i
            break
    for i, line in enumerate(lines_outcar):
        if 'Following reciprocal coordinates:' in line:
            line_kpoints_outcar_reciprocal=i
            break
        
kpoints_tmp=lines_outcar[line_kpoints_outcar+2:line_kpoints_outcar+2+nkpoints]
kpoints_tmp_reciprocal=lines_outcar[line_kpoints_outcar_reciprocal+2:line_kpoints_outcar_reciprocal+2+nkpoints]

for index, l in enumerate(kpoints_tmp):
    peso=(l.split())[3]
    if float(peso) == 0:
        indice_primero_cero = index
        break  # Salimos del bucle al encontrar el primer 0
for index, l in enumerate(kpoints_tmp_reciprocal):
    peso=(l.split())[3]
    if float(peso) == 0:
        indice_primero_cero_reciprocal = index
        break  # Salimos del bucle al encontrar el primer 0

x_kpoints=kpoints_tmp[indice_primero_cero:len(kpoints_tmp)]
x_kpoints_reciprocal=kpoints_tmp_reciprocal[indice_primero_cero_reciprocal:len(kpoints_tmp_reciprocal)]
x_kpoint_float=[]
x_kpoint_float_reciprocal=[]
for linea in x_kpoints:
    elementos = linea.split()[:3]  # Dividir la línea en elementos y tomar los primeros tres
    elementos_flotantes = [float(elemento) for elemento in elementos]  # Convertir a flotantes
    x_kpoint_float.append(elementos_flotantes)
for linea in x_kpoints_reciprocal:
    elementos = linea.split()[:3]  # Dividir la línea en elementos y tomar los primeros tres
    elementos_flotantes = [float(elemento) for elemento in elementos]  # Convertir a flotantes
    x_kpoint_float_reciprocal.append(elementos_flotantes)
distance=[]
for i in range(len(x_kpoint_float) - 1):
    punto1 = x_kpoint_float[i]
    punto2 = x_kpoint_float[i + 1]
    distancia = math.sqrt(sum((p1 - p2) ** 2 for p1, p2 in zip(punto1, punto2)))
    distance.append(distancia)


# Redondear los valores a 5 decimales
distancias_redondeadas = [round(dist, 3) for dist in distance]
diferencias_indices = [i for i in range(1, len(distancias_redondeadas)) 
                       if distancias_redondeadas[i] != distancias_redondeadas[i-1]]

discontinuidades=[]
for i in range(0, len(diferencias_indices), 2):
    discontinuidades.append(diferencias_indices[i])

suma_acumulada = []

# Inicializar una variable para llevar la cuenta de la suma acumulada
suma_total = 0
suma_acumulada.append(suma_total)
# Iterar sobre la lista de distancias y calcular la suma acumulada
for distancia in distance:
    suma_total += distancia
    suma_acumulada.append(suma_total)

head_eigenval=6
# buscar los puntos k y los indices
inicio=head_eigenval + 1
for i in range (inicio, len(lines_eigenval),nbands +2):
    peso=(lines_eigenval[i].split())[3]
    if float(peso) == 0:
        indice_primero_cero = i
        break  # Salimos del bucle al encontrar el primer 0


arc_band=lines_eigenval[indice_primero_cero:]
arc_band.append('\n') 

for j in range(1,nbands+1,1):
    globals()["band_up_"+str(j)]= []
    globals()["band_down_"+str(j)]= []
    globals()["band_up_oc_"+str(j)]= []
    globals()["band_down_oc_"+str(j)]= []

inicio=0
for l in range(int((len(arc_band))/(nbands+2))):
    final=(nbands+2)*(l+1)
    kp=arc_band[inicio:final]
    for k in kp[1:-1]:
        tmp=k.split()
        globals()["band_up_"+str(int(tmp[0]))].append(tmp[1])
        globals()["band_down_"+str(int(tmp[0]))].append(tmp[2])
        globals()["band_up_oc_"+str(int(tmp[0]))].append(tmp[3])
        globals()["band_down_oc_"+str(int(tmp[0]))].append(tmp[4])
    inicio=final

for j in range(1,nbands+1,1):
    globals()["band_up_" + str(j)] = [float(x) for x in globals()["band_up_" + str(j)]]
    globals()["band_down_" + str(j)] = [float(x) for x in globals()["band_down_" + str(j)]]
    globals()["band_up_oc_" + str(j)] = [float(x) for x in globals()["band_up_oc_" + str(j)]]
    globals()["band_down_oc_" + str(j)] = [float(x) for x in globals()["band_down_oc_" + str(j)]]

e_max_up=-100000000000000000000000
e_max_down=-100000000000000000000000
for j in range(1, nbands + 1):
    if globals()["band_up_oc_" + str(j)][0] == 1:
        max_energy = max(globals()["band_up_" + str(j)])
        if e_max_up <= max_energy:
            banda_up=j
            e_max_up= max_energy
    if globals()["band_down_oc_" + str(j)][0] == 1:
        max_energy = max(globals()["band_down_" + str(j)])
        if e_max_down <= max_energy:
            banda_down=j
            e_max_down= max_energy





if e_max_up >= e_max_down:
    mayor_valor = e_max_up
elif e_max_down > e_max_up:
    mayor_valor = e_max_down


e_min_up_cond=100000000000000000000000
e_min_down_cond=100000000000000000000000
banda_up_cond=0
banda_down_cond=0
for j in range(1, nbands + 1):
    if globals()["band_up_oc_" + str(j)][0] == 0:
        max_energy = min(globals()["band_up_" + str(j)])
        if e_min_up_cond >= max_energy:
            banda_up_cond=j
            e_min_up_cond= max_energy
    if globals()["band_down_oc_" + str(j)][0] == 0:
        max_energy = min(globals()["band_down_" + str(j)])
        if e_min_down_cond >= max_energy:
            banda_down_cond=j
            e_min_down_cond= max_energy

index_x_band_up_V = globals()["band_up_" + str(banda_up)].index(e_max_up)
index_x_band_down_V = globals()["band_down_" + str(banda_down)].index(e_max_down)

index_x_band_up_C = globals()["band_up_" + str(banda_up_cond)].index(e_min_up_cond)
index_x_band_down_C = globals()["band_down_" + str(banda_down_cond)].index(e_min_down_cond)



for j in range(1,nbands+1,1):
    globals()["band_up_" + str(j)] = [x - mayor_valor for x in globals()["band_up_" + str(j)]]
    globals()["band_down_" + str(j)] = [x - mayor_valor for x in globals()["band_down_" + str(j)]]


#letras_finales = ['G', 'M', 'K', 'G']

letras_finales=[]
letras=[0]
# Configura la figura
plt.figure(figsize=(14, 12))
for j in range(1, nbands + 1):
    plt.grid(True, color='gray', linestyle='--', alpha=0.3)
    plt.plot(suma_acumulada, globals()[f"band_up_{j}"], color='red', linewidth=2.5)
    plt.plot(suma_acumulada, globals()[f"band_down_{j}"], color='blue', linewidth=2.5)
    plt.ylim(ymin, ymax)
    plt.xlim(0, suma_acumulada[-1])
    for i in discontinuidades:
        plt.axvline(suma_acumulada[i], linestyle='--', lw=2.5, color='black')
    plt.axhline(0.0, linestyle='--', lw=2.5, color='black')
for i in discontinuidades:
    letras.append(i)
letras.append(len(suma_acumulada)-1)
# Crear etiquetas personalizadas para el eje x
xticks = [suma_acumulada[i] for i in letras]# + suma_acumulada[-4:]
#xtick_labels = ['0'] * len(discontinuidades) + letras_finales#
print(f"Need {len(xticks)} kpoint:")
for i in range(len(xticks)):
    k_value = input(f"Enter k value for point {i+1}:")
    if k_value.lower() in ["g", "G", "gama", "Gama", "gamma", "Gamma"]:
        k_value = '\u0393'
    letras_finales.append(k_value)
print('Your path is:',letras_finales)
## Aplicar las etiquetas personalizadas
plt.xticks(ticks=xticks, labels=letras_finales)
plt.ylabel('E - E$_{Fermi}$  //  eV')

# Mostrar el gráfico
plt.savefig('Bands.pdf')
#plt.show()


#Print Information
Band_Gap_up=e_min_up_cond-e_max_up
Band_Gap_down=e_min_down_cond-e_max_down
Band_Gap_total=min(e_min_up_cond,e_min_down_cond)-max(e_max_up,e_max_down)


print("+-------------------------- Summary ----------------------------+")
print("           Spin Channel:     <UP>       <DOWN>      <TOTAL>      ")
print("          Band Gap (eV):   ",round(Band_Gap_up,4),"    ",round(Band_Gap_down,4),"     ",round(Band_Gap_total,4))
print(" Eigenvalue of VBM (eV):   ",round(e_max_up,4),"    ",round(e_max_down,4),"     ",round(max(e_max_up,e_max_down),4))
print(" Eigenvalue of CBM (eV):   ",round(e_min_up_cond,4),"    ",round(e_min_down_cond,4),"     ",round(min(e_min_up_cond,e_min_down_cond),4))
if e_max_up >= e_max_down:
    banda_general_V=banda_up
    index_x_band_G_V=index_x_band_up_V
else:
    banda_general_V=banda_down
    index_x_band_G_V=index_x_band_down_V
print("      Band Index of VBM:   ",banda_up,"     ",banda_down,"        ",banda_general_V)
if e_min_up_cond <= e_min_down_cond:
    banda_general_C=banda_up_cond
    index_x_band_G_C=index_x_band_up_C
else:
    banda_general_C=banda_down_cond
    index_x_band_G_C=index_x_band_down_C
print("      Band Index of CBM:   ",banda_up_cond,"     ",banda_down_cond,"        ",banda_general_C)
print("+-----------Santos morto-------------+")
print("      Fermi Energy (eV):      ",round(max(e_max_up,e_max_down),4))
print("Location of VBM (UP):" ,x_kpoint_float_reciprocal[index_x_band_up_V])
print("Location of CBM (UP):" ,x_kpoint_float_reciprocal[index_x_band_up_C])
print("Location of VBM (Down):" ,x_kpoint_float_reciprocal[index_x_band_down_V])
print("Location of CBM (Down):" ,x_kpoint_float_reciprocal[index_x_band_down_C])
print("Location of VBM (Total):" ,x_kpoint_float_reciprocal[index_x_band_G_V])
print("Location of CBM (Total):" ,x_kpoint_float_reciprocal[index_x_band_G_C])
print("+----------------- Pedro, Daniel and Lucas ---------------------+")