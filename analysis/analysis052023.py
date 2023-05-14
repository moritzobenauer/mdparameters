import MDAnalysis as mda
from MDAnalysis.analysis.dihedrals import Dihedral
from MDAnalysis.analysis.dihedrals import Ramachandran
from MDAnalysis.core.topologyattrs import Atomnames
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from shapely.geometry import Polygon, Point
import sys
import csv
import os
import datetime
from MDAnalysis.analysis import distances
from MDAnalysis.analysis.waterdynamics import MeanSquareDisplacement as MSD
import copy

def smooth_diff(RG, kernel_size):

    data = np.absolute(np.gradient(RG, 1))
    kernel = np.ones(kernel_size) / kernel_size
    smooth_diff = np.convolve(data, kernel, mode="same")

    return smooth_diff

def smooth_data(data, kernel_size):

    kernel = np.ones(kernel_size) / kernel_size
    smooth_data = np.convolve(data, kernel, mode="same")

    return smooth_data

name = 'md_0_1_noPBC'
u = mda.Universe(f"{name}.gro", f"{name}.xtc")
chain_length = int(len(u.atoms.select_atoms('resname PEY'))/3)
print(chain_length)

counter, sheet_counter, helix_counter, hairpin_counter, random_counter = 0, 0, 0, 0, 0

# Radius of Gyration

t, RG = [], []

for ts in u.trajectory:
    ca = u.select_atoms("all")
    rad = ca.radius_of_gyration()
    t.append(ts.frame * 0.01)
    RG.append(rad)
delta_rg = np.asarray(smooth_diff(RG, 5))
norm_factor = np.amax(delta_rg)
delta_rg = delta_rg / norm_factor


A = u.select_atoms("resname PHE and name CA")
B = u.select_atoms("resname PHX and name CA")
out = mda.analysis.distances.dist(A, B)
print("Average PHE(CA) - PHX(CA) Distance")
print(np.average(out[2]))

C = u.select_atoms("resname HISD and name CA")
D = u.select_atoms("resname HISX and name CA")
out = mda.analysis.distances.dist(C, D)
print("Average HISD(CA) - HISX(CA) Distance")
print(np.average(out[2]))

ags = u.select_atoms("resid 2-6")
R = Ramachandran(ags, c_name='C', n_name='N', ca_name='CA',check_protein=False).run()
angles = R.angles

phe, his, secstructure, hairpin, helix, coil = [], [], [], [], [], []
time_c = np.arange(0, len(u.trajectory), 1)
time1 = time_c / 100

for ts in u.trajectory:
    phe.append(np.average(mda.analysis.distances.dist(A, B)[2]))
    his.append(np.average(mda.analysis.distances.dist(C, D)[2]))


for ts in angles:
    for item in ts:
        secstructure_temp = []
        random_temp = []
        if (-180 < item[0] < -60) and (60 < item[1] < 180):
            sheet_counter += 1
            counter += 1
        elif (60< item[0] < 120) and (60 < item[1] < 180):
            hairpin_counter += 1
            counter +=1
        elif (-180 < item[0] < -60) and (-60 < item[1] < 60):
            helix_counter +=1
            counter +=1
        else:
            random_counter += 1
            counter += 1
        secstructure_temp.append(sheet_counter/counter)
        random_temp.append(random_counter/counter)

    secstructure.append(np.average(secstructure_temp))
    coil.append(np.average(random_temp))

print(f"PEO Chain Length: {chain_length}")
print(f"Beta Structure: {round(sheet_counter/counter * 100, 1)}%")
print(f"Helix Structure: {round(helix_counter/counter * 100, 1)}%")
print(f"Hairpin Structure: {round(hairpin_counter/counter * 100, 1)}%")
print(f"Random Coil Structure: {round(random_counter/counter * 100, 1)}%")

sheet = np.round(copy.deepcopy(sheet_counter/counter), 2)
random = np.round(copy.deepcopy(random_counter/counter), 2)

x_bin, y_bin = [], []
for timestep in angles:
    for coordinate in timestep:
        x_bin.append(coordinate[0])
        y_bin.append(coordinate[1])


select = "name CA"
MSD_analysis = MSD(u, select, 0, 1000, 20)
MSD_analysis.run()
time = 0
for msd in MSD_analysis.timeseries:
      print("{time} {msd}".format(time=time, msd=msd))
      time += 1


fig, ((ax1, ax3), (ax4, ax5)) = plt.subplots(2, 2)
ax2 = ax1.twinx()
y=np.vstack([secstructure, coil])
ax1.stackplot(time1, y, colors=["#7f8c8d", "#ecf0f1"], labels=["Beta Sheet Structure", "Rancom Coil Structure"])
ax1.set(ylim=(0.5, 1), yticks=np.arange(0.5, 1, 0.1), xlim=(0, 10))
ax2.plot(time1, phe, label="PHE-PHX", color="#e67e22")
ax2.plot(time1, his, label="HISD-HISX", color="#2980b9")
ax2.legend(loc=0)
ax1.legend(loc="upper center")
ax1.set_xlabel("Time / ns")
ax1.set_ylabel("Secondary Structure")
ax2.set_ylabel(r"<C$\alpha$-C$\alpha$-Distances> / nm")
plt.title(f"Seondary Structure Analysis with N(PEO) = {chain_length} \n <Beta Sheet>: {sheet}, <Random Coil>: {random}", pad=20)
fig.set_size_inches(14, 10)
export = list(zip(time1, phe, his, secstructure, coil, RG, delta_rg))
df = pd.DataFrame(export,columns =['Time / ns', 'C-C-Distance (Phe)', 'C-C-Distance (His)', 'BSS', 'RCS', 'Radius of Gyration', 'Delta RG'])
df.to_csv('Timeseries_Analysis.csv', sep=',')
export=list(zip(x_bin, y_bin))
df = pd.DataFrame(export, columns=["Phi", "Psi"])
df.to_csv('Ramachandran_Analysis.csv', sep=',')
ax3.hist2d(x_bin, y_bin, bins=(120, 120), range=([-180, 180], [-180, 180]), density=True, cmap=plt.cm.jet)
ax3.set_xticks([-180, -120, -60, 0, 60, 120, 180])
ax3.set_yticks([-180, -120, -60, 0, 60, 120, 180])
ax3.grid(visible=True, which="major")
#plt.title(f'Ramachandran 2D Binning of FHFHF-PEO{chain_length}-FHFHF')
ax3.set_xlabel(r'$\Phi$')
ax3.set_ylabel(r'$\Psi$')
ax4.scatter(t, RG, marker='.',color="#20bf6b", linewidths=0.5, alpha=0.8, label="Radius of Gyration")
ax4.set_xlabel("Time / ns")
ax4.set_ylabel("Distance / nm")
ax4.legend()
ax6 = ax4.twinx()
ax6.plot(t, delta_rg, linewidth=1, alpha=0.8, color='#4b6584', label="Smoothed Derivative")
ax5.loglog(range(0,time),MSD_analysis.timeseries)
plt.show()




