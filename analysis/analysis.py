from functions import *

now = datetime.datetime.now()
timenow = str(now.hour) + str(now.minute)
print(timenow)

name = sys.argv[1]
peo = int(sys.argv[2])

prot_start = peo + 11
prot_end = prot_start + 4


topol = 'md_0_1_noPBC'
traj = 'md_0_1_noPBC'
u = mda.Universe(f"{name}.gro", f"{name}.xtc")
print(u)
print(len(u.trajectory))
dir = f'Analysis_{topol}_{timenow}'
os.mkdir(dir)
print(u.select_atoms(f"resid 2-6 or resid {prot_start}-{prot_end}").n_atoms)

# Radius of Gyration

t, RG = [], []

for ts in u.trajectory:
    ca = u.select_atoms("all")
    rad = ca.radius_of_gyration()
    t.append(ts.frame * 0.01)
    RG.append(rad)


delta_rg = smooth_diff(RG, 5)
plt.plot(t, RG)
plt.plot(t, delta_rg)
plt.savefig(f'{dir}/ROGAnalysis.png')

export = {'Timestep': t, 'R_G': RG, 'Delta R_G': delta_rg}
df = pd.DataFrame(export)
df.to_csv(f'{dir}/ROGAnalysis.csv')

# Ramachandran Advanced Analysis
# Import Beta Sheet Reference Data and convert it into Polygon

df = pd.read_csv("betasheet.csv", sep=",", header=None)
vertices = df.values.tolist()
betasheet = Polygon(vertices)

dg = pd.read_csv("helix.csv", sep=",", header=None)
vertices2 = dg.values.tolist()
helix = Polygon(vertices2)

def CheckBetaStructure(point):

    if point.within(betasheet) == True:
        return True
    else:
        return False

def CheckHelixStructure(point):

    if point.within(helix) == True:
        return True
    else:
        return False

def Chi(angles):
    timeframe = 0
    chi = []
    global time
    time = []
    for timestep in angles:
        timeframe += 1
        time.append(timeframe * 0.01)
        overall = len(timestep)
        counter = 0
        for coordinate in timestep:
            geopoint = Point(coordinate[0], coordinate[1])
            if CheckBetaStructure(geopoint) == True:
                counter += 1
            else:
                pass
        chi_temp = counter / overall
    
        chi.append(chi_temp)
    return chi

def Alpha(angles):
    timeframe = 0
    alpha = []
    global time
    time = []
    for timestep in angles:
        timeframe += 1
        time.append(timeframe * 0.01)
        overall = len(timestep)
        counter = 0
        for coordinate in timestep:
            geopoint = Point(coordinate[0], coordinate[1])
            if CheckHelixStructure(geopoint) == True:
                counter += 1
            else:
                pass
        alpha_temp = counter / overall
    
        alpha.append(alpha_temp)
    return alpha

def list_average(list):
    a = sum(list) / len(list)
    b = np.full(shape=1001, fill_value=a)
    return a, b


# Check if backbone for Res 2-6 is beta sheet
ags = u.select_atoms("resid 2-6 ")
R = Ramachandran(ags, c_name='C', n_name='N', ca_name='CA',
                 check_protein=False).run()
angles = R.angles
chi1 = Chi(angles)
alpha1 = Alpha(angles)
fig, ax = plt.subplots(figsize=plt.figaspect(1))

# Check if backbone for Res 21-25 is beta sheet
ags = u.select_atoms(f"resid {prot_start}-{prot_end} ")
R = Ramachandran(ags, c_name='N', n_name='C', ca_name='CA',
                 check_protein=False).run()
angles = R.angles
chi2 = Chi(angles)
alpha2 = Alpha(angles)

chi_avg = [(g+h) / 2 for g, h in zip(chi1, chi2)]
smooth_chi = smooth_data(chi_avg, 5)
average = list_average(chi_avg)[0]

alpha_avg = [(g+h) / 2 for g, h in zip(alpha1, alpha2)]
smooth_alpha = smooth_data(alpha_avg, 5)
average_alpha = list_average(alpha_avg)[0]

plt.plot(time, smooth_data(chi_avg, 5))
plt.plot(time, list_average(chi_avg)[1])
plt.xlabel("Time / ns")
plt.ylabel(r" $\chi$")
plt.title(f"Secondary Structure Analysis of Ac-FHFHF-PEG({peo})-FHFHF-Ac")
plt.savefig(f'{dir}/ChiAnalysis.png')

plt.plot(time, smooth_data(alpha_avg, 5))
plt.plot(time, list_average(alpha_avg)[1])
plt.xlabel("Time / ns")
plt.ylabel(r" $\alpha$")
plt.title(f"Secondary Structure Analysis of Ac-FHFHF-PEG({peo})-FHFHF-Ac")
plt.savefig(f'{dir}/AlphaAnalysis.png')



export = {'Time (ns)': time, 'Raw Chi': chi_avg, 'Smooth Chi': smooth_chi, 'Average Chi': average,
            'Raw Alpha': alpha_avg, 'Smooth Alpha': smooth_alpha, 'Average Alpha': average_alpha}
df = pd.DataFrame(export)
df.to_csv(f'{dir}/AdvancedRamachandranAnalysis.csv')

# Standard Ramachandran Analysis

ags = u.select_atoms("resid 2-6 ")
R = Ramachandran(ags, c_name='C', n_name='N', ca_name='CA',
                 check_protein=False).run()
angles = R.angles

file = open("RamachandranAnalysis.txt", 'w')
file.write("""Phi; Psi \n""")
file.close()

for i in range(len(R.angles)):
    for j in range(3):
        x = R.angles[i][j][0]
        y = R.angles[i][j][1]
        file = open(f"{dir}/ramachandran_analysis.txt", 'a')
        file.write(f"""{x} ; {y} \n """)
        file.close()

x_bin, y_bin = [], []
for timestep in angles:
    for coordinate in timestep:
        x_bin.append(coordinate[0])
        y_bin.append(coordinate[1])

plt.hist2d(x_bin, y_bin, bins=(120, 120), range=([-180, 180], [-180, 180]), density=True, cmap=plt.cm.jet)
plt.xticks([-180, -120, -60, 0, 60, 120, 180])
plt.yticks([-180, -120, -60, 0, 60, 120, 180])
plt.grid(visible=True, which="major")
plt.title(f'Ramachandran 2D Binning of {topol}')
plt.xlabel(r'$\Phi$')
plt.ylabel(r'$\Psi$')
plt.savefig(f'{dir}/RamachandranAnalysis.png')

