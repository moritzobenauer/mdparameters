# Changing the temperature

NewTemperature(){

temp=$1
name="T$temp"
mkdir $name
cp -r "files/"/* "$name/"
cd $name
file_path="nvt.mdp"
sed -i "s/300/$temp/g" "$file_path"
file_path="npt.mdp"
sed -i "s/300/$temp/g" "$file_path"
file_path="md.mdp"
sed -i "s/300/$temp/g" "$file_path"
cd ..
}

RunMD(){
temp=$1
name="T$temp"
cd $name
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt -v
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -deffnm npt -v
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr
gmx mdrun -deffnm md_0_1 -v
gmx trjconv -s md_0_1.tpr -f md_0_1.xtc -o md_0_1_noPBC.xtc -pbc mol -center
gmx trjconv -s md_0_1.tpr -f md_0_1.gro -o md_0_1_noPBC.gro -pbc mol -center
cd ..	
}

source /usr/local/gromacs/bin/GMXRC
cd "files"

gmx pdb2gmx -ff oplsaa-mod2023  -ter  -o raw.gro -water tip4 -v -f FHFHF-PEG10-FHFHF.pdb
gmx editconf -f raw.gro -o raw_newbox.gro -c -d 1.0 -bt dodecahedron
gmx grompp -f preminim.mdp -c raw_newbox.gro -p topol.top -o vac_em.tpr -maxwarn 2
gmx mdrun -deffnm vac_em -v
gmx solvate -cp vac_em.gro -cs tip4p.gro -o raw_solv.gro -p topol.top
gmx grompp -f ions.mdp -c raw_solv.gro -p topol.top -o ions.tpr
gmx genion -s ions.tpr -o raw_solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.01
gmx grompp -f minim.mdp -c raw_solv_ions.gro -p topol.top -o em.tpr -maxwarn 2
gmx mdrun -deffnm em -v

cd ..

NewTemperature 300
NewTemperature 305
NewTemperature 310
NewTemperature 315
NewTemperature 320

RunMD 300
RunMD 305
RunMD 310
RunMD 315
RunMD 320

echo "Erfolgreich abgeschlossen."







