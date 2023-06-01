source /usr/local/gromacs/bin/GMXRC


gmx trjconv -s md_0_1.tpr -f md_0_1.xtc -o md_0_1_noPBC.xtc -pbc mol -center

gmx trjconv -s md_0_1.tpr -f md_0_1.gro -o md_0_1_noPBC.gro -pbc mol -center

cp md_0_1_noPBC.gro analysis/md_0_1_noPBC.gro
cp md_0_1_noPBC.xtc analysis/md_0_1_noPBC.xtc


cd analysis
python3 analysis.py

