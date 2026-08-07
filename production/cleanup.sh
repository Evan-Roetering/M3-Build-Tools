#!/bin/bash

# ============================ Cleanup After Completing Simulation ============================
# import variables
cores=$(<cores.txt)
nodes=$(( (cores + 127) / 128 ))
partition=$(<partition.txt)
account=$(<account.txt)
email_address=$(<email_address.txt)
email_info=$(<email_info.txt)
last_step=$(<last_step.txt)

sbatch << CLEANUP
#!/bin/bash
#SBATCH --job-name=step7_cleanup
#SBATCH --open-mode=append
#SBATCH -o ./slurmlogs/cleanup.out
#SBATCH -e ./slurmlogs/cleanup.err
#SBATCH -A ${account}
#SBATCH -p ${partition}
#SBATCH -t 2-00:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mail-user=${email_address}
#SBATCH --mail-type=${email_info}

module load gromacs

# ================================== Combine .xtc Files =======================================
echo "Concatenating trajectories in the following order:"
ls step7*production.*xtc | sort -V
gmx_mpi trjcat -f \$(ls step7*production.*xtc | sort -V) -o step7_production.xtc

# ============================ Trim Frames to Preserve Memory =================================
if [[ -f step7_production.xtc ]]; then
	gmx_mpi trjconv -f step7_production.xtc -o step7_production_trimmed_10.xtc -skip 10
	gmx_mpi trjconv -f step7_production.xtc -o step7_production_trimmed_25.xtc -skip 25
else
        echo "ERROR: trajectory files were not combined" >&2
        exit 1
fi

# ============================ Move to Destination Directories ================================
if [[ -f step7_production_trimmed_10.xtc ]]; then
	mkdir -p ../analysis
	mkdir -p ../analysis/gro_snapshots
	mkdir -p ../martini_glass
	cp -r ./toppar ../analysis/toppar
	cp -r ./toppar ../martini_glass/toppar
	cp ./system.top ../analysis/system.top 
        cp ./system.top ../martini_glass/system.top
        cp ./step7_production_trimmed_25.xtc ../analysis/trajectory.xtc
        cp ./step7_production_trimmed_10.xtc ../martini_glass/trajectory.xtc
        cp ./step6.7_equilibration.gro ../analysis/topology.gro
        cp ./step6.7_equilibration.gro ../martini_glass/topology.gro
	cp ./step7.${last_step}_production.tpr ../analysis/input.tpr
	cd ../analysis
	gmx_mpi trjconv -s "input.tpr" -f "trajectory.xtc" -o frame_.gro -sep <<TRJCONV
	0
	TRJCONV
	mv frame_*.gro ./gro_snapshots/
else
	echo "ERROR: trajectory files were not trimmed" >&2
	exit 1
fi
CLEANUP
