#!/bin/bash

step_cpus=4
max_cpus=$(<max_cpus.txt)
partition=$(<partition.txt)
account=$(<account.txt)
email_address=$(<email_address.txt)
email_info=$(<email_info.txt)

if (( $max_cpus < $step_cpus )); then
	cpus=$max_cpus
else
	cpus=$step_cpus
fi

nodes=$(( (cpus + 127) / 128 ))

sbatch << MINIMIZATION
#!/bin/bash
#SBATCH --job-name=6.1_minimization
#SBATCH -o step6.1_minimization.out
#SBATCH -e step6.1_minimization.err
#SBATCH -A ${account}
#SBATCH -t 02:00:00
#SBATCH -n ${cpus}
#SBATCH -N ${nodes}
#SBATCH -p ${partition}
#SBATCH --mail-user=${email_address}
#SBATCH --mail-type=${email_info}

module load gromacs
if [[ ! -f "step6.1_minimization.tpr" ]]; then
        gmx_mpi grompp -f step6.1_minimization.mdp -o step6.1_minimization.tpr -c step6.0_minimization.gro -r step5_charmm2gmx.pdb -p system.top -n index.ndx -maxwarn 1
fi

if [[ -f "step6.1_minimization.cpt" ]]; then
        mpirun -np ${cpus} gmx_mpi mdrun -deffnm step6.1_minimization -s step6.1_minimization.tpr -cpi step6.1_minimization.cpt
else
        mpirun -np ${cpus} gmx_mpi mdrun -deffnm step6.1_minimization
fi

if [ -f step6.1_minimization.gro ]; then
	./step6.2_minimization.sh
else
	echo "ERROR: step6.1_minimization failed, adjust and retry" >&2
fi
MINIMIZATION
