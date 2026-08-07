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

sbatch << EQUILIBRATION
#!/bin/bash
#SBATCH --job-name=6.3_equilibration
#SBATCH -o step6.3_equilibration.out
#SBATCH -e step6.3_equilibration.err
#SBATCH -A ${account}
#SBATCH -t 12:00:00
#SBATCH -n ${cpus}
#SBATCH -N ${nodes}
#SBATCH -p ${partition}
#SBATCH --mail-user=${email_address}
#SBATCH --mail-type=${email_info}

module load gromacs
if [[ ! -f "step6.3_equilibration.tpr" ]]; then
        gmx_mpi grompp -f step6.3_equilibration.mdp -o step6.3_equilibration.tpr -c step6.2_minimization.gro -r step5_charmm2gmx.pdb -p system.top -n index.ndx -maxwarn 1
fi

if [[ -f "step6.3_equilibration.cpt" ]]; then
        mpirun -np ${cpus} gmx_mpi mdrun -deffnm step6.3_equilibration -s step6.3_equilibration.tpr -cpi step6.3_equilibration.cpt
else
        mpirun -np ${cpus} gmx_mpi mdrun -deffnm step6.3_equilibration
fi

if [ -f step6.3_equilibration.gro ]; then
        ./step6.4_equilibration.sh
else
        echo "ERROR: step6.3_equilibration failed, adjust and retry" >&2
fi
EQUILIBRATION
