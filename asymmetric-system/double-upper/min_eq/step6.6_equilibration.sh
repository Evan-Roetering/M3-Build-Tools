#!/bin/bash

step_cpus=128
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
#SBATCH --job-name=6.6_equilibration
#SBATCH -o step6.6_equilibration.out
#SBATCH -e step6.6_equilibration.err
#SBATCH -A ${account}
#SBATCH -t 12:00:00
#SBATCH -n ${cpus}
#SBATCH -N ${nodes}
#SBATCH -p ${partition}
#SBATCH --mail-user=${email_address}
#SBATCH --mail-type=${email_info}

module load gromacs
if [[ ! -f "step6.6_equilibration.tpr" ]]; then
        gmx_mpi grompp -f step6.6_equilibration.mdp -o step6.6_equilibration.tpr -c step6.5_equilibration.gro -r step5_charmm2gmx.pdb -p system.top -n index.ndx -maxwarn 1
fi

if [[ -f "step6.6_equilibration.cpt" ]]; then
        mpirun -np 16 gmx_mpi mdrun -deffnm step6.6_equilibration -s step6.6_equilibration.tpr -cpi step6.6_equilibration.cpt
else
        mpirun -np 16 gmx_mpi mdrun -deffnm step6.6_equilibration
fi

if [ -f step6.6_equilibration.gro ]; then
        ./step6.7_equilibration.sh
else
        echo "ERROR: step6.6_equilibration failed, adjust and retry" >&2
fi
EQUILIBRATION
