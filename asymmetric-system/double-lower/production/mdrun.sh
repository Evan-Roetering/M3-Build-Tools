#!/bin/bash

# import variables
cores=$(<cores.txt)
nodes=$(( (cores + 127) / 128 ))
partition=$(<partition.txt)
account=$(<account.txt)
email_address=$(<email_address.txt)
email_info=$(<email_info.txt)
run_steps=$(<run_timesteps.txt)
this_step=$(<next_step.txt)
sim_steps=$(<total_timesteps.txt)
inc_steps=$(<increment_timesteps.txt)
dt=$(<step_size.txt)
write_freq=$(<write_freq.txt)
last_step=$(<last_step.txt)
prev_step=$(( this_step - 1 ))
next_step=$(( this_step + 1 ))

sbatch << MDRUN
#!/bin/bash
#SBATCH --job-name=7.${this_step}_production
#SBATCH --open-mode=append
#SBATCH -o ./slurmlogs/7.${this_step}_production.out
#SBATCH -e ./slurmlogs/7.${this_step}_production.err
#SBATCH -A ${account}
#SBATCH -t 5-00:00:00
#SBATCH -n ${cores}
#SBATCH -N ${nodes}
#SBATCH -p ${partition}
#SBATCH --mail-user=${email_address}
#SBATCH --mail-type=${email_info}

echo "Running simulation for step ${this_step} of ${last_step} in the production run"

module load gromacs
if [[ "$this_step" == "0" ]]; then
	if [[ -f step7.${this_step}_production.cpt ]]; then
		mpirun -np ${cores} gmx_mpi mdrun -deffnm step7.${this_step}_production -s step7.${this_step}_production.tpr -cpi step7.${this_step}_production.cpt
	else
		mpirun -np ${cores} gmx_mpi mdrun -deffnm step7.${this_step}_production
	fi
else
	if [[ -f step7.${this_step}_production.cpt ]]; then
		mpirun -np ${cores} gmx_mpi mdrun -deffnm step7.${this_step}_production -s step7.${this_step}_production.tpr -cpi step7.${this_step}_production.cpt -noappend
	else
		mpirun -np ${cores} gmx_mpi mdrun -deffnm step7.${this_step}_production -s step7.${this_step}_production.tpr -cpi step7.${prev_step}_production.cpt -noappend
	fi
fi

shopt -s nullglob
grofile=(step7.${this_step}_production.*gro)
grocount=\${#grofile[@]}

echo \$grofile
echo \$grocount

if (( grocount == 0 )); then
        echo "Error: no .gro file found. Production run may have failed" >&2
        exit 1
elif (( grocount == 1 )); then
        echo "Production run ran properly, one .gro file found: \${grofile[0]}"
        if [[ "$this_step" == "$last_step" ]]; then
                echo "Completed all production steps - moving to combine .xtc files"
                ./cleanup.sh
        else
                echo "$next_step" > next_step.txt
                ./preprocessing.sh
        fi
else
        echo "Error: multiple .gro files were found (\$grocount):" >&2
        printf '  %s\\n' "\${grofile[@]}"
        exit 1
fi
MDRUN
