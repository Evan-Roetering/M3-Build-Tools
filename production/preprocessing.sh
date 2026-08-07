#!/bin/bash

# ╔════════════════════════════════════════════════════════════════════════════════════════════════════════════════╗
# ║                          Preprocessing For Automated Martini Simulation with Gromacs                           ║
# ╟────────────────────────────────────────────────────────────────────────────────────────────────────────────────╢
# ║ Written By: Evan Roetering                                                                                     ║
# ║             Klauda Lab                                                                                         ║
# ║             University of Maryland                                                                             ║
# ║                                                                                                                ║
# ║ Last Edit:  5/20/2026                                                                                          ║
# ║                                                                                                                ║
# ╚════════════════════════════════════════════════════════════════════════════════════════════════════════════════╝

# ┌────────────────────────────────────────────────────────────────────────────────────────────────────────────────┐
# │                                       Import and Calculate Variables                                           │
# └────────────────────────────────────────────────────────────────────────────────────────────────────────────────┘

# import slurm info
partition=$(<partition.txt)
account=$(<account.txt)
email_address=$(<email_address.txt)
email_info=$(<email_info.txt)

# import production step variables
next_step=$(<next_step.txt)
last_step=$(<last_step.txt)
prev_step=$(( next_step - 1 ))

# import timestep variable
dt=$(<step_size.txt)                            # size of one timestep in ps
sim_steps=$(<total_timesteps.txt)               # full simulation time in timesteps
inc_steps=$(<increment_timesteps.txt)           # increment added for next step of simulation in timesteps
run_steps=$(<run_timesteps.txt)                 # length of simulation up to current point in timesteps
next_run_steps=$(( run_steps + inc_steps ))     # number of timesteps in next run
write_freq=$(<write_freq.txt)                   # distance between frames in timesteps

# convert timestep values to picoseconds
sim_ps=$(echo "$sim_steps * $dt" | bc -l)
add_ps=$(echo "$inc_steps * $dt" | bc -l)
prev_ps=$(echo "$run_steps * $dt" | bc -l)
next_ps=$(echo "$next_run_steps * $dt" | bc -l)

# convert picosecond values to nanoseconds
sim_ns=$(echo "$sim_ps * 1000" | bc -l)
add_ns=$(echo "$add_ps * 1000" | bc -l)
prev_ns=$(echo "$prev_ps * 1000" | bc -l)
next_ns=$(echo "$next_ps * 1000" | bc -l)

# convert nanosecond values to microseconds
sim_us=$(echo "$sim_ns * 1000" | bc -l)
add_us=$(echo "$add_ns * 1000" | bc -l)
prev_us=$(echo "$prev_ns * 1000" | bc -l)
next_us=$(echo "$next_ns * 1000" | bc -l)

# ┌────────────────────────────────────────────────────────────────────────────────────────────────────────────────┐
# │                                      Submit Preprocessing Job Via Slurm                                        │
# └────────────────────────────────────────────────────────────────────────────────────────────────────────────────┘
sbatch << PREPROCESSING
#!/bin/bash
#SBATCH --job-name=7.${next_step}_preprocessing
#SBATCH --open-mode=append
#SBATCH -o ./slurmlogs/7.${next_step}_production.out
#SBATCH -e ./slurmlogs/7.${next_step}_production.err
#SBATCH -A ${account}
#SBATCH -t 0:30:00
#SBATCH -p ${partition}
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mail-user=${email_address}
#SBATCH --mail-type=${email_info}

module load gromacs

if [[ ${next_step} -eq 0 ]]; then
	echo "Generating .tpr file for step ${next_step} of ${last_step} in the production run"
	gmx_mpi grompp -f step7.${next_step}_production.mdp -o step7.${next_step}_production.tpr -c step6.7_equilibration.gro -r step5_charmm2gmx.pdb -p system.top -n index.ndx -maxwarn 1
	if [[ -f step7.${next_step}_production.tpr ]]; then
		./mdrun.sh
	else
        	echo "step7.${next_step}_production.tpr missing" >&2
        	exit 1
	fi
else
	echo "Extending ${prev_ns} ns previous .tpr file for step ${prev_step} of ${last_step} in the production run by ${add_ns} ns"
	echo "New .tpr file will be generated for step ${next_step} of ${last_step} with ${next_ns} ns over ${run_steps} time steps"
	gmx_mpi convert-tpr -s step7.${prev_step}_production.tpr -o step7.${next_step}_production.tpr -extend ${add_ps}
	if [[ -f step7.${next_step}_production.tpr ]]; then
        	echo "$next_run_steps" > run_timesteps.txt
        	cp step7.${prev_step}_production.mdp step7.${next_step}_production.mdp
        	sed -i "s/^nsteps[[:space:]]*=.*/nsteps                   = ${next_run_steps}/" step7.${next_step}_production.mdp
        	rm step7.${prev_step}_production*.trr
        	./mdrun.sh
	else
        	echo "ERROR: step7.${next_step}_production.tpr missing" >&2
        	exit 1
	fi
fi
PREPROCESSING
