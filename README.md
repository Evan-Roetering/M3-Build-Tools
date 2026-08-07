# M3-Build-Tools
### Guide for Building Membrane Simulations with Martini 3 (Work in Progress)
&emsp;&#8202;│  
&emsp;├──🞂 Evan Roetering  
&emsp;├──🞂 Klauda Lab  
&emsp;├──🞂 University of Maryland  
&emsp;└──🞂 Last Edit: 8/7/2026  

## Step by Step Asymmetric Bilayer Build guide:
1. Start with the asymmetric-system template on zaratan
2. Gather upper leaflet and lower leaflet composition in terms of concentration
3. Build and Simulate Symmetric systems in the `asymmetric-system/double-upper` and `asymmetric-system/double-lower` directories  
&emsp; a. navigate to `asymmetric-system/double-\<upper/lower\>/system_build`  
&emsp; b. format make_toppar.sh with the lipids needed for your system
&emsp; c. execute make_toppar.sh  
&emsp; d. assemble symmetric bilayer with my modified version of insane using the syntax:  
&emsp;&emsp; `python3 evan_insane.py -a A -x X -y Y -z Z -o membrane.gro -p insane.top -u NAME1:num1 -u NAME2:num2 -u NAME3:num3 ... -l NAME1:num4 -l NAME2:num5 -l NAME3:num6 ... -sol W`  
&emsp;&emsp;&emsp; A = Estimate of the Average Area per lipid in nm, I use 0.6 typically, but if unsure check CHARMM-GUI for literature values and calculate based on those and the composition  
&emsp;&emsp;&emsp; X = X box dimension in nm (should equal Y)  
&emsp;&emsp;&emsp; Y = Y box dimension in nm (should equal X)  
&emsp;&emsp;&emsp; Z = Z box dimension in nm (usually 10 is a good guess)  
&emsp;&emsp;&emsp; NAME# = 4 letter code for lipid used in make_toppar.sh  
&emsp;&emsp;&emsp; num# = number of the corresponding lipid  
&emsp;&emsp;&emsp; -u = flag indicating lipid is being entered for upper leaflet (should be identical to -l flags for this step)  
&emsp;&emsp;&emsp; -l = flag indicating lipid is being entered for the lower leaflet (should be identical to -u flags for this step)  
&emsp;&emsp;&emsp; -sol W = use water as solvent  
&emsp;&emsp;&emsp; NOTE 1: Insane places in a square grid, so the sum of lipid counts in each leaflet should be a perfect square  
&emsp;&emsp;&emsp; NOTE 2: Insane places a number of lipids based on the ratio of input numbers, to ensure proper lipid counts X and Y should equal sqrt(A*sum(leaflet counts))  
&emsp; e. name system and make `system.top` file using `./make_system_top.sh "\<system name\>"`  
&emsp; f. load gromacs with `module load gromacs`  
&emsp; g. run `gmx_mpi grompp -f ions.mdp -c membrane.gro -p system.top -o ions.tpr -maxwarn 1`  
&emsp; h.  run `gmx_mpi genion -s ions.tpr -o system_ions.gro -p system.top -neutral -pname NA`  
&emsp;&emsp; select waters as the solvent  
&emsp; i. make MEMBRANE and SOLUTE index with `gmx_mpi make_ndx -f system_ions.gro -o index.ndx`  
&emsp;&emsp; i. use `2|3|4|...` to select all lipid species in an index with a new number  
&emsp;&emsp; ii. use `name # MEMBRANE`, where # is the index number for the new group to give it the name MEMBRANE  
&emsp;&emsp; iii. use `!#` to select all atoms not in the MEMBRANE group and save them under a new index number  
&emsp;&emsp; iv. use `name #\* SOLUTE`, where #\* is the new index number, to name the group SOLUTE  
&emsp;&emsp; v. use `q` to save the index file as index.ndx  
&emsp; j. run `gmx_mpi editconf -f system_ions.gro -o step5_charmm2gmx.pdb` to save in the same output format as CHARMMM-GUI  
4. navigate up to `asymmetric-system/double-\<upper/lower\>` and edit `min_eq.start` with your email address and a number of cpus (I recommend 8-16 for small systems, 32-64 for medium to large system, and 128 for very large systems)
5. run it with `./min_eq.start`
6. after minimization and equilibration edit `production.start` with email and number of cores (64 is a good starting point for medium systems, but may be unstable for small systems or slow for large ones)
7. run it with `./production.start`
8. Visualize double leaflet simulation in vmd using files generated in the `asymmetric-system/double-\<upper/lower\>/martiniglass` directory
9. If undulations appear to be present in visualizations use `./undulation_area.sh` in the `asymmetric-system/double-\<upper/lower\>/analysis` directory, otherwise use `flat_area.sh` edit the email address in the script before running
10. Check outputs for area metrics
11. Based on area metrics, quantify the number of lipids per leaflet to cover the same area
12. navigate to the `asymmetric-system/asymmetric/system` build directory
13. construct lipid topology with `./make_toppar.sh`
14. construct asymmetrical bilayer using insane with the same syntax as earlier:  
&emsp; `python3 evan_insane.py -a A -x X -y Y -z Z -o membrane.gro -p insane.top -u NAME1:num1 -u NAME2:num2 -u NAME3:num3 ... -l NAME1:num4 -l NAME2:num5 -l NAME3:num6 ... -sol W`
15. name system and make `system.top` file using `./make_system_top.sh "\<system name\>"`
16. run `gmx_mpi grompp -f ions.mdp -c membrane.gro -p system.top -o ions.tpr -maxwarn 1`
17. run `gmx_mpi genion -s ions.tpr -o system_ions.gro -p system.top -neutral -pname NA`  
&emsp; select waters as the solvent  
18. make MEMBRANE and SOLUTE index with `gmx_mpi make_ndx -f system_ions.gro -o index.ndx`  
&emsp; a. use `2|3|4|...` to select all lipid species in an index with a new number  
&emsp; b. use `name # MEMBRANE`, where # is the index number for the new group to give it the name MEMBRANE  
&emsp; c. use `!#` to select all atoms not in the MEMBRANE group and save them under a new index number  
&emsp; d. use `name #\* SOLUTE`, where #\* is the new index number, to name the group SOLUTE  
&emsp; e. use `q` to save the index file as index.ndx  
19. run `gmx_mpi editconf -f system_ions.gro -o step5_charmm2gmx.pdb` to save in the same output format as CHARMMM-GUI
20. navigate back one directory to `asymmetric-system/asymmetric` and edit `min_eq.start` with email and cpus
21. run `./min_eq.start`
22. after minimization and equilibration, edit `production.start` with email, cores, sim_steps, and n_runs
&emsp; email - email notifications should be sent to
&emsp; cores - number of cpu cores to use for run
&emsp; n_runs - number of pieces to separate simulation into (I recommend keeping them to a minimum of 100 ns per run and no more than 1 microsecond, make sure there are at least 2 or the script wont work properly)
&emsp; sim_steps - number of 20 fs steps to conduct simulation over (default is 25 microseconds aka 1250000000 steps, there are 50000000 steps per microsecond)
23. run with `./production.start`
