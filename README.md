#  oxdna-tile-binding

##  Setting Up oxdna-tile-binding

1. Clone the repository
2. Run `cd oxdna-tile-binding`
3. Run `module load python/3.6.0-anaconda`
4. Run `echo ". /share/apps/python/anaconda3-3.6.0/etc/profile.d/conda.sh" >> ~/.bashrc`
5. Run `echo "conda activate" >> ~/.bashrc`
6. Set the paths in oxdna-tile-binding/config.yml for your system
7. Run `conda activate pyoxdna` (make sure your conda environment is set to pyoxdna every time you use these scripts. Sometimes, if you log off of a server, it may reset your conda environment, and you will have to run this command again)

If you run into any issues running things locally (with errors like `KeyError: 'PYOXDNA_HOME'`),  you may need to set the environment variables by following these instructions:
-   If you want to temporarily set your environment variables (until you exit) run `export ENV_NAME=PATH` (e.x. export `OXDNA_HOME=/home/username/oxDNA`) on each variable
-   If you want them to always be set you can run  
    `echo "export OXDNA_HOME=/home/username/oxDNA" >> ~/.bashrc `
    On each variable OR open `~/.bashrc` with an editor and add these lines at the bottom:  
    ```
    export HOME=/home/username  
    export SIM_HOME=/home/username/simulations  
    export OXDNA_HOME=/home/username/oxDNA  
    export PYOXDNA_HOME=/home/username/oxdna-tile-binding/pyoxdna  
    ```
    (make sure you change the paths to the correct ones for your system)

## Setting up oxDNA
[oxDNA](https://dna.physics.ox.ac.uk/index.php/Main_Page) is a DNA molecular simulator and analysis tool developed at the University of Oxford. Good documentation can be found [here](https://dna.physics.ox.ac.uk/index.php/Documentation). 

### Setup
*Notes before you start*: 
 - You may need to change the version of gcc
 - in select_compute_arch.cmake replace 4 instances of “VERSION_GREATER_EQUAL” with “VERSION_GREATER"

Follow [these](https://dna.physics.ox.ac.uk/index.php/Download_and_Installation) directions. We used these following commands to compile:
```
    cd oxDNA
    mkdir build
    cd build
    module load cuda/9.2 cmake gcc/9.1.1
    cmake -DCUDA=1 ..
    make
```
### Using oxDNA
#### oxDNA File Formats
oxDNA has 4 main files:
 -   An input parameter file (usually named “input” or “input.params” or “input.dat”). This specifies all of the parameters for the simulation, including the type of simulation, the number of simulation steps, the temperature of the simulation, how often to print output, etc.
 -   An input topology file (top file). This specifies the bases and topology of the DNA strands to be simulated. This tells the program which nucleotides are present in the simulation and how they are connected. Documentation [here](https://dna.physics.ox.ac.uk/index.php/Documentation#Topology_file).
-   An input configuration file (conf file). This specifies how the nucleotides defined in the topology file are oriented in space. Documentation [here](https://dna.physics.ox.ac.uk/index.php/Documentation#Configuration_file).
-   An output trajectory file. This is a giant file showing the position and orientation of all the DNA strands at each recorded timestep during the simulation. It’s a plain-text file containing multiple configuration files separated by newlines. You can change how often the configuration is recorded by modifying the print_conf_interval input parameter.

Because each of these file formats is plain text and contains only decimal numbers, it is possible to create and manipulate these files without interacting with external programs, such as scadnano. For example, it’s possible to generate custom (or random) configuration and topology files or to splice configurations together. Examples of this can be found in create_tiles.py and computation_experiment.py. Although it’s a huge file, pyoxdna/analysis/base.py, the python code that comes with oxDNA, is a good reference for file formats and file modification.

A great visualization tool for oxDNA simulations is [Oxdna-viewer](https://sulcgroup.github.io/oxdna-viewer/). To load a simulation, click “open” and then select either a top and conf file OR a top and trajectory file.

## Scripts and Packages in oxdna-tile-binding

### utils/job_launcher.py
This is a python module for easily launching jobs on Razor. Given a command, working_directory, dependencies, job file, etc., this module will automatically create a job script and submit a job to your computing cluster on your behalf. Documentation can be found in job_launcher.py and example usage can be found on line 120 in run_simulations.py.

Many scripts allow the user to input a job_file.

The job input file should have the command the script should be run with in the first line (for example `sbatch` for Slurm or `qsub` for PBS/Torque). In the following lines, users should put the job script they wish to run the program with. This job script should only have the system specific parameters including `[job_name]` and `[out_file]` as flags. All other job specific parameters (such as loading modules or running commands) will be formatted by job_launcher.py, and `[job_name]` and `[out_file]` will be replaced with the job appropriate data.
Example of an input job script:
```
sbatch  
#!/bin/bash  
#SBATCH --job-name=[job_name]  
#SBATCH -p gpu72  
#SBATCH --nodes=1  
#SBATCH --ntasks=6  
#SBATCH -t 72:00:00  
#SBATCH -o [out_file]  
cd "/scratch/$SLURM_JOB_ID"
```
The job.txt that results when used as input to job_launcher.py (which will be run with `sbatch`):
```
#!/bin/bash  
#SBATCH --job-name=example1  
#SBATCH -p gpu72  
#SBATCH --nodes=1  
#SBATCH --ntasks=6  
#SBATCH -t 72:00:00  
#SBATCH -o /home/user/example1/example.out  
cd "/scratch/$SLURM_JOB_ID"  
module purge  
module load mkl/18.0.2  
module list  
export HOME=/home/user  
export PATH=$PATH:/home/user/oxDNA/build/bin  
time python3 /scrfs/storage/user/home/example.py
```

### Pyoxdna
A python package written by Kyle Sadler to interact with oxDNA in python. It has two important modules:
-   **pyoxdna** is the main module for running simulations in oxDNA. Given a directory and a python dictionary containing [oxDNA input options](https://dna.physics.ox.ac.uk/index.php/Documentation#Input_file), pyoxdna will run an oxDNA simulation. Documentation can be found in pyoxdna/pyoxdna.py and example usage can be found on line 159 in run_simulations.py
-   **dna_relaxer** is a module for automatically relaxing DNA configurations. DNA configurations must be relaxed in order to simulate them in oxDNA. [Here](https://docs.google.com/document/d/1zP__47jWaXR0NSNC0wEH4XGCFSHxNVxmBai3VXTGlN0/edit) are some notes on relaxation. Documentation for DNARelaxer can be found in pyoxdna/dna_relaxer.py and example usage can be found on line 18 in run_simulations.py

### analyze.py
This is a script for analyzing tile binding in oxDNA after simulation. Given input, topology, and trajectory files, analyze.trim_strands() creates a topology and trajectory files with ONLY the strands that bind during the simulation. This makes it easier to see strand interaction in a large simulation.

*Run with*: `python analyze.py [input_file] [trajectory_file] [top_file] [job_file]`. If you would like to run locally, run without a job_file.

Since we can compute simulation binding events (e.g. nucleotide1 binds to nucleotide2 at time t), this script can be improved to filter for bonds with x or more bound nucleotides. Therefore, this can be used to automatically count the number of full tile bindings in a simulation.
Credit: this script is based on the work of Michael Sharp who worked on a similar project with Dr. Patitz a few years ago.

### computation_experiment.py 
A script to profile the run time of oxDNA’s simulations using both the GPU and CPU based on the number of nucleotides in the simulations. Results from the experiment are found in [HOME]/results.txt

*Run with*: `python computation_experiment.py run [num_strands] [GPU:True/False][output_dir (optional, default is SIM_HOME])]` to run the experiment with a specific number of strands or `python computation_experiment.py [start_num_strands] [end_num_strands] [GPU:True/False] [job_file]` to run the experiment with values in the range of  start_num_strands to end_num_strands.

We recommend that the user run a range of sizes multiple times to determine which sizes of structures are best suited for GPU or CPU use.

### create_tiles.py
A script to generate a grid of aligned tiles for simulation. This works by “copying and pasting” a conf and top file of two complementary tiles in order to make an x by y grid. Input files for this program can be found in oxdna_files/tiles/original. Output files can be found in oxdna_files/tiles.

*Run with*: `python create_tiles.py run [file_name] [num_tiles] [output_dir]` to run

### run_analysis_multi.py
A small script to run other scripts for multiple iterations.
Uses run_analysis.txt as input

*Contents of run_analysis.txt*: 

&emsp;First line: number of iterations to run

&emsp;Second line: the string to run

*Run with* `python run_analysis_multi.py`

### run_simulation.py
This program runs a simulation given a config file, topology file. There are several options the user can choose for relaxation and simulation. 

*Run with*: `python run_simulation.py -t [top_file] -c [conf_file] -j [job_file] -r` (this command will just preform relaxation)
Only a topology file, a config file, and either the -r or -m flag  are required for run_simulation.py . To run locally, simply exclude the -j flag and the job_file. To see more configuration options, run `python run_simulation.py -h`.

Depending on your system, the output files may end up in a temporary folder if you run this script as a job. Be sure to move those to a more permanent location before they expire.

### tile_binding_auto_monitoring_light.py
This script runs and analyzes tile binding simulations to see if tiles will attach or detach.

*Run with*: `python tile_binding_auto_monitoring_light.py -s [sim_conf_file] -b [bonds_file] -i [num_iterations] -o [out_dir]`. Only the simulation configuration file, bonds file, output directory, and number of iterations to run are required. To see more configuration options, run `python tile_binding_auto_monitoring_light.py -h`

An example of a simulation configuration can be found as a result of running run_simulation.py or in `pyoxdna/configs/molecular-dynamic.conf` (though this file is missing a topology and conf_file entry).

The `bonds_file` must have the .py file extension, and should look like the example as follows:

```
        STARTING_BONDS = [
            [944, 498],
            [945, 497],
            [946, 496],
            [947, 495],
            [948, 494]
        ]
        
        TARGET_BONDS = [
            [949, 576],
            [950, 575],
            [951, 574],
            [952, 573],
            [953, 572],
        
            [954, 456],
            [955, 455],
            [956, 454],
            [957, 453],
            [958, 452],
        
            [959, 639],
            [960, 638],
            [961, 637],
            [962, 636],
            [963, 635],
        ]
```
Where each  `[ID1,ID2]` is a pair of IDs of a nucleotide in the topology file. The program exits when **all** of the `STARTING_BONDS` have dissolved (printing `THE TILE DETACHED!`) or **any** of the `TARGET_BONDS` have formed (printing `TILE IS ATTACHED!`).
