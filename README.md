# oxdna-tile-binding

## Setting Up
1. Clone the repository
2. Run `cd oxdna-tile-binding`
3. Run `module load python/3.6.0-anaconda`
4. Run `echo ". /share/apps/python/anaconda3-3.6.0/etc/profile.d/conda.sh" >> ~/.bashrc`
5. Run `echo "conda activate" >> ~/.bashrc`
6. Set the paths in oxdna-tile-binding/config.yml
7. You may need to set some paths in your environment as well
    - If you want to temporarily set your environment variables (until you exit) run `export ENV_NAME=PATH` (e.x. `export OXDNA_HOME=/home/username/oxDNA`) on each variable
    - If you want them to always be set you can run 
        - `echo "export OXDNA_HOME=/home/username/oxDNA" >> ~/.bashrc` 
    on each variable 
    - OR open ~/.bashrc with an editor and add these lines at the bottom:
        ```
        export HOME=/home/username
        export SIM_HOME=/home/username/simulations
        export OXDNA_HOME=/home/username/oxDNA
        export PYOXDNA_HOME=/home/username/oxdna-tile-binding/pyoxdna
        ```
        (make sure you change the paths to the correct ones for your system)
8. Run `conda activate pyoxdna`

