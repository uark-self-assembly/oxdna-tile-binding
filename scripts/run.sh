# this script makes it easy to run a simulation.
# You can move this to your home directory or add it to your path to make launching a simulation even easier.

# the directory where top and conf files live
OXNDA_FILES="$HOME/oxDNA-simulations/oxdna-tile-binding/oxdna_files"

# change this to the files you want to run the simulation with
PREFIX="$OXNDA_FILES/tiles/10_500"

SIMULATIONS=16 # how many simulations to launch
SIMULATION_STEPS=30000 # how many steps per simulation

# cd; ./python3.sh run_simulation.py "$PREFIX.conf $PREFIX.top $SIMULATIONS $SIMULATION_STEPS"
echo "$PREFIX.conf $PREFIX.top $SIMULATIONS $SIMULATION_STEPS"