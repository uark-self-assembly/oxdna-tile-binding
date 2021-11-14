# this is a script that makes it easy to run python files
# usage: ./python3.sh [file.py] [arguments]
# NOTE: file.py must be a python file in ~/oxDNA-simulations/oxdna-tile-binding/
#       enclose multiple arguments in parentheses. i.e. python file.py "arg1 arg2 arg3"

module load python/3.7.3-anaconda
python3 ~/oxDNA-simulations/oxdna-tile-binding/$1 $2