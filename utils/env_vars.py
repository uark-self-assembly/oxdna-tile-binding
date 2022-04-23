import os
import yaml

""" 
creates a PBS job and runs it on the local machine
    (This is only needed on a mainframe computer, such as Razor)

NOTE: all commands / programs run in a job must use the current working directory as the output directory
specify the output_dir to change where the output is stored 
This is because job_launcher.py uses different temporary and output directories in the backend

python /home/user/oxDNA-simulations/kyle/mycode/launch_job.py

    /home/user (HOME)
        | --> oxdna-code/oxDNA (OXDNA_HOME)
        | --> results.csv
        | --> simulations
            | --> simulation_name (working_dir)
                | --> oxdna outputs and logs
                | --> sim.conf
                | --> top.conf
                | --> job.pbs (job PBS script)
                | --> job.out (simulation output and errors)
        | --> oxDNA-simulations/kyle/mycode
            | --> launch_job.py
            | --> relax_and_sim.py
            | --> pyoxdna (PYOXDNA_HOME)

"""
environment_config = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'config.yml'))
with open(environment_config, 'r') as f:
    cfg_dict = yaml.safe_load(f)

HOME = cfg_dict['HOME']
SIM_HOME = cfg_dict['SIM_HOME']
OXDNA_HOME = cfg_dict['OXDNA_HOME']
PYOXDNA_HOME = cfg_dict['PYOXDNA_HOME']
EMAIL_ADDRESS = cfg_dict['EMAIL_ADDRESS']
DEBUG_LEVEL = cfg_dict['DEBUG_LEVEL']

print(f"Setting environment variables from file '{environment_config}':")
print('\n'.join(f'{k}: {v}' for k,v in cfg_dict.items()))
print()
