import os
import subprocess
from .env_vars import *


def format_pbs(sim_name, command, working_dir, queue, nodes, ppn, walltime, email, modules, env_variables):
    """ 
        everything is a string except:
            modules is a list of strings
            env_variables is a dict of KEY=value environment variables to set
    """
    # job is run in scratch_dir and output is copied to working_dir upon completion
    slurm = os.getenv('$SLURM_JOB_ID')
    print(f'slurm {slurm}')
    scratch_dir = '"/scratch/$SLURM_JOB_ID"' # leading slash for rsync
    modules = ' '.join(modules)
    outfile = os.path.join(working_dir, 'job.out')
    environment_variables = '\n'.join([f'export {key}={value}' for key, value in env_variables.items()])

    return (
        f'#!/bin/bash\n'
        f'#SBATCH --job-name={sim_name}\n'
        f'#SBATCH -p {queue}\n'
        f'#SBATCH --nodes={nodes}\n'
        f'#SBATCH --ntasks={ppn}\n'
        f'#SBATCH -t {walltime}\n'
        #f'#SBATCH -j oe\n' #TODO maybe -o
        f'#SBATCH -o "{outfile}"\n'
        # f'#SBATCH --mail-type=ALL\n'
        # f'#SBATCH --mail-user={email}\n'
        f'cd {scratch_dir}\n'
        f'module purge\n'
        f'module load {modules}\n'
        f'module list\n'
        f'{environment_variables}\n'
        f'time {command}\n'
        f'rsync -av {scratch_dir} "{working_dir}"\n'
    )

def submit_job(filename, contents):
    with open(filename, 'w') as f:
        f.write(contents)

    subprocess.run(['sbatch', filename])

class JobLauncher:
    """ a class that launches PBS jobs on the Razor machine (and possibly other machines) 
        arguments required in either __init__() or launch_job(): 
            sim_name, command, working_dir, queue, nodes, ppn, walltime, email
    """
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def launch_job(self, **kwargs):
        """ args must contain keyword arguments for: sim_name, command, working_dir, queue, nodes, ppn, walltime, email 
            these keyword arguments can either be passed into __init__() or launch_job()
            launch_job() kwargs override __init__() kwargs
        """
        args = {**self.__dict__, **kwargs}
        
        args['modules'] = [
            'mkl/18.0.2', 
            'vmd/1.9.3', 
            'python/3.7.3-anaconda', 
            'arbd/jan18', 
            'cuda/9.2',
            'cmake',
            'gcc/9.1.1',
            'openmpi/4.0.1'
        ]

        args['env_variables'] = {**cfg_dict}
        args['env_variables']['PATH'] = f'$PATH:{OXDNA_HOME}/build/bin'

        from pyoxdna.utils import mkdir
        mkdir(args['working_dir'])
        filename = os.path.join(args['working_dir'], 'job.slurm')
        contents = format_pbs(**args)
        submit_job(filename, contents)
