import os
import subprocess
from .env_vars import *


def submit_job(filename, sim_name, command, working_dir, modules, env_variables, in_file = 'example.job'):
    """ 
        everything is a string except:
            modules is a list of strings
            env_variables is a dict of KEY=value environment variables to set
    """
    # job is run in scratch_dir and output is copied to working_dir upon completion
    modules = ' '.join(modules)
    outfile = os.path.join(working_dir, 'job.out')
    environment_variables = '\n'.join([f'export {key}={value}' for key, value in env_variables.items()])

    # load the following text into the contents string
    contents =\
    f'module purge\n' +\
    f'module load {modules}\n' +\
    f'module list\n' +\
    f'{environment_variables}\n' +\
    f'time {command}\n' +\
    f'rsync -av . "{working_dir}"\n'
    
    with open(in_file, 'r') as f:
        run_command = f.readline().rstrip() # Get the command to run the file with from in_file
        script = f.read()
        # Find and replace [job_name] and [out_file]
        script = script.replace('job_name', sim_name) 
        script = script.replace('out_file', outfile)
        
        with open(filename, 'w') as g:
            g.write(script)
            # write the contents of in_file to filename
            g.write(contents)
        print(run_command)

    subprocess.run([run_command, filename])

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
        filename = os.path.join(args['working_dir'], 'job.txt')
        submit_job(filename, **args)

