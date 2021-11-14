from .utils import *

class pyoxdna:
    """ this class runs and manages oxDNA simulations. PYOXDNA_HOME must be set to run correctly """
    
    def __init__(self, output_dir, output_level=3):
        self.output_dir = output_dir or 'pyoxdna_' + current_time()
        mkdir(self.output_dir)
        
        """ 
        level of output to produce: 
            3+ keeps everything
            2 removes input files
            1 removes input and output conf files
            0 removes everything (including output dir)
        """
        self.output_level = output_level
        
        self.paths_generator = self._yield_paths()

        self.pid_file = os.path.join(self.output_dir, '.pids')

        print_message(5, f"In pyoxdna init, output_level = {output_level}, pid_file = '{self.pid_file}'")

    def relax(self, step_1, step_2):
        self.run(get_relax_configs(step_1, step_2))

    def simulate(self, configs):
        self.run(get_simulation_configs(configs))

    def run(self, configs):
        """ Runs a simulation based on configs
            Confgs may be an oxDNA configuration input dict or a list of configuration dicts to run synchronously
            If configs is a list of configuration dicts, simulations are chained together by feeding the output .conf
                of the first config dictionary into the input of the next

            throws error if oxDNA errors
        """

        if isinstance(configs, list):
            for config in configs:
                self.run_one(config)
        elif isinstance(configs, dict):
            self.run_one(configs)
        else:
            raise ValueError("configs should be a config dict or a list of config dicts")
        
            
#===============================================================================
#         if isinstance(configs, list):
#             if len(configs) == 1:
#                 self.run_one(configs[0])
#             else:
#                 for i, config in enumerate(configs):
#                     if i > 0:
#                         config['conf_file'] = output_path # chain simulations together
#                         input_path, output_path = next(self.paths_generator)
# 
#                     self.run_one(config)
# 
#         else: # only one config dictionary
#             self.run_one(configs)
#===============================================================================

    def run_one(self, config):
        """ writes config to file and runs oxDNA simulation 
            throws error if oxDNA errors
            THIS FUNCTION BLOCKS (runs synchonously, stopping code flow until it is finished)
        """
        
        #input_path, output_path = next(self.paths_generator)
        input_path = os.path.join(self.output_dir, "input.config")
        output_path = os.path.join(self.output_dir, config['lastconf_file'])
        
        for name in ['energy', 'trajectory', 'log']:
            config[name+'_file'] = os.path.join(self.output_dir, name + '.dat')

        self.config = config # save so outside scripts can access variables
        # print(output_path +' oxDNA config:')
        # pprint(config)

        write_config(config, input_path)
        cmd = ['bash', os.path.join(os.environ['PYOXDNA_HOME'], 'simulate.sh'), input_path]
        print_message(5, f"In run_one, about to execute command: '{cmd}'")
        self.process = Popen(cmd)

        with open(self.pid_file, 'a+') as f:
            f.write(str(self.process.pid)+'\n')

        return_code = self.process.wait()

        if self.output_level <= 2:
            rm(input_path)

            if self.output_level <= 1:
                rm(output_path)

        if return_code != 0: # did not exit sucessfully
            self.exit()
            raise CalledProcessError(return_code, cmd)
    
    def _yield_paths(self):
        """ generator yielding input and output file paths """
        i = 0
        while True:
            # yield input_path, output_path
            yield [ os.path.join(self.output_dir, x) for x in [f'{i}.params', f'{i}.conf'] ]
            i += 1

    def exit(self):
        if os.path.exists(self.pid_file):
            kill_pid_file(self.pid_file)
            rm(self.pid_file)
        
        if self.output_level == 0:
            rm(self.output_dir)

    def __del__(self):
        self.exit()
