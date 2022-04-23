
#!/usr/bin/env python

# from .base import *
from .readers import LorenzoReader
import numpy as np
import os.path
import sys
import subprocess
import tempfile
import json


def create_mappers(topologyfile):

	with open(topologyfile,'r') as file:
		content = file.read()

	nucleotides = content.split('\n')

	total_nucl = int(nucleotides[0].split(' ')[0])
	total_stra = int(nucleotides[0].split(' ')[1])

	strand_to_sequence = [None for x in range(total_stra+1)] # add one cause strand numbers start at one
	absolute_to_strand = [None for x in range(total_nucl)]
	strand_to_absolute = [[] for x in range(total_stra+1)]

	index = 1
	strand_num = 1
	first_nucl = 0

	while index < total_nucl:

		seq = ''
		data = nucleotides[index].split(' ')
		first_nucl = index
		
		while strand_num == int(data[0]):
			
			absolute_to_strand[index-1] = [ strand_num, index-first_nucl ]
			strand_to_absolute[strand_num].append( (index-1) )
		
			seq += data[1]
			index += 1
			if index <= total_nucl:
				data = nucleotides[index].split(' ')
			else:
				break
				
		strand_to_sequence[strand_num] = seq
		strand_num += 1
	
	return (strand_to_sequence, absolute_to_strand, strand_to_absolute)

def print_progress( current, total, first=False ):
	message = "Parsing through trajectory file :\t["
	for i in range(total):
		if i < current:
			message += "="
		else:
			message += " "
	message += "] {:.1f}%".format( 100*current/total )
	if not first:
		print("\b" * (200+total), end="")
	print(message, end="")
	if current == total:
		print()
	sys.stdout.flush()


def analyze_bonds(input_file, trajectory_file, topology_file, include_starting_bonds=False, oxDNA_dir=None):
	"""

	input_file -- oxDNA input file

	returns {
		'num_nuc': int number of nucleotides,
		'num_str': int number of strands,
		'str_seq': dict mapping strand numbers to sequences,
		'nuc_str': dict mapping base number to strand number,
		'str_nuc': dict mapping strand number to base number,
		'min_time': int starting time of system,
		'max_time': int ending time of system,
		'time_step': float suggested number of steps for each timestep,
		'bond_events': list of bond events (bonds breaking and forming)
	}
	"""
	#print(f'input_file = {input_file}, trajectory_file = {trajectory_file}, topology_file = {topology_file}')
	#print(f'include_starting_bonds = {include_starting_bonds}')

	with open(trajectory_file, 'r') as f:
		total_frames = f.read().count( 't = ' )

	reader = LorenzoReader(trajectory_file, topology_file)
	system = reader.get_system()

	if system == False:
		print("ERROR : Invalid trajectory file, no timesteps recorded")
		sys.exit(1)

	counter = 0
	
	temp_file = tempfile.NamedTemporaryFile()

	max_time = 0
	bond_events = []
	strand_to_sequence, base_to_strand, strand_to_base = create_mappers(topology_file)

	output = {
		"num_nuc": len(base_to_strand),
		"num_str": len(strand_to_sequence),
		"str_seq": strand_to_sequence,
		"nuc_str": base_to_strand,
		"str_nuc": strand_to_base,
		"min_time": system._time
	}
		
	# print_progress(0, total_frames, True)

	if not oxDNA_dir is None:
		DNAnalysis = f'{oxDNA_dir}/build/bin/DNAnalysis'
	else:
		DNAnalysis = 'DNAnalysis'

	args = [
		DNAnalysis,
		input_file,
		f'trajectory_file={temp_file.name}', 
		'analysis_data_output_1 = { \n name = stdout \n print_every = 1 \n col_1 = { \n type=pair_energy \n} \n}'
	]

	prev_bonds = []

	while system != False:


		system.map_nucleotides_to_strands()
		system.print_lorenzo_output(temp_file.name, '/dev/null')
		
		temp_file.flush()
		
		myinput = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		stdout, stderr = myinput.communicate()
		
		system.read_H_bonds(stdout.decode('utf-8').split('\n')[:-1])
		current_bonds = system.get_H_interactions_nucleotides()

		# if list of bonded nucleoties changes
		if current_bonds != prev_bonds and (counter > 0 or include_starting_bonds):
		
			# append all broken bonds to bond_events
			[ log(system, x, 'BREAK', base_to_strand, bond_events) for x in prev_bonds if x not in current_bonds ]
			
			# append all newly formed bonds to bond_events
			[ log(system, x, 'BINDS', base_to_strand, bond_events) for x in current_bonds if x not in prev_bonds ]
		
		prev_bonds = current_bonds
		max_time = system._time # used to record ending timestamp before system becomes False
		system = reader.get_system()
		counter += 1
		
		# print_progress(counter, total_frames)

	output['max_time'] = max_time
	output['time_step'] = (output['max_time'] - output['min_time']) / counter
	output['bond_events'] = bond_events

	return output



def log(system, bond_pair, action, base_to_strand, bond_events):
	""" records bond event in bond_events 
		only records the event if from a separate strand
	"""
	if system != None and bond_pair != None and base_to_strand != None:
		event = get_event(system, bond_pair, action, base_to_strand)
		if event[2] != event[-1]: # if different strands
			bond_events.append(event)
		# print( '\t\t'.join( [str(x) for x in event] ))


def get_event(system, bond_pair, action, base_to_strand):
	n1 = int(bond_pair[0])
	n2 = int(bond_pair[1])

	# time, nucleotide1, strand1, action, nucleotide2, strand2
	return [ system._time, n1, base_to_strand[n1][0], action, n2, base_to_strand[n2][0] ]





if __name__ == '__main__':

	if (len(sys.argv) < 3):
		print(f'Usage: python {sys.argv[0]} <input> <trajectory> <topology>')
		sys.exit()

	data = analyze_bonds(input_file=sys.argv[1], trajectory_file=sys.argv[2], topology_file=sys.argv[3])

	with open('bond_data.json','w') as f:
		json.dump(data, f)
