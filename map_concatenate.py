#python2 map_concatenate.py Output_files/window_mapping/*

import sys
import numpy as np

in_files = sys.argv[1:]

window = {}
window_diff = {} 
outfile1 = open('%s/mapmatrix.txt' % ('/'.join(in_files[0].split("/")[:-1])), 'w')

for file in in_files:
	with open(file, "r") as opentxt_file:
		for line in opentxt_file:
			line = line.strip().split("\t")
			scaff_pos = line[0]
			if scaff_pos not in window:
				window[scaff_pos] = []
				window_diff[scaff_pos] = []
			freq = line[2]
			window[scaff_pos].append(float(freq))

for scaff_pos in window:
	for scaff_pos_compare in window:
		scaff_pos_array = np.array(window[scaff_pos])
		scaff_pos_compare_array = np.array(window[scaff_pos_compare])
		diff = abs(scaff_pos_array - scaff_pos_compare_array)
		average_diff = sum(diff)/len(diff)
		window_diff[scaff_pos].append([scaff_pos_compare, average_diff])

for window in window_diff:
	outfile1.write('%s' % (window))
	for info in sorted(window_diff[window], key=lambda x: x[1])[1:]:
 		outfile1.write('\t%s' % (info))
 	outfile1.write('\n')
	
