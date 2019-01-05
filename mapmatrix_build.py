#python2 mapmatrix_build.py -t Output_files/window_mapping/mapmatrix.txt 

import sys
import argparse

script_description = """
uses the matrix txt file to build the map
"""

parser = argparse.ArgumentParser(
		description = script_description,
		formatter_class = argparse.RawDescriptionHelpFormatter)
			
parser.add_argument("-t", "--txt_file", required = True)
args = parser.parse_args()

matrix_file = args.txt_file

nopairsleft = []
reference_dict = {}
scaff_pairs_used = []
allscaff_used = []
allscaff_duplicates = []
reciprocal_list = []
reciprocal_list_1 = []
reciprocal_list_2 = []
reciprocal_list_3 = []
reciprocal_list_4 = []
reciprocal_list_5 = []
allscaff_pairs = []
allscaff_used_round1 = []
allscaff_used_round2 = []
allscaff_used_round3 = []
allscaff_used_round4 = []
scaff_used_round1 = []
scaff_used_round2 = []
scaff_used_round3 = []
scaff_used_round4 = []
scaff_unused_round1 = []
scaff_unused_round2 = []
scaff_unused_round3 = []
scaff_unused_round4 = []
round1_emptyranks = []
round2_emptyranks = []
round3_emptyranks = []
round4_emptyranks = []

outfile1 = open('%s/reciprocalpairs.txt' % ('/'.join(matrix_file.split("/")[:-1])), 'w')
outfile2 = open('%s.round2.txt' % ('.'.join(matrix_file.split(".")[:-1])), 'w')
outfile3 = open('%s.round3.txt' % ('.'.join(matrix_file.split(".")[:-1])), 'w')
outfile4 = open('%s.round4.txt' % ('.'.join(matrix_file.split(".")[:-1])), 'w')
outfile5 = open('%s.round5.txt' % ('.'.join(matrix_file.split(".")[:-1])), 'w')

with open(matrix_file, "r") as openmatrix:
	for line in openmatrix:
		line = line.strip().split("\t")
		ref = line[0]
		ref_scaff = "_".join(ref.split("_")[:-1])
		if ref_scaff not in reference_dict:
			reference_dict[ref_scaff] = {}
		ref_orient = ref.split("_")[-1]
		reference_dict[ref_scaff][ref_orient] = line[1:7]

for ref_scaff in reference_dict:
	for ref_orient in reference_dict[ref_scaff]:
		rank_1_ref = reference_dict[ref_scaff][ref_orient][0]
		rank_2_ref = reference_dict[ref_scaff][ref_orient][1]
		ref_scaff_check = "%s_" % (ref_scaff)
		ref_scaff_combine = "%s_%s" % (ref_scaff, ref_orient)
		if ref_scaff_check in rank_1_ref:
			rank_1_ref = rank_2_ref		
 		for compare_scaff in reference_dict:
 			for compare_orient in reference_dict[compare_scaff]:
 				compare_combine = "%s_%s" % (compare_scaff, compare_orient)
 				if compare_combine in rank_1_ref:
 					rank_1_compare = reference_dict[compare_scaff][compare_orient][0]
					rank_2_compare = reference_dict[compare_scaff][compare_orient][1]
					if compare_scaff in rank_1_compare:
						rank_1_compare = rank_2_compare
					if ref_scaff_combine in rank_1_compare:
						if [compare_combine, ref_scaff_combine] in reciprocal_list or [ref_scaff_combine, compare_combine] in reciprocal_list:
							pass
						else:
							if "all" in ref_scaff_combine and "all" in compare_combine:
								if [compare_combine, ref_scaff_combine] not in reciprocal_list:
									allscaff_pairs.append([ref_scaff_combine, compare_combine])  
							reciprocal_list.append([ref_scaff_combine, compare_combine])
							reciprocal_list_1.append([ref_scaff_combine, compare_combine])
							scaff_pairs_used.append([ref_scaff, compare_scaff])
							if "all" not in ref_scaff_combine:
								scaff_used_round1.append(ref_scaff_combine)
							else: 
								allscaff_used_round1.append(ref_scaff_combine)
								allscaff_used.append(ref_scaff_combine)
							if "all" not in compare_combine:
								scaff_used_round1.append(compare_combine)
							else: 
								allscaff_used_round1.append(compare_combine)
								allscaff_used.append(compare_combine)
					else: 
						scaff_unused_round1.append(ref_scaff_combine)

for pairs in reciprocal_list_1:
	outfile1.write('%s\n' % (pairs))	
					
round2_dict = {}

with open(matrix_file, "r") as openmatrix:
	for line in openmatrix:
		line = line.strip().split("\t")
		column_1 = line[0]
		column_1_scaff = "_".join(column_1.split("_")[:-1])
		column_1_scaff_check = "%s_" % (column_1_scaff)
		column_1_orient = column_1.split("_")[-1]
		if column_1 in scaff_used_round1:
			pass
		elif column_1 in allscaff_used_round1:
			all_count = allscaff_used.count(column_1)
			if all_count == 2:
				pass
			else:
				if column_1_scaff not in round2_dict:
					round2_dict[column_1_scaff] = {}
				round2_dict[column_1_scaff][column_1_orient] = []
				for info in line[1:6]:
					info_split = info.split(",")[0][2:-1]
					if info_split not in scaff_used_round1:
						if column_1_scaff not in info_split:
							if [column_1, info_split] in allscaff_pairs or [info_split, column_1] in allscaff_pairs:
								pass
							else:
								round2_dict[column_1_scaff][column_1_orient].append(info)
				allscaff_duplicates.append(column_1_scaff)
 		else:
 			if column_1_scaff not in round2_dict:
 				round2_dict[column_1_scaff] = {}
			round2_dict[column_1_scaff][column_1_orient] = []
			for info in line[1:7]:
				info_split = info.split(",")[0][2:-1]
				if info_split not in scaff_used_round1:
					if column_1_scaff not in info_split:
						if "all" in info_split:
							all_count = allscaff_used.count(info_split)
							if all_count == 2:
								pass
							else:
								round2_dict[column_1_scaff][column_1_orient].append(info)
						else:
							round2_dict[column_1_scaff][column_1_orient].append(info)

							

for key in round2_dict:
	for key_2 in round2_dict[key]:
		outfile2.write("%s_%s\t%s\n" % (key, key_2, "\t".join(round2_dict[key][key_2])))
outfile2.close()

for round2_scaff in round2_dict:
	for round2_orient in round2_dict[round2_scaff]:
		round2_scaff_combine = "%s_%s" % (round2_scaff, round2_orient)
		if round2_dict[round2_scaff][round2_orient] != []:
			rank_1_round2 = round2_dict[round2_scaff][round2_orient][0]
		else:
			round1_emptyranks.append(round2_scaff)
			nopairsleft.append(round2_scaff_combine)
			break
 		for round2_compare_scaff in round2_dict:
 			for round2_compare_orient in round2_dict[round2_compare_scaff]:
 				round2_compare_combine = "%s_%s" % (round2_compare_scaff, round2_compare_orient)
 				if round2_dict[round2_compare_scaff][round2_compare_orient] != []:
					rank_1_round2_compare = round2_dict[round2_compare_scaff][round2_compare_orient][0]
				else:
					pass
 				if round2_compare_combine in rank_1_round2:
 					if round2_scaff_combine in rank_1_round2_compare:
						if [round2_compare_combine, round2_scaff_combine] in reciprocal_list or [round2_scaff_combine, round2_compare_combine] in reciprocal_list:
							pass
						elif [round2_scaff, round2_compare_scaff] in scaff_pairs_used or [round2_compare_scaff, round2_scaff] in scaff_pairs_used:
							if len(round2_dict[round2_scaff][round2_orient]) <= 1:
								pass
							else: 
								rank_1_round2 = round2_dict[round2_scaff][round2_orient][1]
								round2_scaff_combine = "%s_%s" % (round2_scaff, round2_orient)
								for round2_compare_scaff in round2_dict:
									for round2_compare_orient in round2_dict[round2_compare_scaff]:
										round2_compare_combine = "%s_%s" % (round2_compare_scaff, round2_compare_orient)
										if round2_dict[round2_compare_scaff][round2_compare_orient] != []:
											rank_1_round2_compare = round2_dict[round2_compare_scaff][round2_compare_orient][0]
										else:
											pass
										if round2_compare_combine in rank_1_round2:
											if round2_scaff_combine in rank_1_round2_compare:
												if [round2_compare_combine, round2_scaff_combine] in reciprocal_list or [round2_scaff_combine, round2_compare_combine] in reciprocal_list:
													pass
												else: 
													if "all" in round2_scaff_combine and "all" in round2_compare_combine:
														if [round2_compare_combine, round2_scaff_combine] not in reciprocal_list:
															allscaff_pairs.append([round2_scaff_combine, round2_compare_combine])  
													reciprocal_list.append([round2_scaff_combine, round2_compare_combine])
													reciprocal_list_2.append([round2_scaff_combine, round2_compare_combine])
													scaff_pairs_used.append([round2_scaff, round2_compare_scaff])	
													if "all" not in round2_scaff_combine:
														scaff_used_round2.append(round2_scaff_combine)
													else: 
														allscaff_used_round2.append(round2_scaff_combine)
														allscaff_used.append(round2_scaff_combine)
													if "all" not in round2_compare_combine:
														scaff_used_round2.append(round2_compare_combine)
													else: 
														allscaff_used_round2.append(round2_compare_combine)
														allscaff_used.append(round2_compare_combine)																		
						else: 
							if "all" in round2_scaff_combine and "all" in round2_compare_combine:
								if [round2_compare_combine, round2_scaff_combine] not in reciprocal_list:
									allscaff_pairs.append([round2_scaff_combine, round2_compare_combine])  
							reciprocal_list.append([round2_scaff_combine, round2_compare_combine])
							reciprocal_list_2.append([round2_scaff_combine, round2_compare_combine])
							scaff_pairs_used.append([round2_scaff, round2_compare_scaff])
							if "all" not in round2_scaff_combine:
								scaff_used_round2.append(round2_scaff_combine)
							else: 
								allscaff_used_round2.append(round2_scaff_combine)
								allscaff_used.append(round2_scaff_combine)
							if "all" not in round2_compare_combine:
								scaff_used_round2.append(round2_compare_combine)
							else: 
								allscaff_used_round2.append(round2_compare_combine)
								allscaff_used.append(round2_compare_combine)
					else: 
						scaff_unused_round2.append(round2_scaff_combine)

round3_dict = {}

outfile2_read = open('%s.round2.txt' % ('.'.join(matrix_file.split(".")[:-1])), 'r')

for pairs in reciprocal_list_2:
	outfile1.write('%s\n' % (pairs))	

for line in outfile2_read:
	line = line.strip().split("\t")
	column_1 = line[0]
	column_1_scaff = "_".join(column_1.split("_")[:-1])
	column_1_scaff_check = "%s_" % (column_1_scaff)
	column_1_orient = column_1.split("_")[-1]
	if column_1 in scaff_used_round2 or column_1 in nopairsleft:
		pass
	elif column_1 in allscaff_used:
		all_count = allscaff_used.count(column_1)
		if all_count == 2:
			pass
		else:
			if column_1_scaff not in round3_dict:
				round3_dict[column_1_scaff] = {}
			round3_dict[column_1_scaff][column_1_orient] = []
			for info in line[1:7]:
				info_split = info.split(",")[0][2:-1]
				if info_split not in scaff_used_round2 and info_split not in nopairsleft:
					if column_1_scaff not in info_split:
						if [column_1, info_split] in allscaff_pairs or [info_split, column_1] in allscaff_pairs:
							pass
						else:
							round3_dict[column_1_scaff][column_1_orient].append(info)
			allscaff_duplicates.append(column_1_scaff)
	else:
		if column_1_scaff not in round3_dict:
			round3_dict[column_1_scaff] = {}
		round3_dict[column_1_scaff][column_1_orient] = []
		for info in line[1:7]:
			info_split = info.split(",")[0][2:-1]
			if info_split not in scaff_used_round2 and info_split not in nopairsleft:
				if column_1_scaff not in info_split:
					if "all" in info_split:
						all_count = allscaff_used.count(info_split)
						if all_count == 2:
							pass
						else:
							round3_dict[column_1_scaff][column_1_orient].append(info)
					else:
						round3_dict[column_1_scaff][column_1_orient].append(info)


for key in round3_dict:
	for key_2 in round3_dict[key]:
		outfile3.write("%s_%s\t%s\n" % (key, key_2, "\t".join(round3_dict[key][key_2])))
outfile3.close()

for round3_scaff in round3_dict:
	for round3_orient in round3_dict[round3_scaff]:
		round3_scaff_combine = "%s_%s" % (round3_scaff, round3_orient)
		if round3_dict[round3_scaff][round3_orient] != []:
			rank_1_round3 = round3_dict[round3_scaff][round3_orient][0]
		else:
			round2_emptyranks.append(round3_scaff)
			nopairsleft.append(round3_scaff_combine)
			break
 		for round3_compare_scaff in round3_dict:
 			for round3_compare_orient in round3_dict[round3_compare_scaff]:
 				round3_compare_combine = "%s_%s" % (round3_compare_scaff, round3_compare_orient)
 				if round3_dict[round3_compare_scaff][round3_compare_orient] != []:
					rank_1_round3_compare = round3_dict[round3_compare_scaff][round3_compare_orient][0]
				else:
					break
 				if round3_compare_combine in rank_1_round3:
 					if round3_scaff_combine in rank_1_round3_compare:
						if [round3_compare_combine, round3_scaff_combine] in reciprocal_list or [round3_scaff_combine, round3_compare_combine] in reciprocal_list:
							pass
						elif [round3_scaff, round3_compare_scaff] in scaff_pairs_used or [round3_compare_scaff, round3_scaff] in scaff_pairs_used:
							rank_1_round3 = round3_dict[round3_scaff][round3_orient][1]
							round3_scaff_combine = "%s_%s" % (round3_scaff, round3_orient)
							for round3_compare_scaff in round3_dict:
								for round3_compare_orient in round3_dict[round3_compare_scaff]:
									round3_compare_combine = "%s_%s" % (round3_compare_scaff, round3_compare_orient)
									if round3_dict[round3_compare_scaff][round3_compare_orient] != []:
										rank_1_round3_compare = round3_dict[round3_compare_scaff][round3_compare_orient][0]
									else:
										pass
									if round3_compare_combine in rank_1_round3:
										if round3_scaff_combine in rank_1_round3_compare:
											if [round3_compare_combine, round3_scaff_combine] in reciprocal_list or [round3_scaff_combine, round3_compare_combine] in reciprocal_list:
												pass
											else: 
												if "all" in round3_scaff_combine and "all" in round3_compare_combine:
													if [round3_compare_combine, round3_scaff_combine] not in reciprocal_list:
														allscaff_pairs.append([round3_scaff_combine, round3_compare_combine])  
												reciprocal_list.append([round3_scaff_combine, round3_compare_combine])
												reciprocal_list_3.append([round3_scaff_combine, round3_compare_combine])
												scaff_pairs_used.append([round3_scaff, round3_compare_scaff])
												if "all" not in round3_scaff_combine:
													scaff_used_round3.append(round3_scaff_combine)
												else: 
													allscaff_used_round3.append(round3_scaff_combine)
													allscaff_used.append(round3_scaff_combine)
												if "all" not in round3_compare_combine:
													scaff_used_round3.append(round3_compare_combine)
												else: 
													allscaff_used_round3.append(round3_compare_combine)
													allscaff_used.append(round3_compare_combine)
						else: 
							if "all" in round3_scaff_combine and "all" in round3_compare_combine:
								if [round3_compare_combine, round3_scaff_combine] not in reciprocal_list:
									allscaff_pairs.append([round3_scaff_combine, round3_compare_combine])  
							reciprocal_list.append([round3_scaff_combine, round3_compare_combine])
							reciprocal_list_3.append([round3_scaff_combine, round3_compare_combine])
							scaff_pairs_used.append([round3_scaff, round3_compare_scaff])
							if "all" not in round3_scaff_combine:
 								scaff_used_round3.append(round3_scaff_combine)
 							else: 
								allscaff_used_round3.append(round3_scaff_combine)
								allscaff_used.append(round3_scaff_combine)
							if "all" not in round3_compare_combine:
								scaff_used_round3.append(round3_compare_combine)
							else: 
								allscaff_used_round3.append(round3_compare_combine)
								allscaff_used.append(round3_compare_combine)
					else: 
						scaff_unused_round3.append(round3_scaff_combine)

for pairs in reciprocal_list_3:
	outfile1.write('%s\n' % (pairs))

round4_dict = {}

outfile3_read = open('%s.round3.txt' % ('.'.join(matrix_file.split(".")[:-1])), 'r')

for line in outfile3_read:
	line = line.strip().split("\t")
	column_1 = line[0]
	column_1_scaff = "_".join(column_1.split("_")[:-1])
	column_1_scaff_check = "%s_" % (column_1_scaff)
	column_1_orient = column_1.split("_")[-1]
	if column_1 in scaff_used_round3 or column_1 in nopairsleft:
		pass
	elif column_1 in allscaff_used:
		all_count = allscaff_used.count(column_1)
		if all_count == 2:
			pass
		else:
			if column_1_scaff not in round4_dict:
				round4_dict[column_1_scaff] = {}
			round4_dict[column_1_scaff][column_1_orient] = []
			for info in line[1:7]:
				info_split = info.split(",")[0][2:-1]
				if info_split not in scaff_used_round3 and info_split not in nopairsleft:
					if column_1_scaff not in info_split:
						if [column_1, info_split] in allscaff_pairs or [info_split, column_1] in allscaff_pairs:
							pass
						else:
							round4_dict[column_1_scaff][column_1_orient].append(info)
			allscaff_duplicates.append(column_1_scaff)
	else:
		if column_1_scaff not in round4_dict:
			round4_dict[column_1_scaff] = {}
		round4_dict[column_1_scaff][column_1_orient] = []
		for info in line[1:7]:
			info_split = info.split(",")[0][2:-1]
			if info_split not in scaff_used_round3 and info_split not in nopairsleft:
				if column_1_scaff not in info_split:
					if "all" in info_split:
						all_count = allscaff_used.count(info_split)
						if all_count == 2:
							pass
						else:
							round4_dict[column_1_scaff][column_1_orient].append(info)
					else:
						round4_dict[column_1_scaff][column_1_orient].append(info)

for key in round4_dict:
	for key_2 in round4_dict[key]:
		outfile4.write("%s_%s\t%s\n" % (key, key_2, "\t".join(round4_dict[key][key_2])))
outfile4.close()

for round4_scaff in round4_dict:
	for round4_orient in round4_dict[round4_scaff]:
		round4_scaff_combine = "%s_%s" % (round4_scaff, round4_orient)
		if round4_dict[round4_scaff][round4_orient] != []:
			rank_1_round4 = round4_dict[round4_scaff][round4_orient][0]
		else:
			round3_emptyranks.append(round4_scaff)
			nopairsleft.append(round4_scaff_combine)
			break
 		for round4_compare_scaff in round4_dict:
 			for round4_compare_orient in round4_dict[round4_compare_scaff]:
 				round4_compare_combine = "%s_%s" % (round4_compare_scaff, round4_compare_orient)
 				if round4_dict[round4_compare_scaff][round4_compare_orient] != []:
					rank_1_round4_compare = round4_dict[round4_compare_scaff][round4_compare_orient][0]
				else:
					break
 				if round4_compare_combine in rank_1_round4:
 					if round4_scaff_combine in rank_1_round4_compare:
						if [round4_compare_combine, round4_scaff_combine] in reciprocal_list:
							pass
						elif [round4_scaff, round4_compare_scaff] in scaff_pairs_used or [round4_compare_scaff, round4_scaff] in scaff_pairs_used:
							if len(round4_dict[round4_scaff][round4_orient]) == 1:
								pass
							else: 
								rank_1_round4 = round4_dict[round4_scaff][round4_orient][1]
								round4_scaff_combine = "%s_%s" % (round4_scaff, round4_orient)
								for round4_compare_scaff in round4_dict:
									for round4_compare_orient in round4_dict[round4_compare_scaff]:
										round4_compare_combine = "%s_%s" % (round4_compare_scaff, round4_compare_orient)
										if round4_dict[round4_compare_scaff][round4_compare_orient] != []:
											rank_1_round4_compare = round4_dict[round4_compare_scaff][round4_compare_orient][0]
										else:
											pass
										if round4_compare_combine in rank_1_round4:
											if round4_scaff_combine in rank_1_round4_compare:
												if [round4_compare_combine, round4_scaff_combine] in reciprocal_list:
													pass
												else: 
													if "all" in round4_scaff_combine and "all" in round4_compare_combine:
														if [round4_compare_combine, round4_scaff_combine] not in reciprocal_list:
															allscaff_pairs.append([round4_scaff_combine, round4_compare_combine])  
													reciprocal_list.append([round4_scaff_combine, round4_compare_combine])
													reciprocal_list_4.append([round4_scaff_combine, round4_compare_combine])
													scaff_pairs_used.append([round4_scaff, round4_compare_scaff])
													if "all" not in round4_scaff_combine:
														scaff_used_round4.append(round4_scaff_combine)
													else: 
														allscaff_used_round4.append(round4_scaff_combine)
														allscaff_used.append(round4_scaff_combine)
													if "all" not in round4_compare_combine:
														scaff_used_round4.append(round4_compare_combine)
													else: 
														allscaff_used_round4.append(round4_compare_combine)
														allscaff_used.append(round4_compare_combine)
						else: 
							if "all" in round4_scaff_combine and "all" in round4_compare_combine:
								if [round4_compare_combine, round4_scaff_combine] not in reciprocal_list:
									allscaff_pairs.append([round4_scaff_combine, round4_compare_combine])  
							reciprocal_list.append([round4_scaff_combine, round4_compare_combine])
							reciprocal_list_4.append([round4_scaff_combine, round4_compare_combine])
							scaff_pairs_used.append([round4_scaff, round4_compare_scaff])
							if "all" not in round4_scaff_combine:
 								scaff_used_round4.append(round4_scaff_combine)
 							else: 
								allscaff_used_round4.append(round4_scaff_combine)
								allscaff_used.append(round4_scaff_combine)
							if "all" not in round4_compare_combine:
								scaff_used_round4.append(round4_compare_combine)
							else: 
								allscaff_used_round4.append(round4_compare_combine)
								allscaff_used.append(round4_compare_combine)
					else: 
						scaff_unused_round4.append(round4_scaff_combine)

for pairs in reciprocal_list_4:
	outfile1.write('%s\n' % (pairs))

round5_dict = {}

outfile4_read = open('%s.round4.txt' % ('.'.join(matrix_file.split(".")[:-1])), 'r')

for line in outfile4_read:
	line = line.strip().split("\t")
	column_1 = line[0]
	column_1_scaff = "_".join(column_1.split("_")[:-1])
	column_1_scaff_check = "%s_" % (column_1_scaff)
	column_1_orient = column_1.split("_")[-1]
	if column_1 in scaff_used_round4 or column_1 in nopairsleft:
		pass
	elif column_1 in allscaff_used:
		all_count = allscaff_used.count(column_1)
		if all_count == 2:
			pass
		else:
			if column_1_scaff not in round5_dict:
				round5_dict[column_1_scaff] = {}
			round5_dict[column_1_scaff][column_1_orient] = []
			for info in line[1:6]:
				info_split = info.split(",")[0][2:-1]
				if info_split not in scaff_used_round4 and info_split not in nopairsleft:
					if column_1_scaff not in info_split:
						if [column_1, info_split] in allscaff_pairs or [info_split, column_1] in allscaff_pairs:
							pass
						else:
							round5_dict[column_1_scaff][column_1_orient].append(info)
			allscaff_duplicates.append(column_1_scaff)
	else:
		if column_1_scaff not in round5_dict:
			round5_dict[column_1_scaff] = {}
		round5_dict[column_1_scaff][column_1_orient] = []
		for info in line[1:6]:
			info_split = info.split(",")[0][2:-1]
			if info_split not in scaff_used_round4 and info_split not in nopairsleft:
				if column_1_scaff not in info_split:
					if "all" in info_split:
						all_count = allscaff_used.count(info_split)
						if all_count == 2:
							pass
						else:
							round5_dict[column_1_scaff][column_1_orient].append(info)
					else:
						round5_dict[column_1_scaff][column_1_orient].append(info)

for key in round5_dict:
	for key_2 in round5_dict[key]:
		outfile5.write("%s_%s\t%s\n" % (key, key_2, "\t".join(round5_dict[key][key_2])))
for key in nopairsleft:
	outfile5.write("%s\tnopairslist\n" % (key))


outfile1.close()
outfile2.close()
outfile3.close()
outfile4.close()
outfile5.close()











