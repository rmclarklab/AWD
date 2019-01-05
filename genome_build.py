#python2 genome_build.py -e Output_files/window_mapping/mapmatrix.round5.txt -t Output_files/window_mapping/reciprocalpairs.txt -d Output_files/genome_order

import sys
import argparse

script_description = """
uses the recipricol txt file from the matrix build code to assemble the genome
"""

parser = argparse.ArgumentParser(
		description = script_description,
		formatter_class = argparse.RawDescriptionHelpFormatter)

parser.add_argument("-e", "--endchrom_file", required = True)			
parser.add_argument("-t", "--txt_file", required = True)
parser.add_argument("-d", "--directory", required = True)
args = parser.parse_args()

directory = args.directory
endchrom_file = args.endchrom_file
recippairs_file = args.txt_file

chromends_list = []
chromends_justscaff_list = []
used_chromend_list = []
pair_dict = {}
all_fix = []

outfile1 = open('%s/%s.txt' % (directory, "genomebuild"), 'w')
outfile2 = open('%s/%s.txt' % (directory, "genomebuildeasyread"), 'w')

with open(endchrom_file,"r") as openendchrom_file:
	for line in openendchrom_file:
		line = line.strip().split("\t")
		column_1 = line[0]
		column_1_justscaff = "_".join(column_1.split("_")[:-1])
		if "all" in column_1:
			column_1 = "_".join(column_1.split("_")[:-1])
			column_1 = "%s_left" % (column_1)
			all_fix.append(column_1)
		chromends_list.append(column_1)
		chromends_justscaff_list.append(column_1_justscaff)
		
with open(recippairs_file, "r") as open_pairs_file:
	for line in open_pairs_file:
		line = line.strip().split(",")
		pair_1 = line[0][2:-1]
		pair_2 = line[1][2:-2]
		pair_1_scaff = "_".join(pair_1.split("_")[:-1])
		pair_2_scaff = "_".join(pair_2.split("_")[:-1])
		if "all" in pair_1:
			pair_1_check = "%s_left" % (pair_1_scaff)
			if pair_1_check not in all_fix:
				pair_1 = pair_1_check
				all_fix.append(pair_1)
			else:
				pair_1 = "%s_right" % (pair_1_scaff)
				all_fix.append(pair_1)
		if "all" in pair_2:
			pair_2_check = "%s_left" % (pair_2_scaff)
			if pair_2_check not in all_fix:
				pair_2 = pair_2_check
				all_fix.append(pair_2)
			else:
				pair_2 = "%s_right" % (pair_2_scaff)
				all_fix.append(pair_2)
		if pair_1_scaff not in pair_dict:
			pair_dict[pair_1_scaff] = {}
		pair_dict[pair_1_scaff][pair_1] = pair_2
		if pair_2_scaff not in pair_dict:
			pair_dict[pair_2_scaff] = {}
		pair_dict[pair_2_scaff][pair_2] = pair_1


for chromend in chromends_list:
	chromend_justscaff = "_".join(chromend.split("_")[:-1])
	if chromend_justscaff in used_chromend_list:
		pass
	else:
		#print used_chromend_list
		chrom_build_list = []
		used_chromend_list.append(chromend_justscaff)
		half_1 = chromend
		half_justscaff = "_".join(half_1.split("_")[:-1])
		chromosome_building = True
		#print "outside_while", pair, pair_justscaff
		while chromosome_building == True:
			if half_justscaff not in pair_dict:
				print "does not exist"
				outfile1.write("%s\n" % (half_justscaff))
				outfile2.write("%s\n" % (half_justscaff))
				chromosome_building = False
			else:
				for half in pair_dict[half_justscaff]:
					if half == half_1:
						pass
					else:
						half_2 = half
						half_1_orient = half_1.split("_")[-1]
						half_2_orient = half_2.split("_")[-1]
						if "left" in half_1_orient or "left" in half_2_orient:
							whole = "%s_all" % (half_justscaff)
							chrom_build_list.append(whole)
						elif "beg" == half_1_orient and "end" == half_2_orient:
							whole = "%s_f" % (half_justscaff)
							chrom_build_list.append(whole)
						elif "end" == half_1_orient and "beg" == half_2_orient:
							whole = "%s_r" % (half_justscaff)
							chrom_build_list.append(whole)	
						match = pair_dict[half_justscaff][half_2]
						half_1 = match
						half_justscaff = "_".join(half_1.split("_")[:-1])
						if half_justscaff in chromends_justscaff_list:
							half_1_orient = half_1.split("_")[-1]
							if "left" in half_1_orient or "right" in half_1_orient:
								whole = "%s_all" % (half_justscaff)
								chrom_build_list.append(whole)
								used_chromend_list.append(half_justscaff)
								chromosome_building = False
							elif "beg" == half_1_orient:
								whole = "%s_f" % (half_justscaff)
								chrom_build_list.append(whole)
								used_chromend_list.append(half_justscaff)
								chromosome_building = False
							elif "end" == half_1_orient:
								whole = "%s_r" % (half_justscaff)
								chrom_build_list.append(whole)
								used_chromend_list.append(half_justscaff)
								chromosome_building = False
						break
		for entry in chrom_build_list:
			outfile1.write("%s\t" % entry)
			easyread = "_".join(entry.split("_")[1:])
			outfile2.write("%s - " % easyread)
		outfile1.write("end \n")
		outfile2.write("end \n")
		
outfile1.close()
outfile2.close()
