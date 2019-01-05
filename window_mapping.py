# updates 4Jan18
# find accompanying shell script for command line example
import sys 
import argparse

script_description = """
Code slides across a vcf and grabs the windows at the beginning and end of 
a scaffold to calculate the frequency.  It then compares it to all other windows
at the beginning and end of scaffolds to determine the closest frequencies.
"""

parser = argparse.ArgumentParser(
		description = script_description,
		formatter_class = argparse.RawDescriptionHelpFormatter)
			
parser.add_argument("-v", "--vcf", required = True, help = "provide vcf file for SNP information")
parser.add_argument("-1", "--parent_1", required = True)
parser.add_argument("-2", "--parent_2", required = True)
parser.add_argument("-s", "--selected", required = True)
parser.add_argument("-u", "--unselected", required = True)
parser.add_argument("-c", "--undercoverage", required = True, default = 15)
parser.add_argument("-o", "--overcoverage", required = True, default = 150)
parser.add_argument("-w", "--window", required = True, default = 150000)
parser.add_argument("-m", "--move", required = True, default = 15000)
parser.add_argument("-i", "--SNPsinwindow", required = True, default = 20)
parser.add_argument("-d", "--directory", required = False) # default is the same folder as the VCF
parser.add_argument("-e", "--end_scaffold", required = True)
parser.add_argument("-n", "--no_call_scaffold")
args = parser.parse_args()

directory = args.directory
stop_scaffold = args.end_scaffold
vcf = args.vcf 
parent = args.parent_1 
parent2 = args.parent_2 
sel = args.selected 
unsel = args.unselected 
undercoverage = args.undercoverage #this is a low coverage filter based on the average number of reads per sample with the number being a percentage 10 which is 10%
overcoverage = args.overcoverage
slider_window = args.window #window size in bp Ex. 150,000bp
slider_move = args.move #window move in bp Ex. 50,000bp 
SNPsinwindow = args.SNPsinwindow
undercoverage = int(undercoverage) * .01
overcoverage = int(overcoverage) * .01
slider_window = int(slider_window)
slider_move = int(slider_move)
SNPsinwindow = int(SNPsinwindow)
vcf_split = '.'.join(vcf.split('.')[:-1])
stop_scaffold = float(stop_scaffold)
no_call_scaffold = args.no_call_scaffold
if no_call_scaffold != None:
	no_call_scaffold = float(no_call_scaffold)

parent_coverage = []
parent2_coverage = []
sel_coverage = []
unsel_coverage = []

with open(vcf, "r") as openvcf:
	for line in openvcf:
		if line[0:2] == "##": #this contains the lengths of the scaffold for the original pass I don't care about it. 
			pass
		elif line[0:6] == "#CHROM": #need to store all the variable names for the samples and which columns correspond to it. 
			header = line.strip().split("\t")
			header_strains = header[9:]
		else:
			if line != "#": #need average coverage for parent1, parent2, sel, and unsel
				SNP_line = line.strip().split("\t")	
				SNP_score = float(SNP_line[5])
				SNP_scaff = SNP_line[0]
				SNP_pos = SNP_line[1]
				ref_allele = SNP_line[3]
				alt_alleles = SNP_line[4].split(",")
				SNP_samples = SNP_line[9:]
				if len(ref_allele) == 1 and SNP_score > 100:
					if all(len(alt_allele) == 1 for alt_allele in alt_alleles):
						haplotypes = SNP_line[9:]
 						parent_num = header_strains.index(parent) 
 						parent_num2 = header_strains.index(parent2)
 						sel_num = header_strains.index(sel) 
 						unsel_num = header_strains.index(unsel)
 						
 						if ("./" in haplotypes[parent_num]) == False:
 							parent_reads = haplotypes[parent_num].split(":")[1].split(",")
 							parent_cov = sum([int(i) for i in parent_reads])
 							parent_coverage.append(parent_cov)
 						if ("./" in haplotypes[parent_num2]) == False:
							parent2_reads = haplotypes[parent_num2].split(":")[1].split(",")
							parent2_cov = sum([int(i) for i in parent2_reads])
							parent2_coverage.append(parent2_cov)
						if ("./" in haplotypes[sel_num]) == False:
							sel_reads = haplotypes[sel_num].split(":")[1].split(",")
							sel_cov = sum([int(i) for i in sel_reads])
							sel_coverage.append(sel_cov)
						if ("./" in haplotypes[unsel_num]) == False:
							unsel_reads = haplotypes[unsel_num].split(":")[1].split(",")
							unsel_cov = sum([int(i) for i in unsel_reads])
							unsel_coverage.append(unsel_cov)

							
parent_avg_coverage = int(sum(parent_coverage)/len(parent_coverage))
parent2_avg_coverage = int(sum(parent2_coverage)/len(parent2_coverage))
sel_avg_coverage = int(sum(sel_coverage)/len(sel_coverage))
unsel_avg_coverage = int(sum(unsel_coverage)/len(unsel_coverage))

over_parent = int(parent_avg_coverage * overcoverage)
over_parent2 = int(parent2_avg_coverage * overcoverage)
over_sel = int(sel_avg_coverage * overcoverage)
over_unsel = int(unsel_avg_coverage * overcoverage)

under_parent = int(parent_avg_coverage * undercoverage)
under_parent2 = int(parent2_avg_coverage * undercoverage)
under_sel = int(sel_avg_coverage * undercoverage)
under_unsel = int(unsel_avg_coverage * undercoverage)

with open(vcf, "r") as openvcf: 
	if directory == None: 
		outfile3 = open('%s.%s.%s.%s.map.txt' % ('.'.join(vcf_split.split('.')[:-1]), sel, slider_move, slider_window), 'w') 
		outfile4 = open('%s.%s.%s.%s.map.txt' % ('.'.join(vcf_split.split('.')[:-1]), unsel, slider_move, slider_window), 'w')
	else: 
		outdir = directory + "/"
		before = before = '.'.join(vcf.split('.')[:-1]).split("/")[-1]
		after = after = outdir + before
		outfile3 = open('%s.%s.%s.%s.map.txt' % (after, sel, slider_move, slider_window), 'w') 
		outfile4 = open('%s.%s.%s.%s.map.txt' % (after, unsel, slider_move, slider_window), 'w')
	scaff_length_dict = {} 
	scaff_length_dict = {} 
	sel_allele_SNP_dict = {} 
	sel_totcov_SNP_dict = {}
	unsel_allele_SNP_dict = {}
	unsel_totcov_SNP_dict = {}
	SNP_sel_slide_dict = {} #1st key is median SNP position in window, 2nd key is mean SNP frequency difference for the window
	SNP_unsel_slide_dict = {}
	scaff_range_dict = {}
	map_windows_sel_1 = {}
	map_windows_sel_2 = {}
	map_windows_unsel_1 = {}
	map_windows_unsel_2 = {}
	count = 0
	header_order_list = [] #header needs to get the order of the out scaffs and store it line by line in the list for later
	for line in openvcf:
		if line[0:2]=="##":  
			info = line.strip().split("=")  #line parsing - check
			scaff = info[2].split(",")[0]  # - check
			if 'mt' not in scaff: # if vcf contains mitochondrial DNA this will get rid of it. 
				scaff_number = float(scaff.split("_")[1])
				if scaff_number <= stop_scaffold and scaff_number != no_call_scaffold:
					header_order_list.append(scaff) #this will help print the outfile later. 
					length = int(info[3][:-1]) # - check
					scaff_length_dict[scaff] = int(length)
					sel_allele_SNP_dict[scaff] = {} 
					sel_totcov_SNP_dict[scaff] = {}
					unsel_allele_SNP_dict[scaff] = {}
					unsel_totcov_SNP_dict[scaff] = {}
					SNP_sel_slide_dict[scaff] = {}
					SNP_unsel_slide_dict[scaff] = {}
					scaff_range_dict[scaff] = {}
		elif line[0:6]=="#CHROM":
			header = line.strip().split("\t")
			header_strains = header[9:] 
		else:
			if line[0]!="#":
				SNP_line = line.strip().split("\t")
				SNP_scaff = SNP_line[0]
				if 'mt' not in SNP_scaff:
					SNP_scaff_number = float(SNP_scaff.split("_")[1])
					if SNP_scaff_number <= stop_scaffold and SNP_scaff_number != no_call_scaffold:
						SNP_pos = int(SNP_line[1])
						SNP_pos_bin = int(SNP_pos/slider_window)
						if SNP_pos_bin not in sel_allele_SNP_dict[SNP_scaff]:
							sel_allele_SNP_dict[SNP_scaff][SNP_pos_bin] = {}
							sel_totcov_SNP_dict[SNP_scaff][SNP_pos_bin] = {}
							unsel_allele_SNP_dict[SNP_scaff][SNP_pos_bin] = {}
							unsel_totcov_SNP_dict[SNP_scaff][SNP_pos_bin] = {}
						ref_allele = SNP_line[3]
						alt_alleles = SNP_line[4].split(",")
						SNP_score = float(SNP_line[5])
						if len(ref_allele) == 1 and SNP_score > 100: 
							if all(len(alt_allele) == 1 for alt_allele in alt_alleles): # per Robert
								haplotypes = SNP_line[9:]
								parent_num = header_strains.index(parent) #find the appropriate column in header for the parent
								parent_num2 = header_strains.index(parent2)
								sel_num = header_strains.index(sel) #find the appropriate column in header for the selected F1
								unsel_num = header_strains.index(unsel) #find the appropriate column in header for the unselected F1
								if ("./" in haplotypes[parent_num]) == False and ("./" in haplotypes[parent_num2]) == False and ("./" in haplotypes[unsel_num])==False and ("./" in haplotypes[sel_num])==False: 
									parent_reads = haplotypes[parent_num].split(":")[1].split(",") #reads into a list 
									parent2_reads = haplotypes[parent_num2].split(":")[1].split(",")
									sel_reads = haplotypes[sel_num].split(":")[1].split(",")
									unsel_reads = haplotypes[unsel_num].split(":")[1].split(",")
									parent_cov = sum([int(i) for i in parent_reads]) #sum reads to get total coverage for the SNP
									parent2_cov = sum([int(i) for i in parent2_reads])
									sel_cov = sum([int(i) for i in sel_reads])
									unsel_cov = sum([int(i) for i in unsel_reads])
									if parent_cov <= over_parent and parent_cov >= under_parent: # filter SNPs based on sample coverage 
										if parent2_cov <= over_parent2 and parent2_cov >= under_parent2:
											if sel_cov <= over_sel and sel_cov >= under_sel:
												if unsel_cov <= over_unsel and unsel_cov >= under_unsel:
													par_genotype = haplotypes[parent_num].split(":")[0].split("/")						
													par_genotype2 = haplotypes[parent_num2].split(":")[0].split("/")
													par_allele_1 = par_genotype[0] # value is either 0 or 1
													par_allele_2 = par_genotype[1]
													par2_allele_1 = par_genotype2[0] # value is either 0 or 1
													par2_allele_2 = par_genotype2[1] # value is either 0 or 1
													if par_allele_1 == par_allele_2 and par2_allele_1 == par2_allele_2 and par_allele_1 != par2_allele_1: #filter to check for fixation and that there isn't copy variation in parent
														sel_coverage = haplotypes[sel_num].split(":")[1].split(",")
														sel_allele_coverage = sel_coverage[int(par_allele_1)]
														sel_total_coverage = sum([int(i) for i in sel_coverage])
														sel_allele_SNP_dict[SNP_scaff][SNP_pos_bin][SNP_pos] = sel_allele_coverage
														sel_totcov_SNP_dict[SNP_scaff][SNP_pos_bin][SNP_pos] = sel_total_coverage
														unsel_coverage = haplotypes[unsel_num].split(":")[1].split(",")
														unsel_allele_coverage = unsel_coverage[int(par_allele_1)]
														unsel_total_coverage = sum([int(i) for i in unsel_coverage])
														unsel_allele_SNP_dict[SNP_scaff][SNP_pos_bin][SNP_pos] = unsel_allele_coverage
														unsel_totcov_SNP_dict[SNP_scaff][SNP_pos_bin][SNP_pos] = unsel_total_coverage	

small_scaffs_2 = []
for scaff in header_order_list:
	window_start = 0 + slider_move 
	window_end = slider_window + slider_move
	needed_length = (slider_window + slider_move + slider_move)
	if needed_length >= scaff_length_dict[scaff]:
		small_scaffs_2.append(scaff)
	else:
		window_median = (window_start + window_end)/2
		window_start_bin = int(window_start/slider_window)
		window_end_bin = int(window_end/slider_window)
		window_sel_allele_list = []
		window_sel_totcov_list = []
		window_unsel_allele_list = []
		window_unsel_totcov_list = []
		for SNP_pos_bin in sel_allele_SNP_dict[scaff]:
			if window_start_bin <= SNP_pos_bin <= window_end_bin:
				for SNP_pos in sel_allele_SNP_dict[scaff][SNP_pos_bin]:
					if window_start <= SNP_pos <= window_end:
						window_sel_allele_list.append(int(sel_allele_SNP_dict[scaff][SNP_pos_bin][SNP_pos]))
						window_sel_totcov_list.append(int(sel_totcov_SNP_dict[scaff][SNP_pos_bin][SNP_pos]))
						window_unsel_allele_list.append(int(unsel_allele_SNP_dict[scaff][SNP_pos_bin][SNP_pos]))
						window_unsel_totcov_list.append(int(unsel_totcov_SNP_dict[scaff][SNP_pos_bin][SNP_pos]))
		if len(window_sel_allele_list) > SNPsinwindow and len(window_unsel_allele_list) > SNPsinwindow:
			sel_SNPsinwindow = len(window_sel_allele_list)
			unsel_SNPsinwindow = len(window_unsel_allele_list)
			window_sel_allele_sum = sum(window_sel_allele_list)
			window_sel_totcov_sum = sum(window_sel_totcov_list)
			window_unsel_allele_sum = sum(window_unsel_allele_list)
			window_unsel_totcov_sum = sum(window_unsel_totcov_list)
			window_sel_freq = float(window_sel_allele_sum)/window_sel_totcov_sum
			window_unsel_freq = float(window_unsel_allele_sum)/window_unsel_totcov_sum
			map_windows_sel_2['%s_beg' % (scaff)] = ['%s_beg' % (scaff), window_median, window_sel_freq]
			map_windows_unsel_2['%s_beg' % (scaff)] = ['%s_beg' % (scaff), window_median, window_unsel_freq]
			outfile3.write('%s\t%s\t%s\t%s\n' % ('%s_beg' % (scaff), window_median, window_sel_freq, sel_SNPsinwindow))
			outfile4.write('%s\t%s\t%s\t%s\n' % ('%s_beg' % (scaff), window_median, window_unsel_freq, unsel_SNPsinwindow))
	
for scaff in header_order_list:
	if scaff in small_scaffs_2:
		pass 
	else:
		window_start = (scaff_length_dict[scaff] - slider_window - slider_move)
		window_end = (scaff_length_dict[scaff] - slider_move)
		window_median = (window_start + window_end)/2
		window_start_bin = int(window_start/slider_window)
		window_end_bin = int(window_end/slider_window)
		window_sel_allele_list = []
		window_sel_totcov_list = []
		window_unsel_allele_list = []
		window_unsel_totcov_list = []
		for SNP_pos_bin in sel_allele_SNP_dict[scaff]:
			if window_start_bin <= SNP_pos_bin <= window_end_bin:
				for SNP_pos in sel_allele_SNP_dict[scaff][SNP_pos_bin]:
					if window_start <= SNP_pos <= window_end:
						window_sel_allele_list.append(int(sel_allele_SNP_dict[scaff][SNP_pos_bin][SNP_pos]))
						window_sel_totcov_list.append(int(sel_totcov_SNP_dict[scaff][SNP_pos_bin][SNP_pos]))
						window_unsel_allele_list.append(int(unsel_allele_SNP_dict[scaff][SNP_pos_bin][SNP_pos]))
						window_unsel_totcov_list.append(int(unsel_totcov_SNP_dict[scaff][SNP_pos_bin][SNP_pos]))
		if len(window_sel_allele_list) > SNPsinwindow and len(window_unsel_allele_list) > SNPsinwindow:
			sel_SNPsinwindow = len(window_sel_allele_list)
			unsel_SNPsinwindow = len(window_unsel_allele_list)
			window_sel_allele_sum = sum(window_sel_allele_list)
			window_sel_totcov_sum = sum(window_sel_totcov_list)
			window_unsel_allele_sum = sum(window_unsel_allele_list)
			window_unsel_totcov_sum = sum(window_unsel_totcov_list)
			window_sel_freq = float(window_sel_allele_sum)/window_sel_totcov_sum
			window_unsel_freq = float(window_unsel_allele_sum)/window_unsel_totcov_sum
			map_windows_sel_2['%s_end' % (scaff)] = ['%s_end' % (scaff), window_median, window_sel_freq]
			map_windows_unsel_2['%s_end' % (scaff)] = ['%s_end' % (scaff), window_median, window_unsel_freq]
			outfile3.write('%s\t%s\t%s\t%s\n' % ('%s_end' % (scaff), window_median, window_sel_freq, sel_SNPsinwindow))
			outfile4.write('%s\t%s\t%s\t%s\n' % ('%s_end' % (scaff), window_median, window_unsel_freq, unsel_SNPsinwindow))

for scaff in small_scaffs_2:
	window_start = 0
	window_end = scaff_length_dict[scaff]
	window_median = (window_start + window_end)/2
	window_start_bin = int(window_start/slider_window)
	window_end_bin = int(window_end/slider_window)
	window_sel_allele_list = []
	window_sel_totcov_list = []
	window_unsel_allele_list = []
	window_unsel_totcov_list = []
	for SNP_pos_bin in sel_allele_SNP_dict[scaff]:
		if window_start_bin <= SNP_pos_bin <= window_end_bin:
			for SNP_pos in sel_allele_SNP_dict[scaff][SNP_pos_bin]:
				if window_start <= SNP_pos <= window_end:
					window_sel_allele_list.append(int(sel_allele_SNP_dict[scaff][SNP_pos_bin][SNP_pos]))
					window_sel_totcov_list.append(int(sel_totcov_SNP_dict[scaff][SNP_pos_bin][SNP_pos]))
					window_unsel_allele_list.append(int(unsel_allele_SNP_dict[scaff][SNP_pos_bin][SNP_pos]))
					window_unsel_totcov_list.append(int(unsel_totcov_SNP_dict[scaff][SNP_pos_bin][SNP_pos]))
	if len(window_sel_allele_list) > SNPsinwindow and len(window_unsel_allele_list) > SNPsinwindow:
		sel_SNPsinwindow = len(window_sel_allele_list)
		unsel_SNPsinwindow = len(window_unsel_allele_list)
		window_sel_allele_sum = sum(window_sel_allele_list)
		window_sel_totcov_sum = sum(window_sel_totcov_list)
		window_unsel_allele_sum = sum(window_unsel_allele_list)
		window_unsel_totcov_sum = sum(window_unsel_totcov_list)
		window_sel_freq = float(window_sel_allele_sum)/window_sel_totcov_sum
		window_unsel_freq = float(window_unsel_allele_sum)/window_unsel_totcov_sum
		map_windows_sel_2['%s_all' % (scaff)] = ['%s_all' % (scaff), window_median, window_sel_freq]
		map_windows_unsel_2['%s_all' % (scaff)] = ['%s_all' % (scaff), window_median, window_unsel_freq]
		outfile3.write('%s\t%s\t%s\t%s\n' % ('%s_all' % (scaff), window_median, window_sel_freq, sel_SNPsinwindow))
		outfile4.write('%s\t%s\t%s\t%s\n' % ('%s_all' % (scaff), window_median, window_unsel_freq, unsel_SNPsinwindow))

outfile3.close() 
outfile4.close()
