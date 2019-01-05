import argparse
import sys

script_description = """
Code will cord transform a fasta file based on an input file with
cord transform information.  The code has three functions listed below.
    1. Original fasta --> transformed fasta
        - default setting if nothing else is provided
    2. Original fasta --> transformed vcf file
        - vcf provided default
    3. Transformed cord --> original cord 
        - cords provided default
"""

parser = argparse.ArgumentParser(
         description = script_description,
         formatter_class = argparse.RawDescriptionHelpFormatter)
         
parser.add_argument("-f", "--fasta", required = True, help = "original fasta file from original genome being transformed. Format shoud be > scaffold_name following nucleotide sequence")
parser.add_argument("-i", "--input", required = True, help = "input file provides cord transform order.  It has 5 columns with tabs: out_scaff in_scaff direction beginning end")
parser.add_argument("-n", "--number_n", default = 0, help = "number of N's that go between the joins. Typically 0 for vcf transforms. Out fasta depends on users specification")
parser.add_argument("-c", "--cord", default = None, help = "provide transformed scaffold name of the position you wish to convert to the original fasta position")            
parser.add_argument("-p", "--position", default = None, help = "provide position on the out_scaff of interest to be converted to the original position in fasta")    
parser.add_argument("-v", "--vcf", default = None, help = "provide vcf file to transform SNP scaffold and position into new cord scheme")
parser.add_argument("-o", "--outdir", default = None, help = "provide directory for writing the vcf file")


args = parser.parse_args()
                
fasta = args.fasta                                                                    
input_file = args.input                                                                    
Numb_N = int(args.number_n)                                                                
out_scaff_cord = args.cord                                                            
out_scaff_cord_pos = args.position
vcf = args.vcf
outdir = args.outdir                    

order_dict = {}                                                                         
out_scaff_seq_dict = {}                                                                 
out_scaff_seq_write_dict = {}
out_scaff_order_list = []                                                                         
in_scaff_list = []
in_scaff_cord_dict = {}                                                                     
gap_seq_dict = {}                                                                             
gap_seq_write_dict = {}
vcf_input_info = {}
vcf_dict = {}
vcf_gap_dict = {}
vcf_uncalled_dict = {}
header_length_dict = {}
header_gap_dict = {}
header_gap_length_dict = {}
header_uncalled_dict = {}

def error(message, line):                                                             
    '''error function will occur at anytime the input file is in the incorrect format'''                                                                                 
    sys.exit("%s\n%s" % (message, line))
    
def warning(message, line):                                                             
    '''Error funtion for the scaffolds that are split up with parts left out; however it will not end the code it'll just give you a warning and continue'''
    sys.stderr.write("%s\n%s" % (message, line))

fasta_dict = {}
fasta_order_list= []
scaffold = None

with open (fasta, "r") as fasta_file:
    for line in fasta_file:
        line = line.strip()
        if len(line) > 0:
            if line[0] == ">":
                if scaffold != None:
                    fasta_dict[scaffold] = "".join(seq_list)
                scaffold = line.split()[0][1:]
                fasta_order_list.append(scaffold)
                seq_list = []
            else:
                seq_list.append(line)
    fasta_dict[scaffold] = "".join(seq_list)

with open(input_file, "r") as input_file:                                                                                        
    for line in input_file:                                                                 
        if line[0] != "#":                                                                     
            line = line.strip().split("\t")                                                 
            if len(line) != 5:                                                                 
                error("Error: 5 columns needed", "\t".join(line))
            else:                                                                             
                out_scaff = line[0]                                                         
                in_scaff = line[1]                                                             
                direction = line[2]                                                         
                start = line[3]                                                             
                end = line[4]                                                                 
                if direction not in ["f", "r"]:                                             
                    error("Unkown letter in direction column", "\t".join(line))
                if start == "b":                                                             
                    start = 1
                elif start.isdigit():                                                         
                    start = int(start)
                else:                                                                         
                    error("Error in the start column", "\t".join(line))
                if end == "e":                                                                  
                    end = len(fasta_dict[in_scaff])
                elif end.isdigit():
                    end = int(end)
                else:                                                                         
                    error("Andre you fucked up end", "\t".join(line))
            if out_scaff not in order_dict:                                                 
                order_dict[out_scaff] = []                                                     
                out_scaff_seq_dict[out_scaff] = []                                             
                out_scaff_seq_write_dict[out_scaff] = []
                out_scaff_order_list.append(out_scaff)                                              
                if vcf != None:
                    vcf_dict[out_scaff] = {}
                    header_length_dict[out_scaff] = {}
            order_dict[out_scaff].append([in_scaff, direction, start, end])                    
            if in_scaff not in in_scaff_cord_dict:                                             
                in_scaff_cord_dict[in_scaff] = []
                gap_seq_dict[in_scaff] = []
                gap_seq_write_dict[in_scaff] = []
                in_scaff_list.append(in_scaff)
                if vcf != None:
                    vcf_gap_dict[in_scaff] = {}
                    header_gap_dict[in_scaff] = {}
                    header_gap_length_dict[in_scaff] = 0
            in_scaff_cord_dict[in_scaff].append([start, end])                                 

for in_scaff in in_scaff_cord_dict:
    in_scaff_cord_dict[in_scaff].sort(key=lambda x: x[0])                                    
    prevend = 0
    for cord in in_scaff_cord_dict[in_scaff]:                                                
        start = cord[0]                                                                        
        if start <= prevend:                                                                
            error("Part of scaffolds overlap", "%s %s %s %s %s" % (in_scaff, "start =", start, "previous end =", prevend))
        if start > prevend + 1:
            warning("parts of a scaffold/s missing", 
            "%s%s\n" % ("broken scaffold=", in_scaff))
            mis_start = prevend +1
            mis_end = start -1
            mis_seq = fasta_dict[in_scaff][mis_start-1:mis_end]
            gap_seq_dict[in_scaff].append(mis_seq)
        prevend = cord[1]
    length = len(fasta_dict[in_scaff])
    if prevend != length:
        warning("parts of a scaffold/s missing", 
            "%s%s\n" % ("broken scaffold=", in_scaff))
        mis_end_seq = fasta_dict[in_scaff][prevend:length]
        gap_seq_dict[in_scaff].append(mis_end_seq)

if out_scaff_cord != None and out_scaff_cord_pos != None:
    out_scaff_cord_pos = int(out_scaff_cord_pos)
    length = 0
    remainder = out_scaff_cord_pos
    len = 0
    for in_scaff_info in order_dict[out_scaff_cord]:
        scaff_name = in_scaff_info[0]
        beg = in_scaff_info[2]
        end = in_scaff_info[3]
        if length < out_scaff_cord_pos:
            remainder = remainder - len 
            len = int((end+1) - beg)
            length += len
            cord_transform_scaff = scaff_name
            cord_transform_direction = in_scaff_info[1]
            cord_transform_beg = in_scaff_info[2]
            cord_transform_end = in_scaff_info[3]
    if cord_transform_direction == "f":
        cord = (cord_transform_beg + remainder) - 1
    else: 
        cord = (cord_transform_end - remainder) + 1
    print cord_transform_scaff, cord

elif vcf != None: 
    for scaffold in fasta_order_list:
        if scaffold not in in_scaff_cord_dict:
            vcf_uncalled_dict[scaffold] = {}
    if outdir == None:
         outfile_vcf = open('%s.%s.vcf' % (vcf.split(".")[0].split("/")[-1], "transformed"), 'w')
    else:
         outfile_vcf = open('%s/%s.%s.vcf' % (outdir, vcf.split(".")[0].split("/")[-1], "transformed"), 'w')         
    for out_scaff in out_scaff_order_list:
        length = 0
        for in_scaff_info in order_dict[out_scaff]:
            length_beg = in_scaff_info[2]
            length_end = in_scaff_info[3]
            length += ((length_end+1)-length_beg)
        header_length = "".join(["##contig=<ID=", out_scaff, ",length=", str(length), ">"])
        outfile_vcf.write("%s\n" % (header_length))
    for in_scaff in in_scaff_cord_dict:    
        prevend = 0
        for cord in in_scaff_cord_dict[in_scaff]:
            start = cord[0]    
            if start > prevend + 1:
                mis_start = prevend +1
                mis_end = start 
                header_gap_length_dict[in_scaff] += (mis_end - mis_start)
            prevend = cord[1]
        length = len(fasta_dict[in_scaff])
        if prevend != length:
            header_gap_length_dict[in_scaff] += ((length+1) - prevend)
        if header_gap_length_dict[in_scaff] != 0 :
            header_gap = "".join(["##contig=<ID=", in_scaff, ",length=", str(header_gap_length_dict[in_scaff]), ">"])
            outfile_vcf.write("%s\n" % (header_gap))
    for scaffold in fasta_order_list:
        if scaffold not in in_scaff_cord_dict:
            uncalled = len(fasta_dict[scaffold])     
            uncalled_header = "".join(["##contig=<ID=", scaffold, ",length=", str(uncalled), ">"])
            outfile_vcf.write("%s\n" % (uncalled_header))
    with open(vcf, "r") as openvcf: 
         for line in openvcf:
             if line[0:6]=="#CHROM":
                header = line
                outfile_vcf.write("%s" % (header))
             if line[0] != "#":
                SNP_line = line.strip().split("\t")
                SNP_scaff = SNP_line[0]
                SNP_pos = int(SNP_line[1])
                for out_scaff in out_scaff_order_list:
                    for in_scaff_info in order_dict[out_scaff]:
                        input_scaff = in_scaff_info[0]
                        input_direction = in_scaff_info[1]
                        input_beg = in_scaff_info[2]
                        input_end = in_scaff_info[3]
                        if SNP_scaff == input_scaff:
                            if input_beg <= SNP_pos <= input_end:
                                length = 0
                                for in_scaff_info in order_dict[out_scaff]:
                                    length_scaff = in_scaff_info[0]
                                    length_direction = in_scaff_info[1]
                                    length_beg = in_scaff_info[2]
                                    length_end = in_scaff_info[3]
                                    if length_scaff == SNP_scaff:
                                        if length_beg <= SNP_pos <= length_end:
                                            if length_direction == "f":
                                                SNP_transform = length + ((SNP_pos+1) - length_beg)
                                            else:
                                                SNP_transform = length + ((length_end+1) - SNP_pos)
                                    length += ((length_end+1)-length_beg)
                                SNP_line[0] = out_scaff
                                SNP_line[1] = str(SNP_transform)
                                transform_line = "\t".join(SNP_line)
                                vcf_dict[out_scaff][SNP_transform] = transform_line
                        if SNP_scaff not in in_scaff_cord_dict:
                            transform_line = "\t".join(SNP_line)
                            vcf_uncalled_dict[SNP_scaff][SNP_pos] = transform_line
                for in_scaff in in_scaff_cord_dict:
                    if in_scaff == SNP_scaff:
                        prevend = 0
                        for cord in in_scaff_cord_dict[in_scaff]:
                            start = cord[0]    
                            if start > prevend + 1:
                                mis_start = prevend +1
                                mis_end = start
                                if mis_start < SNP_pos < mis_end:
                                    transform_line = "\t".join(SNP_line)
                                    vcf_gap_dict[SNP_scaff][SNP_pos] = transform_line
                            prevend = cord[1]
                        length = len(fasta_dict[in_scaff])
                        if prevend != length:
                            if prevend < SNP_pos < length:
                                transform_line = "\t".join(SNP_line)    
                                vcf_gap_dict[SNP_scaff][SNP_pos] = transform_line            
    print "printing file"
    for out_scaff in out_scaff_order_list:
        for positions in sorted(vcf_dict[out_scaff]):
            outfile_vcf.write("%s\n" % vcf_dict[out_scaff][positions])
    for scaffold in fasta_order_list:
        if scaffold in in_scaff_list:
            if vcf_gap_dict[scaffold] != []:
                for positions in sorted(vcf_gap_dict[scaffold]):
                    outfile_vcf.write("%s\n" % vcf_gap_dict[scaffold][positions])
        if scaffold not in in_scaff_cord_dict:
            for positions in sorted(vcf_uncalled_dict[scaffold]):
                outfile_vcf.write("%s\n" % vcf_uncalled_dict[scaffold][positions])
    outfile_vcf.close()

else:
    from Bio.Seq import Seq
    Numb_N = int(Numb_N)
    for out_scaff in out_scaff_order_list:                                                         
        for in_scaff_info in order_dict[out_scaff]:                                             
            seq = fasta_dict[in_scaff_info[0]][in_scaff_info[2]-1:in_scaff_info[3]]     
            if "f" == in_scaff_info[1]:                                                         
                out_scaff_seq_dict[out_scaff].append(seq)                     
            else:
                out_scaff_seq_dict[out_scaff].append(str(Seq(seq).reverse_complement()))                                
        out_seq = (Numb_N*"N").join(out_scaff_seq_dict[out_scaff])                          
        out_scaff_seq_write_dict[out_scaff] = out_seq
    
    for in_scaff in gap_seq_dict:
        if gap_seq_dict[in_scaff] != []:
            out_gap_seq = (Numb_N*"N").join(gap_seq_dict[in_scaff])
            gap_seq_write_dict[in_scaff] = out_gap_seq
    if outdir == None:
        outfile_fasta = open('%s_%s.fa' % ('.'.join(fasta.split('.')[:-1]), "sample"), 'w')
    else:
        outfile_fasta = open('%s/%s_%s.fa' % (outdir, '.'.join(fasta.split("/")[-1].split('.')[:-1]), "sample"), 'w')
    for out_scaff in out_scaff_order_list:
        outfile_fasta.write(">%s\n" % out_scaff)
        for i in range(0,len(out_scaff_seq_write_dict[out_scaff]),60):
            outfile_fasta.write("%s\n" % out_scaff_seq_write_dict[out_scaff][i:i+60])
    for in_scaff in gap_seq_dict:
        if gap_seq_write_dict[in_scaff] != []:
            outfile_fasta.write(">%s\n" % in_scaff)
            for i in range(0,len(gap_seq_write_dict[in_scaff]),60):
                outfile_fasta.write("%s\n" % gap_seq_write_dict[in_scaff][i:i+60])
    for scaffold in fasta_order_list:
        if scaffold not in in_scaff_cord_dict:
            outfile_fasta.write(">%s\n" % scaffold)
            for i in range(0,len(fasta_dict[scaffold]),60):
                outfile_fasta.write("%s\n" % fasta_dict[scaffold][i:i+60])
    outfile_fasta.close()