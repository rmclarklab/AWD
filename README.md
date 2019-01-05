## AWD: Construction of superscaffolds using population allele frequency data. 

This set of python scripts links genomic scaffolds into super scaffolds using the average window distance (AWD) metric. 

---

### Installation

Make sure version 2.7 of [python](https://www.python.org/download/releases/2.7/) is installed.

---

### Analysis from Wybouw et. al. 2019.

Download and unzip. Make the folder the working directory in command line. In addition to the scripts, make sure the working directory contains (1) Reference genome assembly in FASTA format (download from FILL THIS IN), (2) VCF file with variant predictions (download from FILL THIS IN), (3) a file that includes information on potential misassemblies (`sangerbreaks.txt` in this case). 
```
cd Desktop/AWD-master/
```
Run the master shell script.
```
./master_script.sh
```
If the permission is denied, then make the script executable using command line. 
```
chmod a+x master_script.sh
```
Print statements will appear in terminal designating which portion is running.

When completed, a new folder "Output_Files" will contain two folders "window_mapping" and "genome_order". The main output is "genome_order" folder, but all the intermediate files from the shell script can be found in "window_mapping" folder. 

### Usage. 

Currently, the `master_script.sh` runs through the following steps. 

#### VCF transformation

Using the `coordinate_code.py`, the variant call file (VCF) is transformed to include misassembly information, which breaks scaffolds into sub-scaffolds at specified break points. 

The script requires the following inputs:
- -f / --fasta: Fasta file used to create the VCF
- -i / --input: Input file providing coordinates to transform. The `sangerbreaks.txt` file was used for the *T.urticae* analysis.
- -v / --vcf: VCF file intended to transform the scaffold and SNP position. 
Here's an example using the *T.urticae* data. 
```
python2 coordinate_code.py -f T_urticae_2009.09.28.fasta -i sangerbreaks.txt -v Tomato_Spirodiclofen_Joint.3.6-0g89b7209.vcf
```

#### Calculate allele frequencies at terminal scaffold ends

Using the `window_mapping.py`, the user loads VCF data and indicates the experimental samples to calculate allele frequencies at terminal scaffold ends (two windows per scaffold). The script requires the following input:
- -v / --vcf: VCF file containing samples of interest

Control of the filters and calculations is through the following arguments:
- -1 / --parent_1
- -2 / --parent_2
- -s / --selected
- -u / --unselected
  - hint: the first four arguments are the sample names which match the VCF header. 
- -c / --undercoverage
- -o / --overcoverage
- -w / --window
- -m / --move
- -i / --SNPsinwindow
- -d / --directory
- -e / --end_scaffold
- -n / --no_call_scaffold

An example of this command line is given in `window_mapping_commands.sh`

#### AWD calculations and super scaffold map construction 

Using `map_concatenate.py`, the average window distance (AWD) is calculated for every possible scaffold end pairing. To calculate an AWD, we (1) determine the window's allele frequency (previous script), (2) calculate the absolute value of the difference in the two frequency values between two non-overlapping windows on a per sample basis, and (3) average the resulting values across all samples. For genome construction, we (1) calculate allele frequency for only the terminal windows of each scaffold, (2) calculate AWD for all possible scaffold end pairings. For each scaffold end, a list is produced with all non-self AWD comparisons containing the scaffold end and sorted in ascending order. The five smallest AWD comparisons are retained for downstream scripts. Using `mapmatrix_build.py`, reciprical smallest AWD values matches two scaffold ends together. When this occurs, the scaffold ends are removed from the lists and excluded from the remaining rounds of matching. This process is iterative until there are no more rankings to compare. The remaining unmatched scaffold ends are used as chromosome end anchors in `genome_build.py`, which constructs the final superscaffold map by placing and ordering scaffolds according to the catalog of repicricol best hits. 


