## AWD: Construction of superscaffolds using population allele frequency data. 

This collection of python programs links genomic scaffolds in draft genome sequences together to form superscaffolds using the average window distance (AWD) metric. This metric, and its use for joining together scaffolds in the draft genome sequence of the two-spotted spider mite (*Tetranychus urticae*) to produce three pseudochromosomes, is described in the following publication:

- Wybouw, N., Kosterlitz, O., Kurlovs, A. H., Bajda, S., Greenhalgh, R., Snoeck, S., Bui, H., Bryon, A., Dermauw, W., Van Leeuwen, T., and Clark, R. M. 2019. Long-term population studies uncover the genome structure and genetic basis of xenobiotic and host plant adaptation in the herbivore *Tetranychus urticae*. Under review (for a preprint, see bioRxiv 474064; doi: https://doi.org/10.1101/474064)

The collection of programs, and associated documentation, used to produce Table S4 and Table S5 in Wyboux et al. (2019) are provided (the instructions for how to use the software to generate Table S4 is given below). Although the programs were written for and tested on draft genome sequences from *T. urticae*, in concerte with population allele frequency data generated as part of the Wyboux et al. (2019) study, they were designed to work for related genome projects for which comparable data sets and input files are available.

For the data sets needed to test installation, and to replicate Table S4, please contact Richard Clark (richard.m.clark@utah.edu). The data sets used as input will be available for public download as soon as the manuscript is accepted and appears online.

The assembler.py program was written by Olivia Kosterlitz. She is currently a graduate student at the University of Washington (livkost@uw.edu). Please contact either Olivia Kosterlitz or Richard Clark (University of Utah, richard.m.clark@utah.edu) with questions about the software.

---

### Installation

Make sure version 2.7 of [python](https://www.python.org/download/releases/2.7/) is installed.

---

### Analysis for Wybouw et. al. 2019.

Download and make the resulting folder the working directory. In addition to the program files in this repository, ensure that the working directory contains (1) a reference genome assembly in FASTA format, (2) variant calls for parents of a cross and the segregating population data (as a VCF file), and (3) a file that includes information on potential misassemblies. 

The code was designed for and validated with draft Sanger genome sequence of the two-spotted spider mite (*Tetranychus urticae*) and replicated experimental population data, and the output is an ordered concatenation of Sanger scaffolds into superscaffolds. 

Specifically, the draft genome sequence for *T. urticae* can be downloaded as described in [Grbic *et. al.* (2011)](https://www.nature.com/articles/nature10640). The VCF file can be dowloaded from the [preprint](https://doi.org/10.1101/474064) that describes the study for which this software was developed. A file that includes potential misassemblies in the draft Sanger assembly of *T. urticae* is provided in this repository (`sangerbreaks.txt`). Please see the [preprint](https://doi.org/10.1101/474064) for additional information. 

```
cd Desktop/AWD-master/
```
Run the master shell script.
```
./master_script.sh
```
If the permission is denied, make the script executable. 
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

Examples of this command line are given in `window_mapping_commands.sh`

#### AWD calculations and super scaffold map construction 

Using `map_concatenate.py`, the average window distance (AWD) is calculated for every possible scaffold end pairing. To calculate an AWD, we (1) determine the window's allele frequency (previous script), (2) calculate the absolute value of the difference in the two frequency values between two non-overlapping windows on a per sample basis, and (3) average the resulting values across all samples. For genome construction, we (1) calculate allele frequency for only the terminal windows of each scaffold, (2) calculate AWD for all possible scaffold end pairings. For each scaffold end, a list is produced with all non-self AWD comparisons containing the scaffold end and sorted in ascending order. The five smallest AWD comparisons are retained for downstream scripts. Using `mapmatrix_build.py`, reciprical smallest AWD values matches two scaffold ends together. When this occurs, the scaffold ends are removed from the lists and excluded from the remaining rounds of matching. This process is iterative until there are no more rankings to compare. The remaining unmatched scaffold ends are used as chromosome end anchors in `genome_build.py`, which constructs the final superscaffold map by placing and ordering scaffolds according to the catalog of repicricol best hits. 


