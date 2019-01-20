echo "Transforming VCF"
python2 coordinate_code.py -f T_urticae_2009.09.28.fasta -i sangerbreaks.txt -v Tomato_Spirodiclofen_Joint.3.6-0-g89b7209.vcf 
echo "VCF transformed"
mkdir -p -- "Output_files/window_mapping"
echo "starting allele frequency calculations for each sample"
chmod a+x window_mapping_commands.sh
./window_mapping_commands.sh
echo "finished allele frequency calculations"
echo "building matrix"
python2 map_concatenate.py Output_files/window_mapping/*
python2 mapmatrix_build.py -t Output_files/window_mapping/mapmatrix.txt 
echo "writing genome order output"
mkdir -p -- "Output_files/genome_order"
python2 genome_build.py -e Output_files/window_mapping/mapmatrix.round5.txt -t Output_files/window_mapping/reciprocalpairs.txt -d Output_files/genome_order
