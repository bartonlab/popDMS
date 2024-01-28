The following instruction is incorporated with the example of the replicate #1 of YAP1 as the target protein.

To use the epistasis inference, you can run the bash file run_epistasis.sh under the folder of ./epistasis_inference.

Some explanations about the the codes in bash file:

"g++ -std=c++11 -lgslcblas -lgsl -I ./eigen-3.4.0/ get_freq.cpp -o get_freq": compile the C++ script about frequency extraction from the raw data files. 

"./get_freq [Target-protein] [genoype counts data file name] [target protein indexing file name] [raw sequence file]", which in this analysis is "./get_freq YAP1 ./YAP1_haplotype_count_rep1.csv index_matrix.csv ./YAP1_Reference_AA_sequence.txt.": collect intermediate frequency data files. This command will generate the following files: [Target-protein]_freq_[replicate num]_[generation_num].csv and [Target-protein]_multiple_allele_[replicate_num]_[generation_num].csv within the same directory you executed the ./get_freq command line.

“g++ -std=c++11 -lgslcblas -lgsl -I ./eigen-3.4.0/ inversion.cpp -o inversion”: compiling the inversion.cpp to estimate the fitness for single replicate with intermediate files.

“./inversion [Target-protein] [genotype counts data file] [target protein indexing file name] [raw sequence file] [regularization magnitude]”, which in this analysis is “./inversion YAP1 ./YAP1_haplotype_count_rep1.csv index_matrix.csv ./YAP1_Reference_AA_sequence.txt 100”. This will generate the epistasis output in  [Target_protein]_epistasis_rep[replicate_num].txt within ./epistasis_inference/ directory. In this example, the output file name is YAP1_epistasis_rep1.txt.

“g++ -std=c++11 -lgslcblas -lgsl -I ./eigen-3.4.0/ inversion_joint.cpp -o inversion_joint”: compiling the inversion.cpp to estimate the fitness for multiple replicates with intermediate files.

"./inversion_joint [Target-protein] [target protein indexing file name] [replicate number] [optimized regularization magnitude]", which in this analysis is "./inversion_joint YAP1 ./index_matrix.csv 2 100" to have the joint epistasis inference. It will output YAP1_epistasis_joint.txt within the same directory.
