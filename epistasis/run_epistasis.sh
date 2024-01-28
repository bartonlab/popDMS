python3 variant2count.py ./YAP1.json

g++ -std=c++11 -lgslcblas -lgsl -I ../eigen-3.4.0/ get_freq.cpp -o get_freq

g++ -std=c++11 -lgslcblas -lgsl -I ../eigen-3.4.0/ inversion.cpp -o inversion

g++ -std=c++11 -lgslcblas -lgsl -I ../eigen-3.4.0/ inversion_joint.cpp -o inversion_joint

./get_freq YAP1 ./YAP1_haplotype_count_rep1.csv index_matrix.csv ./YAP1_Reference_AA_sequence.txt

./inversion YAP1 ./YAP1_haplotype_count_rep1.csv index_matrix.csv ./YAP1_Reference_AA_sequence.txt 100

./get_freq YAP1 ./YAP1_haplotype_count_rep2.csv index_matrix.csv ./YAP1_Reference_AA_sequence.txt

./inversion YAP1 ./YAP1_haplotype_count_rep2.csv index_matrix.csv ./YAP1_Reference_AA_sequence.txt 100

./inversion_joint YAP1 ./index_matrix.csv 2 100