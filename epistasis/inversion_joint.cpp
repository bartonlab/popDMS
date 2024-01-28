#include <algorithm>
#include <array>
#include <iterator>
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <cstring>
#include <unordered_map>
#include <sstream>
#include <stdio.h>
#include <string>
#include <vector>
#include <chrono>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/Dense>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
typedef Eigen::Triplet<double> Trip;
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::SparseVector<double> SpVec;
typedef Eigen::MatrixXd dMat;
template <class OutputIt>
OutputIt safe_tokenizer(const std::string & s, char token, OutputIt out){
    std::string::size_type pos = 0, f;  
    while((f = s.find(token, pos)) != std::string::npos){    
        *out++ = s.substr(pos, f - pos);
        pos = f + 1;
    }
    if(pos < s.size())
        *out++ = s.substr(pos);

    return out;
}

struct Variant_Table{
    int replicate;
    int generation;
    std::vector<int> variant;
    int counts;
    double frequency;
};

typedef std::pair<int, int> pairs;
typedef struct Variant_Table Struct; 

void get_lines_num(FILE* file_head, int *line_num){
    char line[128];

    std::string delimiter = ",";

    std::string token;

    while (fgets(line, 128, file_head)){
        *line_num = *line_num + 1;
    }
}

std::vector<std::array<int, 3>> get_lines_info(FILE* file_head){
    char line[128];

    std::vector<std::array<int, 3>> Info_vec;
    std::string delimiter = ",";
    int i = 0;
    int rep_tmp;
    int gen_tmp;
    int cnt_tmp;
    int vec_list[4];

    std::string substr("total_reads");
    fgets(line, 128, file_head);
    std::string token;

    while (fgets(line, 128, file_head)){

        std::string str(line);
        if(str.find(substr)!=std::string::npos){
            //std::cout << line;
            size_t pos = 0;
            i = 0;
            while ((pos = str.find(delimiter)) != std::string::npos) {
                token = str.substr(0, pos);
                if(i != 2){
                    vec_list[i] = stoi(token);
                }
                //std::cout << token <<std::endl;
                i++;
                str.erase(0, pos + delimiter.length());
            }
            vec_list[i] = stoi(str);
            //std::cout << vec_list[0] <<" "<< vec_list[1] <<" "<<vec_list[3] <<" "<<std::endl; 
            Info_vec.push_back({vec_list[0], vec_list[1], vec_list[3]});
        }
    }

    return Info_vec;  //rep gen counts
}

std::vector<std::string> string_split(std::string s, const char delimiter){
    size_t start = 0;
    size_t end   = s.find_first_of(delimiter);
    std::vector<std::string> output;
    
    while (end <= std::string::npos)
    {
        output.emplace_back(s.substr(start, end - start));
        if (end == std::string::npos)
            break;
        start = end + 1;
        end   = s.find_first_of(delimiter, start);
    }
    
    return output;
}

//['*', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
int get_AA_idx(char AA){
    std::map<char, int> AA_idx = {{'*', 0}, {'A', 1}, {'C', 2}, {'D', 3}, {'E', 4}, {'F', 5}, {'G', 6}, 
                                  {'H', 7}, {'I', 8}, {'K', 9}, {'L', 10}, {'M', 11}, {'N', 12}, {'P', 13}, 
                                  {'Q', 14}, {'R', 15}, {'S', 16}, {'T', 17}, {'V', 18}, {'W', 19}, {'Y', 20}};
    int idx = AA_idx[AA];
    return idx;
}


char codon_trans(std::string codon){
 //-------------------------TXX
    if(codon[0] == 'T'){
        if(codon[1] == 'T'){
            if(codon[2] == 'T' || codon[2] == 'C'){
                return 'F';
            }
            else if(codon[2] == 'A' || codon[2] == 'G'){
                return 'L';
            }
        }
        else if(codon[1] == 'A'){
            if(codon[2] == 'T' || codon[2] == 'C'){
                return 'Y';
            }
            else if(codon[2] == 'A' || codon[2] == 'G'){
                return '*';
            }
        }
        else if(codon[1] == 'C'){
            return 'S';
        }
        else if(codon[1] == 'G'){
            if(codon[2] == 'T' || codon[2] == 'C'){
                return 'C';
            }
            else if(codon[2] == 'G'){
                return 'W';
            }
            else if(codon[2] == 'A'){
                return '*';
            }
        }
    }
 //----------------------------------CXX
    else if(codon[0] == 'C'){
        if(codon[1] == 'T'){
            return 'L';
        }
        else if(codon[1] == 'C'){
            return 'P';
        }
        else if(codon[1] == 'A'){
            if(codon[2] == 'T' || codon[2] == 'C'){
                return 'H';
            }
            else if(codon[2] == 'G' || codon[2] == 'A'){
                return 'Q';
            }
        }
        else if(codon[1] == 'G'){
            return 'R';
        }
    }
 //----------------------------------AXX
    else if(codon[0] == 'A'){
        if(codon[1] == 'T'){
            if (codon[2] == 'G'){
                return 'M';
            }
            else{
                return 'I';
            }
        }
        else if(codon[1] == 'C'){
            return 'T';
        }
        else if(codon[1] == 'A'){
            if(codon[2] == 'T' || codon[2] == 'C'){
                return 'N';
            }
            else if(codon[2] == 'G' || codon[2] == 'A'){
                return 'K';
            }      
        }
        else if(codon[1] == 'G'){
            if(codon[2] == 'T' || codon[2] == 'C'){
                return 'S';
            }
            else if(codon[2] == 'G' || codon[2] == 'A'){
                return 'R';
            }     
        }
    }
 //-----------------------------------------GXX
    else if(codon[0] == 'G'){
        if(codon[1] == 'T'){
            return 'V';
        }
        else if(codon[1] == 'C'){
            return 'A';
        }
        else if(codon[1] == 'A'){
            if(codon[2] == 'T' || codon[2] == 'C'){
                return 'D';
            }
            else if(codon[2] == 'G' || codon[2] == 'A'){
                return 'E';
            }      
        }
        else if(codon[1] == 'G'){
            return 'G';
        }    
    }
    return '*';
}


int mut2idx(char sin_AA, int site_idx){
    int vector_idx;
    int AA_idx = get_AA_idx(sin_AA);
    return (site_idx - 1) * 21 + AA_idx;
}


int get_site(std::string item) {
   std::stringstream variant;
   variant << item; //convert the string s into stringstream
   std::string temp_str;
   int temp_int;
   while(!variant.eof()) {
      variant >> temp_str; //take words into temp_str one by one
      if(std::stringstream(temp_str) >> temp_int) { //try to convert string to int
         //return temp_int;
      }
      //temp_str = ""; //clear temp string
   }
   return temp_int;
}


void printSelectionCoefficients(FILE *output, const std::vector<double> &s) {
    
    for (int i=0;i<s.size();i++) fprintf(output,"%.6e\n",s[i]);
    fflush(output);

}


int main(int argc, char* argv[]) {
    const char* protein_name = argv[1];
    const char* Index_file = argv[2];
    const char* regular = argv[4];
    const char* replicate = argv[3];
    const char* covariance_matrix_file;
    const char* frequency_change_file;
    double regularization = atof(regular);
    int replicates = atoi(replicate);
    FILE* get_line = fopen(Index_file, "r");
    int matrix_dim = 0;
    get_lines_num(get_line, &matrix_dim);
    std::cout << "Dimension of matrix = " << matrix_dim << std::endl;
    double *dx       = new double[matrix_dim];                           // difference between start and end allele frequencies
    double *totalCov = new double[matrix_dim*matrix_dim];                // accumulated allele covariance matrix
    // double *freq_vect_full = new double[replicates*matrix_dim];
    for (int a=0; a < matrix_dim;a++)            dx[a]       = 0;
    for (int a=0; a < matrix_dim*matrix_dim;a++) totalCov[a] = 0;

    double *dx_each_rep       = new double[matrix_dim];                           // difference between start and end allele frequencies
    double *totalCov_each_rep = new double[matrix_dim*matrix_dim];                // accumulated allele covariance matrix
    char line[128];
    double freq;
    int row_num, col_num;
    std::string delimiter = ",";
    std::string token;
    std::string file_name;

    for(int rep = 0; rep < replicates; rep++){
        // for (int a=0; a < matrix_dim; a++)            dx_each_rep[a]       = 0;
        // for (int a=0; a < matrix_dim*matrix_dim; a++) totalCov_each_rep[a] = 0;
        file_name = std::string(protein_name) + "_freq_change_rep" + std::to_string(rep+1) + ".txt";
        frequency_change_file = file_name.c_str();
        get_line = fopen(frequency_change_file, "r");
        while(fgets(line, 128, get_line)){
            size_t pos = 0;
            int col_idx = 0;
            std::string s = line;
            while ((pos = s.find(delimiter)) != std::string::npos) {
                token = s.substr(0, pos);
                if(col_idx == 0){
                    row_num = std::stoi(token);
                    //std::cout << token<< "\t";
                }
                s.erase(0, pos + delimiter.length());
                col_idx++;
            }
            freq = std::stod(s);
            dx[row_num] += freq;
        }
        // for (int a=0; a < matrix_dim; a++)   dx[a] += dx_each_rep[a];

        file_name = std::string(protein_name) + "_cov_matrix_rep" + std::to_string(rep+1) + ".txt";
        covariance_matrix_file = file_name.c_str();
        get_line = fopen(covariance_matrix_file, "r");
        while(fgets(line, 128, get_line)){
            // std::cout<<line;
            size_t pos = 0;
            int col_idx = 0;
            std::string s = line;
            while (col_idx<3) {
                pos = s.find(delimiter);
                token = s.substr(0, pos);
                // std::cout<<token<<std::endl;
                if(col_idx == 0){
                    row_num = std::stoi(token);
                }
                if(col_idx == 1){
                    col_num = std::stoi(token);
                }
                if(col_idx == 2){
                    // std::cout<<token<<std::endl;
                    freq = std::stod(token);
                }                 
                col_idx++;
                s.erase(0, pos + delimiter.length());
            }
            // std::cout<<"---"<<freq<<" "<<row_num<<" "<<col_num<<std::endl;
            if(row_num == col_num) totalCov[col_num*matrix_dim+row_num] += freq;
            else{
                
                totalCov[col_num*matrix_dim+row_num] += freq;
                totalCov[row_num*matrix_dim+col_num] += freq;
            }
        }

    }

    // for(int i =0;i<matrix_dim;i++){
    //     if (dx[i]!=0) std::cout<<i<<" "<<dx[i]<<std::endl;
    //     for(int j =i;j<matrix_dim;j++){
    //         if (totalCov[i*matrix_dim+j]!=0)
    //             std::cout<<i<<" "<<j<<" "<<totalCov[i*matrix_dim+j]<<std::endl;
    //     }
    // }
    for(int i =0;i<matrix_dim;i++){
        totalCov[i*matrix_dim+i] += regularization;
        // std::cout<<dx[i]<<std::endl;
    }

    int status;
    std::vector<double> sMAP(matrix_dim,0);

    std::cout<<"prepare vectors"<<std::endl;
    gsl_matrix_view _cov = gsl_matrix_view_array(totalCov, matrix_dim, matrix_dim);   // gsl covariance + Gaussian regularization
    gsl_vector_view  _dx = gsl_vector_view_array(dx, matrix_dim);            // gsl dx vector
    std::cout<<"prepare space"<<std::endl;
    gsl_vector    *_sMAP = gsl_vector_alloc(matrix_dim);                     // maximum a posteriori selection coefficients for each allele
    gsl_permutation  *_p = gsl_permutation_alloc(matrix_dim);
    std::cout<<"decomposition"<<std::endl;
    gsl_linalg_LU_decomp(&_cov.matrix, _p, &status);
    std::cout<<"solve"<<std::endl;
    gsl_linalg_LU_solve(&_cov.matrix, _p, &_dx.vector, _sMAP);
    
    for (int a=0;a<matrix_dim;a++) {
        sMAP[a] = gsl_vector_get(_sMAP, a);
    }

    FILE *fp;  
    // fp = fopen(std::string(protein_name)+"_epistasis_rep"+std::to_string(Info_vec[0][0])+".txt", "w");
    std::string output_file_name = std::string(protein_name)+"_epistasis_joint.txt";
    const char *output_file_name_char = output_file_name.c_str();
    fp = fopen(output_file_name_char, "w");
    std::cout<<"write file"<<std::endl;
    for (int i=0;i<sMAP.size();i++) {
        fprintf(fp,"%.6e\n",sMAP[i]);
    }
    return 0;
}


