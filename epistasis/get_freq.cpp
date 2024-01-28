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

typedef Eigen::Triplet<double> Trip;
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::SparseVector<double> SpVec;
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

std::vector<std::array<int, 3>> get_lines_info(FILE* file_head, int *line_num){
    char line[128];

    std::vector<std::array<int, 3>> Info_vec;
    std::string delimiter = ",";
    int i = 0;
    int rep_tmp;
    int gen_tmp;
    int cnt_tmp;
    int vec_list[4];

    std::string substr("total_reads");
    // std::cout << "Start!-\n";
    fgets(line, 128, file_head);
    // std::cout << "Start!\n";
    std::string token;

    while (fgets(line, 128, file_head)){

        *line_num = *line_num + 1;
        std::string str(line);
        if(str.find(substr)!=std::string::npos){
            // std::cout << line;
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
void display(const std::set<pairs>& s)
{
    bool found = false;
  
    // range-based for loop
    for (auto const &x : s) {
        found = true;
        std::cout << "(" << x.first << ", "
             << x.second << ")"
             << " ";
    }
  
    if (not found) {
        std::cout << "No valid pair\n";
    }
}

Struct extract_idx(std::string line, std::string delimiter, int rep, int gen, 
                   std::vector<int> var, int cnt, std::vector<std::array<int, 3>> Info_vec, 
                   std::set<int> &idx_list, std::vector<int> &idx_raw, 
                   std::set<pairs> &double_idx_set){

    size_t pos = 0;
    Struct var_tmp;
    std::string s = line;
    std::string token;
    std::string sin_var;
    std::string tt_rd = {"total_reads"};
    char sin_AA;
    int col_idx = 0;
    int site_idx;
    int mut_idx;
    std::vector<int> all_variant;
    //std::cout << Info_vec[0].size() << std::endl;

    while ((pos = s.find(delimiter)) != std::string::npos) {
        token = s.substr(0, pos);

        if(col_idx == 0){
            rep = std::stoi(token);
            //std::cout << rep << std::endl;
        }
        if(col_idx == 1){
            gen = std::stoi(token);
        }  
        if(col_idx == 2){
            if(token != tt_rd){
                std::vector<std::string> container = string_split(token, ' ');
                std::vector<int> mutidx_vec;
                for(int ii=0; ii < container.size(); ii++){
                    sin_var = container.at(ii).substr(container.at(ii).length() - 3);
                    sin_AA = codon_trans(sin_var);
                    site_idx = get_site(container.at(ii));
                    mut_idx = mut2idx(sin_AA, site_idx);
                    //std::cout << container[0] << '\n';
                    mutidx_vec.emplace_back(mut_idx);
                    idx_list.insert(mut_idx);
                    //std::cout << sin_var << site_idx << '\n';
                    //std::cout << get_AA_idx(sin_AA) << ' ' << sin_AA << ' ' << site_idx << ' ' <<mut_idx << '\n';
                }
                var = mutidx_vec;
                all_variant = idx_raw;

                for(int i = 0; i < var.size(); i++){
                    all_variant[var[i]/21] = var[i];
                }

                for(int i =0;i<all_variant.size();i++){
                    double_idx_set.insert(std::make_pair(all_variant[i], all_variant[i]));
                    for(int j=i+1;j<all_variant.size();j++){
                        double_idx_set.insert(std::make_pair(all_variant[i], all_variant[j]));
                    }
                }
                
            } 
            else{
                var = {-1};
            }
        }
        if(col_idx == 3){
            cnt = std::stoi(token);
        }
        col_idx++;
        s.erase(0, pos + delimiter.length());
    }
    var_tmp.replicate = rep;
    var_tmp.variant = var;
    var_tmp.generation = gen;
    var_tmp.counts = std::stoi(s);
    //std::cout << Info_vec.size() << std::endl;
    //std::cout << Info_vec[0][2] << std::endl;
    for(int vec_idx = 0; vec_idx < Info_vec.size(); vec_idx++){
        if(Info_vec[vec_idx][0] == var_tmp.replicate && Info_vec[vec_idx][1] == var_tmp.generation){
            var_tmp.frequency = (double)var_tmp.counts/(double)Info_vec[vec_idx][2];
        }
    }
    //std::cout << var_tmp.replicate << ' ' << var_tmp.generation << ' ' << var_tmp.counts << ' ' << var_tmp.frequency << std::endl;
    return var_tmp;
}


std::vector<Trip> initial_sin(int raw_seq_idx[], int seq_length){
    std::vector<Trip> tripletList;
    //std::map<std::tuple<int, int>, int> sin2dou;
    for(int sin_idx = 0; sin_idx < seq_length; sin_idx++){  // single site single frequency
        //std::cout << raw_seq_idx[sin_idx]<<std::endl;
        tripletList.push_back(Trip(0, raw_seq_idx[sin_idx], 1));
    }

    return tripletList;
}


/*
    std::vector<Trip> tripletList;
    tripletList.reserve(6);
    int i = 0;
    int j = 0;
    for(; i<10 && j<10; i++, j++){
        tripletList.push_back(Trip(i, j, 0.1*i));
    }


    SpMat A_mat(10, 10);
    A_mat.setFromTriplets(tripletList.begin(), tripletList.end());
*/


void cov_init(int seq_len, std::vector<Trip> &cov_matr_tmp, std::vector<int> idx_raw, SpMat sin_idx, SpMat dou_idx){
    std::cout<<"Raw sequence index: ";
    for(int i = 0; i<idx_raw.size();i++){
        std::cout<<idx_raw[i]<<" ";
    }
    std::cout<<std::endl;
    for(int site1 = 0; site1<seq_len; site1++){
        for(int site2 = site1; site2<seq_len; site2++){
            if(site1 != site2){ //AB
                cov_matr_tmp.emplace_back(Trip(sin_idx.coeff(0,idx_raw[site1]), sin_idx.coeff(0,idx_raw[site2]), 1));
                cov_matr_tmp.emplace_back(Trip(dou_idx.coeff(idx_raw[site1], idx_raw[site2]), dou_idx.coeff(idx_raw[site1], idx_raw[site2]), 1));
                //std::cout<<idx_raw[site1]<<" "<<idx_raw[site2]<<" "<<sin_idx.coeff(0, idx_raw[site1])<<std::endl;
            }
            else{ //AA
                cov_matr_tmp.emplace_back(Trip(sin_idx.coeff(0,idx_raw[site1]), sin_idx.coeff(0,idx_raw[site1]), 1));
            }
            for(int site3 = site2; site3<seq_len; site3++){
                if(site1 != site2){
                    if(site2 != site3){ //ABC
                        cov_matr_tmp.emplace_back(Trip(sin_idx.coeff(0,idx_raw[site1]), dou_idx.coeff(idx_raw[site2], idx_raw[site3]), 1));
                        cov_matr_tmp.emplace_back(Trip(sin_idx.coeff(0,idx_raw[site2]), dou_idx.coeff(idx_raw[site1], idx_raw[site3]), 1));
                        cov_matr_tmp.emplace_back(Trip(sin_idx.coeff(0,idx_raw[site3]), dou_idx.coeff(idx_raw[site1], idx_raw[site2]), 1));
                        //std::cout << idx_dou[std::make_tuple(idx_raw[site1], idx_raw[site2])]<<","<<dou_idx.coeff(idx_raw[site1], idx_raw[site2])<<","<< idx_sin[idx_raw[site3]]<<","<< sin_idx.coeff(0, idx_raw[site3]) <<std::endl;
                    }
                    else{ //ABB
                        cov_matr_tmp.emplace_back(Trip(sin_idx.coeff(0,idx_raw[site2]), dou_idx.coeff(idx_raw[site1], idx_raw[site3]), 1));
                    }
                }
                else{ //AAB
                    if(site2 != site3){
                        cov_matr_tmp.emplace_back(Trip(sin_idx.coeff(0,idx_raw[site1]), dou_idx.coeff(idx_raw[site2], idx_raw[site3]), 1));
                    }
                    else{} //AAA
                }
                for(int site4 = site3; site4<seq_len; site4++){
                    if(site1 == site2){
                        if(site2 == site3){} //AAAX
                        else{
                            if(site3 == site4){} //AABB
                            else{ //AABC
                                //cov_matr_tmp.emplace_back(Trip(idx_sin[idx_raw[site1]], idx_dou[std::make_tuple(idx_raw[site3], idx_raw[site4])], 1));
                                cov_matr_tmp.emplace_back(Trip(dou_idx.coeff(idx_raw[site1], idx_raw[site3]), dou_idx.coeff(idx_raw[site2], idx_raw[site4]), 1));
                            }
                        }
                    }
                    else{
                        if(site2 == site3){
                            if(site3 == site4){ //ABBB
                                //cov_matr_tmp.emplace_back(Trip(idx_sin[idx_raw[site2]], idx_dou[std::make_tuple(idx_raw[site1], idx_raw[site3])], 1));
                            }
                            else{ // ABBC
                                //cov_matr_tmp.emplace_back(Trip(idx_sin[idx_raw[site2]], idx_dou[std::make_tuple(idx_raw[site1], idx_raw[site4])], 1));
                                cov_matr_tmp.emplace_back(Trip(dou_idx.coeff(idx_raw[site1], idx_raw[site2]), dou_idx.coeff(idx_raw[site3], idx_raw[site4]), 1));
                            }
                        }
                        else{ 
                            if(site3 == site4){ //ABCC
                                //cov_matr_tmp.emplace_back(Trip(idx_sin[idx_raw[site3]], idx_dou[std::make_tuple(idx_raw[site1], idx_raw[site2])], 1));
                                cov_matr_tmp.emplace_back(Trip(dou_idx.coeff(idx_raw[site1], idx_raw[site3]), dou_idx.coeff(idx_raw[site2], idx_raw[site4]), 1));
                            }
                            else{ //ABCD
                                cov_matr_tmp.emplace_back(Trip(dou_idx.coeff(idx_raw[site1], idx_raw[site2]), dou_idx.coeff(idx_raw[site3], idx_raw[site4]), 1));
                                cov_matr_tmp.emplace_back(Trip(dou_idx.coeff(idx_raw[site1], idx_raw[site3]), dou_idx.coeff(idx_raw[site2], idx_raw[site4]), 1));
                                cov_matr_tmp.emplace_back(Trip(dou_idx.coeff(idx_raw[site1], idx_raw[site4]), dou_idx.coeff(idx_raw[site2], idx_raw[site3]), 1));
                            }
                        }
                    }
                }
            }
        }
    }
}

void cov_trip(std::tuple<int, int> rep_gen, int site1, int site2, int site3, int var1, int var2, int var3, int wt1, int wt2, int wt3, 
              std::map<int, int> idx_sin, std::map<std::tuple<int, int>, int> idx_dou, 
              std::vector<Trip> &cov_matx, double freq, SpMat sin_idx, SpMat dou_idx){

    int idx1, idx2, idx1_wt, idx2_wt;
    //std::cout<<wt1<<" "<<wt2<<" "<<wt3<<std::endl;
    //std::cout<<var1<<" "<<var2<<" "<<var3<<std::endl;
    if(var1 == var2){ //AAB
        idx1 = sin_idx.coeff(0, var2);
        idx1_wt = sin_idx.coeff(0, wt2);
        idx2 = dou_idx.coeff(var1, var3);
        idx2_wt = dou_idx.coeff(wt1, wt3);
        cov_matx.emplace_back(Trip(idx1, idx2, freq));
        cov_matx.emplace_back(Trip(idx1_wt, idx2_wt, -freq));
    }
    else if(var2 == var3){ // ABB
        idx1 = sin_idx.coeff(0, var2);
        idx1_wt = sin_idx.coeff(0, wt2);
        idx2 = dou_idx.coeff(var1, var3);
        idx2_wt = dou_idx.coeff(wt1, wt3);
        cov_matx.emplace_back(Trip(idx1, idx2, freq));
        cov_matx.emplace_back(Trip(idx1_wt, idx2_wt, -freq));
    }
    else{ //ABC
        idx1 = sin_idx.coeff(0, var1);
        idx1_wt = sin_idx.coeff(0, wt1);
        idx2 = dou_idx.coeff(var2, var3);
        idx2_wt = dou_idx.coeff(wt2, wt3);
        //idx2 = dou_idx.coeff(var2, var3);
        cov_matx.emplace_back(Trip(idx1, idx2, freq));
        cov_matx.emplace_back(Trip(idx1_wt, idx2_wt, -freq));

        idx1 = sin_idx.coeff(0, var2);
        //idx2 = dou_idx.coeff(var1, var3);
        idx1_wt = sin_idx.coeff(0, wt2);        
        //idx2_wt = dou_idx.coeff(wt1, wt3);
        idx2 = dou_idx.coeff(var1, var3);
        idx2_wt = dou_idx.coeff(wt1, wt3);
        cov_matx.emplace_back(Trip(idx1, idx2, freq));
        cov_matx.emplace_back(Trip(idx1_wt, idx2_wt, -freq));

        idx1 = sin_idx.coeff(0, var3);
        //idx2 = dou_idx.coeff(var1, var2);
        idx1_wt = sin_idx.coeff(0, wt3);
        //idx2_wt = dou_idx.coeff(wt1, wt2);
        idx2 = dou_idx.coeff(var1, var2);
        idx2_wt = dou_idx.coeff(wt1, wt2);
        cov_matx.emplace_back(Trip(idx1, idx2, freq));
        cov_matx.emplace_back(Trip(idx1_wt, idx2_wt, -freq));
    }
}

void cov_quar(std::tuple<int, int> rep_gen, int site1, int site2, int site3, int site4, int var1, int var2, int var3, int var4, int wt1, int wt2, int wt3, int wt4,
              std::map<int, int> idx_sin, std::map<std::tuple<int, int>, int> idx_dou, 
              std::vector<Trip> &cov_matx, double freq, SpMat sin_idx, SpMat dou_idx, int variant_num){

    int idx1, idx2, idx1_wt, idx2_wt;

    //std::cout<< 

    //12 34
    idx1 = dou_idx.coeff(var1, var2);
    idx2 = dou_idx.coeff(var3, var4);
    if(idx1 == idx2){}
    else if(idx1 < variant_num || idx2 < variant_num){}
    else{
        idx1_wt = dou_idx.coeff(wt1, wt2);
        idx2_wt = dou_idx.coeff(wt3, wt4);
        if(idx1<idx2){
            cov_matx.emplace_back(Trip(idx1, idx2, freq));
            cov_matx.emplace_back(Trip(idx1_wt, idx2_wt, -freq));  
        }
        else{
            cov_matx.emplace_back(Trip(idx2, idx1, freq));
            cov_matx.emplace_back(Trip(idx2_wt, idx1_wt, -freq));              
        }
    }

    //13 24
    idx1 = dou_idx.coeff(var1, var3);
    idx2 = dou_idx.coeff(var2, var4);
    if(idx1 == idx2){}
    else if(idx1 < variant_num || idx2 < variant_num){}
    else{
        idx1_wt = dou_idx.coeff(wt1, wt3);
        idx2_wt = dou_idx.coeff(wt2, wt4);  
        if(idx1<idx2){
            cov_matx.emplace_back(Trip(idx1, idx2, freq));
            cov_matx.emplace_back(Trip(idx1_wt, idx2_wt, -freq));  
        }
        else{
            cov_matx.emplace_back(Trip(idx2, idx1, freq));
            cov_matx.emplace_back(Trip(idx2_wt, idx1_wt, -freq));              
        }        
    }

    //14 23
    idx1 = dou_idx.coeff(var1, var4);
    idx2 = dou_idx.coeff(var2, var3);
    if(idx1 == idx2){}
    else if(idx1 < variant_num || idx2 < variant_num){}
    else{
        idx1_wt = dou_idx.coeff(wt1, wt4);
        idx2_wt = dou_idx.coeff(wt2, wt3);  
        if(idx1<idx2){
            cov_matx.emplace_back(Trip(idx1, idx2, freq));
            cov_matx.emplace_back(Trip(idx1_wt, idx2_wt, -freq));  
        }
        else{
            cov_matx.emplace_back(Trip(idx2, idx1, freq));
            cov_matx.emplace_back(Trip(idx2_wt, idx1_wt, -freq));              
        }      
    }
}


void MergeSortedVectors(std::vector<int>& v1, std::vector<int>& v2, std::vector<int>& v3) {

    // v3 is the output vector
    // it will store the merged vector obtained by merging v1 and v2

    int i, j, n, m;
    i = j = 0;
    n = v1.size();
    m = v2.size();


    // traverse each elemenst of v1 and v2
    while (i < n && j < m) {

        // comparing v1[i] and v2[j]
        // if v1[i] is smaller than v2[j]
        // push v1[i] to v3 and increment i
        // if v[i] is less than v2[j]
        // push v2[j] to v3 and increment j
        if (v1[i] <= v2[j]) {
            v3.push_back(v1[i]);
            ++i;
        }
        else {
            v3.push_back(v2[j]);
            ++j;
        }
    }

    // push the elements left in v1 to v3
    while (i < n) {
        v3.push_back(v1[i]);
        ++i;
    }

    // push the elements left in v2 to v3
    while (j < m) {
        v3.push_back(v2[j]);
        ++j;
    }
}

void build_cov( const int &seq_len, int dimension, int regularization,
                std::vector<Variant_Table> variant_record, 
                std::vector<std::array<int, 3>> Info_vec,  //Info_vec: rep gen counts
                std::vector<int> idx_raw, //raw alleles idx
                std::set<int> idx_list, //all alleles idx
                std::unordered_map<int, std::unordered_map<int, int>> double_idx,
                std::map<std::tuple<int, int>, SpMat> &frq_vect_full,
                std::map<std::tuple<int, int>, SpMat> &cov_matx_full,
                const char* protein_name){

    int rep;
    int gen;
    int inter_len = int(seq_len*(seq_len-1)/2) + seq_len;
    //int dimension = double_idx.size();
    //std::cout<<"Dimension= "<<dimension<<std::endl;
    std::vector<Trip> cov_matx_tmp;
    std::vector<Trip> fre_vect_tmp;
    std::unordered_map<int, std::unordered_map<int, std::unordered_map< int, std::unordered_map<int, double>>>> matrix_element;
    std::unordered_map<int, std::unordered_map<int, std::unordered_map< int, std::unordered_map<int, double>>>> frequency_element;
    std::vector<int> all_site;
    std::vector<int> tmp_idx_sin;
    std::vector<int> tmp_idx_dou;
    int raw_seq[seq_len];

    for(int i = 0; i < seq_len; i++){
        all_site.emplace_back(i);
    }

    SpMat init_cov_mat(dimension, dimension);
    SpMat init_fre_vec(1, dimension);

    std::cout<<"Raw sequence index: ";
    for(int i = 0; i<idx_raw.size();i++){
        std::cout<<idx_raw[i]<<" ";
        raw_seq[i] = idx_raw[i];
    }
    std::cout<<std::endl;


    std::cout<<"All sequence index: ";
    for(auto element:idx_list){
        std::cout<<element<< " ";
        //tmp_idx_sin.emplace_back(element);
    }
    std::cout<<std::endl;
/*
    for(auto const &ent1 : double_idx) {
        // ent1.first is the first key
        for(auto const &ent2 : ent1.second) {
            // ent2.first is the second key, ent2.second is the data
            //std::cout<< ent1.first <<  " "<<  ent2.first << " "<<ent2.second <<std::endl;
            tmp_idx_dou.emplace_back(ent2.second);
        }
    }
*/
    //std::cout<<"==="<<std::endl;
    int tmp1, tmp2;

    for(int i = 0; i < Info_vec.size(); i++){
        rep = Info_vec[i][0];
        gen = Info_vec[i][1];

        cov_matx_full[std::make_tuple(rep, gen)] = init_cov_mat;
        frq_vect_full[std::make_tuple(rep, gen)] = init_fre_vec;
/*
        for(int j =0; j<tmp_idx_sin.size(); j++){
            tmp1 = double_idx[tmp_idx_sin[j]][tmp_idx_sin[j]];
            matrix_element[rep][gen][tmp1][tmp1] = 0;
            frequency_element[rep][gen][0][tmp1] = 0;
            for(int k =j+1; k<tmp_idx_sin.size(); k++){
                tmp2 = double_idx[tmp_idx_sin[j]][tmp_idx_sin[k]];
                matrix_element[rep][gen][tmp1][tmp2] = 0;
                if(tmp1>tmp2){
                    std::cout<<tmp1<<" "<<tmp2<<std::endl;
                }
            }
        }*/
    }




    int var_num;//, wt_num;
    double freq;
    //int site1, site2, site3, site4;
    //int mut_var1, wt_var1, mut_var2, wt_var2, mut_var3, wt_var3, mut_var4, wt_var4;
    //int var1, var2, var3, var4, wt1, wt2, wt3, wt4;

    int line_num = variant_record.size();
    //line_num = 1000;

    int generations = Info_vec[Info_vec.size()-1][1] + 1;

    float mut_haplotype_freq[generations];

    std::cout << "# of lines = " << line_num << std::endl;
    
    for(int line = 0; line < line_num; line++){

        if(line%10 ==0){
            std::cout << line << std::endl;
        }
        
        rep = variant_record[line].replicate;
        gen = variant_record[line].generation;

        if(variant_record[line].variant[0] == -1){}
        else{
            var_num = variant_record[line].variant.size();
            freq = variant_record[line].frequency;
            mut_haplotype_freq[gen] += freq;
            // std::cout<<freq<<std::endl; 
            //std::vector<int> all_variant = idx_raw;
            //std::vector<int> variant_idx;
            std::vector<int> all_variant_sd;

            int var_seq[seq_len];
            //int var_seq_sd[inter_len];
            memcpy(var_seq, raw_seq, sizeof(raw_seq));

            for(int i = 0; i < var_num; i++){
                var_seq[variant_record[line].variant[i]/21] = variant_record[line].variant[i];
                // std::cout<<variant_record[line].variant[i]<<std::endl; 
            }  
            //std::cout<<"==="<<std::endl; 
            for(int i = 0; i<seq_len; i++){
                // std::cout<<var_seq[i]<<"\t"; 
                //var_seq_sd[i*seq_len+i]=double_idx[var_seq[i]][var_seq[i]];
                for(int j = i; j<seq_len; j++){
                    all_variant_sd.emplace_back(double_idx[var_seq[i]][var_seq[j]]);
                    //std::cout<<var_seq[i]<<" "<<var_seq[j]<<" "<<double_idx[var_seq[i]][var_seq[j]]<<std::endl;
                }
            }  

            for(int i = 0; i<all_variant_sd.size(); i++){  
                matrix_element[rep][gen][all_variant_sd[i]][all_variant_sd[i]] += freq;
                frequency_element[rep][gen][0][all_variant_sd[i]] += freq;
                for(int j = i+1; j<all_variant_sd.size(); j++){
                    matrix_element[rep][gen][all_variant_sd[i]][all_variant_sd[j]] += freq;
                    // std::cout << tmp1<<","<<tmp2 << "|\t";
                }
            }
            // exit(0);
        } 
    } 

    for(int wt_gen = 0; wt_gen < generations; wt_gen++){
        std::vector<int> all_variant_sd;
        int var_seq[seq_len];
        memcpy(var_seq, raw_seq, sizeof(raw_seq));
        for(int i = 0; i<seq_len; i++){
            // std::cout<<var_seq[i]<<"\t"; 
            //var_seq_sd[i*seq_len+i]=double_idx[var_seq[i]][var_seq[i]];
            for(int j = i; j<seq_len; j++){
                all_variant_sd.emplace_back(double_idx[var_seq[i]][var_seq[j]]);
                // std::cout<<var_seq[i]<<" "<<var_seq[j]<<" "<<double_idx[var_seq[i]][var_seq[j]]<<std::endl;
            }
        } 
        for(int i = 0; i < all_variant_sd.size(); i++){  
            // std::cout << all_variant_sd[i] << "|\t";
            matrix_element[rep][wt_gen][all_variant_sd[i]][all_variant_sd[i]] += 1 - mut_haplotype_freq[wt_gen];
            // std::cout<<frequency_element[rep][wt_gen][0][all_variant_sd[i]]<<std::endl;
            frequency_element[rep][wt_gen][0][all_variant_sd[i]] += 1 - mut_haplotype_freq[wt_gen];
            // std::cout<<1 - mut_haplotype_freq[wt_gen]<<" "<<frequency_element[rep][wt_gen][0][all_variant_sd[i]] <<std::endl;
            for(int j = i+1; j<all_variant_sd.size(); j++){
                matrix_element[rep][wt_gen][all_variant_sd[i]][all_variant_sd[j]] += 1 - mut_haplotype_freq[wt_gen];
                // std::cout << tmp1<<","<<tmp2 << "|\t";
            }
        }
    }

    // std::cout<<mut_haplotype_freq[0]<<std::endl;
    // std::cout<<mut_haplotype_freq[1]<<std::endl;
    // std::cout<<mut_haplotype_freq[2]<<std::endl;
    // std::cout<<mut_haplotype_freq[3]<<std::endl;


    std::ofstream outputFile;
    std::string filename;

    for(int i = 0; i < Info_vec.size(); i++){

        rep = Info_vec[i][0];
        gen = Info_vec[i][1];
        filename = std::string(protein_name)+"_multiple_allele_"+std::to_string(rep)+"_"+std::to_string(gen)+".csv";
        outputFile.open(filename);
        //std::cout<<rep<<","<<gen<<"----"<<std::endl;
        for(auto const &ent1 : matrix_element[rep][gen]) {
            // ent1.first is the first key
            //cov_matx_tmp.emplace_back(Trip(ent1.first, ent1.first, matrix_element[rep][gen][ent1.first][ent1.first]));
            for(auto const &ent2 : ent1.second) {
                // ent2.first is the second key
                // ent2.second is the data
                outputFile <<ent1.first<<","<< ent2.first << ","<<ent2.second<< std::endl;
                //std::cout << ent1.first<<" "<< ent2.first << ""ent2.second<<std::endl;
                //cov_matx_tmp.emplace_back(Trip(ent1.first, ent2.first, matrix_element[rep][gen][ent1.first][ent2.first]));
                //cov_matx_tmp.emplace_back(Trip(ent2.first, ent1.first, matrix_element[rep][gen][ent1.first][ent2.first]));
            }
        }
        outputFile.close();
    }
    

    for(int i = 0; i < Info_vec.size(); i++){
        rep = Info_vec[i][0];
        gen = Info_vec[i][1];
        filename = std::string(protein_name)+"_freq_"+std::to_string(rep)+"_"+std::to_string(gen)+".csv";
        outputFile.open(filename);

        for(auto const &ent1 : frequency_element[rep][gen]) {
            // ent1.first is the first key
            //cov_matx_tmp.emplace_back(Trip(ent1.first, ent1.first, matrix_element[rep][gen][ent1.first][ent1.first]));
            for(auto const &ent2 : ent1.second) {
                // ent2.first is the second key
                // ent2.second is the data
                outputFile << ent2.first << ","<<ent2.second<< std::endl;
                // std::cout << ent1.first<<" "<< ent2.first << ""ent2.second<<std::endl;
                //cov_matx_tmp.emplace_back(Trip(ent1.first, ent2.first, matrix_element[rep][gen][ent1.first][ent2.first]));
                //cov_matx_tmp.emplace_back(Trip(ent2.first, ent1.first, matrix_element[rep][gen][ent1.first][ent2.first]));
            
            }
        }
        outputFile.close();
    }

}



int main(int argc, char* argv[]) {

    clock_t start, end;
    double cpu_time_used; 
    start = clock();
    // std::cout << "Start!--\n";
    //std::string raw_sequence = "FLIMVD";
    // std::string raw_sequence = "DVPLPAGWEMAKTSSGQRYFLNHIDQTTTWQDPR";

    const char* raw_sequence_file = argv[4];
    std::ifstream t(raw_sequence_file);
    std::stringstream buffer;
    buffer << t.rdbuf();
    std::string raw_sequence = buffer.str();

    std::ofstream outputFile;
    // create a name for the file output
    std::string filename;
    int seq_length = raw_sequence.length();
    int regularization = 1;
    // const char* FileName = "YAP1_genotype_count_rep2.csv";
    const char* protein_name = argv[1];
    const char* FileName = argv[2];
    // int current_rep = argv[2];
    const char *raw_seq = raw_sequence.c_str();
    
    int raw_seq_idx[seq_length];
    for(int tmp_idx = 0; tmp_idx < seq_length; tmp_idx++){
        raw_seq_idx[tmp_idx] = tmp_idx * 21 + get_AA_idx(raw_sequence[tmp_idx]);
    }

// read how many lines in data file & Info_vec = (rep,gen,tot_cnt)
    FILE* get_line = fopen(FileName, "r");
    int line_num = 0;
    // std::cout << "# of lines = " << line_num << std::endl;
    std::vector<std::array<int, 3>> Info_vec = get_lines_info(get_line, &line_num);

    
    //std::cout << Info_vec[0].size() << std::endl;
    //std::cout << "# of lines = " << raw_data_lines << std::endl;

    // construct a data structure array for the file
    //struct Variant_Table record[raw_data_lines-1];
    //std::cout << "size of empty data struct = " << sizeof(record) << std::endl;

// read csv file, parse the strings by ','
    FILE* stream = fopen(FileName, "r");
    char line[128];
    std::string delimiter = ",";

// read the first index line and discard it
    fgets(line, 128, stream);
    std::cout << line << std::endl;

// loop through the file
    int struct_idx = 0;
    std::vector<Variant_Table> variant_record;

// initialize the index vector
    std::set<int> idx_list;
    for(int tmp_i = 0; tmp_i<seq_length; tmp_i++){
        //std::cout << mut2idx(raw_sequence[tmp_i], tmp_i+1) << std::endl;
        idx_list.insert(mut2idx(raw_sequence[tmp_i], tmp_i+1));
    }

    std::vector<int> idx_raw;
    for(auto element:idx_list){
        idx_raw.emplace_back(element);
    }

    std::unordered_map<int, std::unordered_map<int, int>> double_idx;
    std::set<pairs> double_idx_set;
    for(int i =0;i<idx_raw.size();i++){
        double_idx_set.insert(std::make_pair(idx_raw[i], idx_raw[i]));
        for(int j=i+1;j<idx_raw.size();j++){
            double_idx_set.insert(std::make_pair(idx_raw[i], idx_raw[j]));
        }
    }
    // display(double_idx_set);
    std::cout<<std::endl;

    //std::cout << idx_list[0] <<idx_list[1]<< std::endl;

    while (fgets(line, 128, stream)){
        Struct variant_cell;
        int rep = -1;
        int gen = -1;
        int cnt = -1;
        std::vector<int> var;
        variant_cell = extract_idx(line, delimiter, rep, gen, var, cnt, Info_vec, idx_list, idx_raw, double_idx_set);
        variant_record.emplace_back(variant_cell);
        struct_idx++;
    }
    int dimension = double_idx_set.size();

// remove duplicate idx
    //display(double_idx_set);
    std::cout<<"Total index counts: "<<dimension<<std::endl;
    //return 0;
    int idx_num=0;
    for (auto elem : double_idx_set)
    {
        //std::cout << elem.first << ","<<elem.second<<","<<idx_num<<std::endl;
        double_idx[elem.first][elem.second] = idx_num;
        idx_num++;
    }    
    // filename = "index_matrix.csv";
    filename = argv[3];

    outputFile.open(filename);
    for(auto const &ent1 : double_idx) {
        // ent1.first is the first key
        //cov_matx_tmp.emplace_back(Trip(ent1.first, ent1.first, matrix_element[rep][gen][ent1.first][ent1.first]));
        for(auto const &ent2 : ent1.second) {
            // ent2.first is the second key
            // ent2.second is the data
            outputFile <<ent1.first<<","<< ent2.first << ","<<ent2.second<< std::endl;
            //std::cout << ent1.first<<" "<< ent2.first << ""ent2.second<<std::endl;
            //cov_matx_tmp.emplace_back(Trip(ent1.first, ent2.first, matrix_element[rep][gen][ent1.first][ent2.first]));
            //cov_matx_tmp.emplace_back(Trip(ent2.first, ent1.first, matrix_element[rep][gen][ent1.first][ent2.first]));
        }
    }
    outputFile.close();
    //recon_idxlist(idx_list);
    //return 0;
//get idx reference, idx_sin{1:0, 3:1, ...}, idx_dou{(1,5): 14, (1,22):15, ...}
/*    
    for (int k=0; k < dou_idx.outerSize(); ++k){
        for (Eigen::SparseMatrix<double, Eigen::ColMajor>::InnerIterator it(dou_idx,k); it; ++it){
            std::cout << "(" << it.row() << ","<< it.col() << ","<< it.value()<< ")\t"; // row index
        }
    }
*/
    //std::cout << dou_idx.coeff(1, 5) <<","<<idx_dou[std::make_tuple(1, 5)] << std::endl;

// get single/pair frequencies and covariance matrix, tuple = {rep, gen}, Trip = {row, col, val}
    std::vector<Trip> tripletList;

    std::map<std::tuple<int, int>, SpMat> cov_matx;
    std::map<std::tuple<int, int>, SpMat> frq_vect_full;
    
    //int variant_num = triplet_sin.size();
    const int seq_len = idx_raw.size();

    build_cov(seq_len, dimension, regularization, variant_record, Info_vec, idx_raw, idx_list, double_idx, frq_vect_full, cov_matx, protein_name);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    std::cout << cpu_time_used << std::endl;
    return 0;
}

