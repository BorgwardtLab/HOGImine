#include <getopt.h>
#include "hogimine.hpp"
#include "hogimine_wy.hpp"
#include "data.hpp"
#include "cmhhandler.hpp"
#include "timekeeping.h"

using namespace std;




int main(int argc, char** argv){
    double tic = second();
    // X Y snps cov map edge
    string data_file;
    string data_file2;
    string snp_file;
    string covar_file;
    string labels_file;
    string map_file;
    string edge_file;
    string out_file;

    bool flag_input_file = false;
    bool flag_input_file2 = false;
    bool flag_snp_file = false;
    bool flag_covar_file = false;
    bool flag_labels_file = false;
    bool flag_map_file = false;
    bool flag_edge_file = false;
    bool flag_fwer = false;
    bool flag_out = false;
    bool print_pv = false;
    double target_fwer = 0.05;
    int n_threads = 1;
    int n_perm = 0;
    int max_interval_dim = 0;

    char opt;
    while( (opt = getopt(argc, argv, "i:h:s:c:l:m:e:f:n:p:d:o:v")) != -1){
        switch(opt){
            case 'i':
                data_file = string(optarg);
                flag_input_file = true;
                break;
            case 'h':
                data_file2 = string(optarg);
                flag_input_file2 = true;
                break;
            case 's':
                snp_file = string(optarg);
                flag_snp_file = true;
                break;
            case 'c':
                covar_file = string(optarg);
                flag_covar_file = true;
                break;
            case 'l':
                labels_file = string(optarg);
                flag_labels_file = true;
                break;
            case 'm':
                map_file = string(optarg);
                flag_map_file = true;
                break;
            case 'e':
                edge_file = string(optarg);
                flag_edge_file = true;
                break;
            case 'f':
                target_fwer = atof(optarg);
                flag_fwer = true;
                break;
            case 'n':
                n_threads = atof(optarg);
                break;
            case 'p':
                n_perm = atof(optarg);
                break;
            case 'd':
                max_interval_dim = atof(optarg);
                break;
            case 'o':
                out_file = string(optarg);
                flag_out = true;
                break;
            case 'v':
                print_pv = true;
                break;
            default:
                cout << "Incorrect argument " << opt << endl;
                exit(1);
        }
    }
    if(!flag_input_file || !flag_labels_file || !flag_edge_file || !flag_snp_file || !flag_map_file || !flag_covar_file || !flag_out){
        cout << "Missing input " << endl;
        exit(1);
    }
    if(!flag_input_file2){

    }
    else{

    }

    Data data = Data(data_file, data_file2, labels_file, snp_file, covar_file, map_file, edge_file);
    if(n_perm > 0){
        HOGImineWY algo = HOGImineWY(&data, target_fwer, max_interval_dim, n_perm);
        cout << "Init: " << second() - tic << " seconds " << endl;
        algo.process_edges();
        algo.write_significant(out_file, print_pv);
        cout << "Total: " << second() - tic << " seconds " << endl;
        cout << "N. processed patterns: " << algo.cnt_process << endl;
        cout << "Significance threshold: " << algo.cmh.corr_threshold() << " n. significant patterns: " << algo.cnt_signif << endl;
    }
    else{
        HOGImine algo = HOGImine(&data, target_fwer, max_interval_dim);
        cout << "Init: " << second() - tic << " seconds " << endl;
        algo.process_edges();
        algo.write_significant(out_file, print_pv);
        cout << "Total: " << second() - tic << " seconds " << endl;
        cout << "N. processed patterns: " << algo.cnt_process << endl;
        cout << "Significance threshold: " << algo.cmh.corr_threshold() << " n. significant patterns: " << algo.cnt_signif << endl;
    }

}
