#ifndef data_h
#define data_h

#include "utils.hpp"

using namespace std;



class Data{
public:
    int n_trans;
    int n_snps;
    vector<Bitset> matrix;
    vector<Bitset> matrix2;
    Bitset labels;
    vector<int> covars;
    vector<string> snp_names;

    int n_genes;
    vector<string> gene_names;
    vector<vector<int>> gene2snps;
    hash_map<string, int> gene2int;
    hash_map<string, int> snp2int;

    int positive;
    int num_covariates;
    vector<int> trans_per_class;
    vector<int> positive_per_class;
    vector<int> covariate_start;

    //vector<vector<int>> adjList;
    vector<vector<int>> edgeList;

    Data(){}
    ~Data(){}
    Data(string data_f, string labels_f, string snps_f, string cov_f, string map_f, string edge_f){
        //load snp information
        pair<int, int> pr = get_dims(data_f);
        n_snps = pr.first; n_trans = pr.second;
        positive = 0;
        num_covariates = 1;

        vector<int> sort_id = vector<int>(n_trans);

        matrix = vector<Bitset>(n_snps);
        matrix2 = vector<Bitset>(n_snps);
        labels = Bitset(n_trans);
        covars = vector<int>(n_trans);
        snp_names = vector<string>(n_snps);

        load_covars(cov_f, sort_id);
        load_labels(labels_f, sort_id);
        load_matrix(data_f, sort_id);
        load_matrix2(data_f, sort_id);

        load_snps(snps_f);

        // make sure 1 is minority class
        if(positive > n_trans/2){
            positive = n_trans - positive;
            for(int i=0; i<n_trans; i++)
                labels.flip(i);
        }

        trans_per_class = vector<int>(num_covariates);
        positive_per_class = vector<int>(num_covariates);
        for(int i=0; i<n_trans; i++){
            trans_per_class[covars[i]]++;
            if(labels.get(i)) positive_per_class[covars[i]]++;
        }


        // load gene information
        load_mapping(map_f);

        // load graph information
        //adjList = vector<vector<int>>(n_genes);
        load_edges(edge_f);
    }



    pair<int, int> get_dims(string namef){
        int n_rows = 0;
        int n_cols = 0;
        double idx;
        ifstream f(namef);
        string line;
        stringstream ss_line;

        while(getline(f, line)){
            if (n_rows == 0) {
              ss_line << line;
              while(ss_line >> idx) n_cols ++;
            }
        n_rows++;
        }
        f.close();

        return make_pair(n_rows, n_cols); // n_snps, n_trans
    }

    void load_covars(string namef, vector<int>& sort_id){
        string line;
        stringstream ss_line;

        int value;
        int row_idx = 0;
        int maxcov = 0;
        vector<pair<int, int>> tmp_covars(n_trans);

        ifstream f(namef);
        while(getline(f, line)){
            ss_line << line;
            ss_line >> value;
            tmp_covars[row_idx] = make_pair(value, row_idx);
            maxcov = max(maxcov, value);
            ss_line.clear();
            row_idx ++;
        }
        f.close();
        if(row_idx != n_trans) cerr << "Error covars" << endl;

        num_covariates = maxcov+1;

        // make covariates contiguous
        covariate_start = vector<int>(num_covariates+1);
        covariate_start[num_covariates] = n_trans;
        int prev_c = -1;
        sort(tmp_covars.begin(), tmp_covars.end());
        for(int i=0; i<n_trans; i++){
            sort_id[tmp_covars[i].second] = i;
            covars[i] = tmp_covars[i].first;
            if(covars[i] != prev_c){
                prev_c = covars[i];
                covariate_start[prev_c] = i;
            }
        }
    }

    void load_matrix(string namef, vector<int>& sort_id){
        string line;
        stringstream ss_line;

        int value;
        int row_idx = 0;
        int col_idx = 0;

        ifstream f(namef);
        while(getline(f, line)){
            matrix[row_idx] = Bitset(n_trans);
            ss_line << line;
            while(ss_line >> value){
                if(value == 2) matrix[row_idx].set( sort_id[col_idx] ); // set matrix[snp][trans] = 1
                col_idx ++;
            }
            ss_line.clear();
            row_idx ++;
            col_idx = 0;
        }
        f.close();
    }

    void load_matrix2(string namef, vector<int>& sort_id){
        string line;
        stringstream ss_line;

        int value;
        int row_idx = 0;
        int col_idx = 0;

        ifstream f(namef);
        while(getline(f, line)){
            matrix2[row_idx] = Bitset(n_trans);
            ss_line << line;
            while(ss_line >> value){
                if(value >= 1) matrix2[row_idx].set( sort_id[col_idx] ); // set matrix[snp][trans] = 1
                col_idx ++;
            }
            ss_line.clear();
            row_idx ++;
            col_idx = 0;
        }
        f.close();
    }

    void load_labels(string namef, vector<int>& sort_id){
        string line;
        stringstream ss_line;

        int value;
        int row_idx = 0;

        ifstream f(namef);
        while(getline(f, line)){
            ss_line << line;
            ss_line >> value;
            if(value) {
                labels.set( sort_id[row_idx] );
                positive++;
            }
            ss_line.clear();
            row_idx ++;
        }
        f.close();
        if(row_idx != n_trans) cerr << "Error labels" << endl;
    }



    void load_snps(string namef){
        string line;
        stringstream ss_line;

        string value;
        int row_idx = 0;

        ifstream f(namef);
        while(getline(f, line)){
            ss_line << line;
            ss_line >> value;
            snp_names[row_idx] = value;
            snp2int[value] = row_idx;
            ss_line.clear();
            row_idx ++;
        }
        f.close();
        if(row_idx != n_snps) cerr << "Error snps" << endl;
    }

    void load_mapping(string namef){
        string line;
        stringstream ss_line;

        int col_idx;
        int row_idx = 0;
        string value;
        string gene_name;
        vector<int> snp_vec;

        ifstream f(namef);
        while(getline(f, line)){
            col_idx = 0;
            ss_line << line;
            while(ss_line >> value){
                if (col_idx == 0){
                    gene_names.push_back(value);
                    gene2int[value] = row_idx;
                }
                else{
                    auto it = snp2int.find(value);
                    if(it != snp2int.end())
                        snp_vec.push_back(snp2int[value]);
                }
                col_idx += 1;
            }
            gene2snps.push_back(snp_vec);
            snp_vec.clear();
            ss_line.clear();
            row_idx++;
        }
        n_genes = row_idx;
        f.close();
    }
    /*
    void load_edges2(string namef){
        string line, value;
        stringstream ss_line;
        vector<unordered_set<int>> tmpList(n_genes);

        ifstream f(namef);
        while(getline(f, line)){
            ss_line << line;
            ss_line >> value;
            int u = -1, v = -1;
            auto it = gene2int.find(value);
            if(it != gene2int.end())
                u = gene2int[value];
            ss_line >> value;
            it = gene2int.find(value);
            if(it != gene2int.end())
                v = gene2int[value];

            if((u>=0) && (v>=0) && (u != v)){ // no self-loops
                tmpList[u].insert(v);
                tmpList[v].insert(u);
            }

            ss_line.clear();
        }
        for(int i=0; i<n_genes; i++)
            adjList[i] = vector<int>(tmpList[i].begin(), tmpList[i].end());

        f.close();
    }
    //*/
    void load_edges(string namef){
        string line, value;
        stringstream ss_line;
        unordered_set<vector<int>, VectorHash> tmpList;

        ifstream f(namef);
        while(getline(f, line)){
            ss_line << line;
            set<int> edge_set; //ordered
            vector<int> edge;
            while(ss_line >> value){
                int tmp = -1;
                if(gene2int.find(value) != gene2int.end())
                    tmp = gene2int[value];
                edge_set.insert(tmp);
            }
            for(int val : edge_set){
                edge.push_back(val);
            }
            if(edge.size() == 0 || edge[0] == -1){ // malformed edge
                // do nothing for now
            }
            else{
                tmpList.insert(edge);
            }
            ss_line.clear();
        }
        edgeList = vector<vector<int>>(tmpList.begin(), tmpList.end());
        f.close();
    }

};

#endif
