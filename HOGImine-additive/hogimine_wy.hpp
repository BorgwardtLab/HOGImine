#ifndef hogimine_wy
#define hogimine_wy

#include "cmhhandler.hpp"
#include "data.hpp"
#include "utils.hpp"


typedef vector<vector<Bitset>> IntervalMap;
typedef hash_map<vector<int>, tuple<Bitset, double, double>, VectorHash> FrontierMap;

class HOGImineWY{
public:
    Data* data;
    double alpha;
    int maxlen;
    string outfilename;

    CMHPermutationHandler cmh;
    vector<vector<IntervalMap>> gene_intervals;

    vector<tuple<vector<int>, double>> testable_patterns;

    int cnt_process = 0;
    int cnt_signif = 0;

    int onlyClosed = 1;

    HOGImineWY(){}
    ~HOGImineWY(){}
    HOGImineWY(Data* data_, double alpha_, int maxlen_, int nperm_){
        data = data_;
        //map = map_;
        //alpha = alpha_;
        maxlen = maxlen_;
        //outfilename = outfilename_;

        cmh = CMHPermutationHandler(data_, alpha_, nperm_);
        process_genes();
    }

    void process_edges(){
        //ofstream out_stream("pval.txt");
        for(vector<int> edge : data->edgeList){
            vector<vector<IntervalMap>> intv_maps(2, vector<IntervalMap>(edge.size()));
            vector<int> genes(edge.size());
            for(int i=0; i<edge.size(); i++){
                intv_maps[0][i] = gene_intervals[0][edge[i]];
                intv_maps[1][i] = gene_intervals[1][edge[i]];
                genes[i] = edge[i];
                //cout << data->gene_names[genes[i]] << " ";
            }
            //cout << endl;
            //*
            vector<tuple<vector<int>, double>> patterns_tmp = test_interval_combinations_bfs(intv_maps, genes); // only interval info
            for(auto& a : patterns_tmp){
                (get<0>(a)).insert((get<0>(a)).end(), genes.begin(), genes.end()); // add gene info
            }
            testable_patterns.insert(testable_patterns.end(), patterns_tmp.begin(), patterns_tmp.end());
            //*/
            //test_interval_combinations_dfs_helper(intv_maps, genes, out_stream);
        }
        //cout << cnt_process << endl;
    }
    /*
    void process_edges2(){
        ofstream out_stream(outfilename);
        for(int u=0; u<data->n_genes; u++){
            for(int v : data->adjList[u]){
                if(v < u) continue;

                //cout << "Edge " << data->gene_names[u] << " " << data->gene_names[v] << endl;
                vector<IntervalMap> v1(2);
                vector<int> v2(2);
                v1[0] = gene_intervals[u]; v1[1] = gene_intervals[v];
                v2[0] = u; v2[1] = v;
                //test_interval_combinations_helper(v1, v2, out_stream);
                //test_interval_pairs(v1, v2, out_stream);
                vector<tuple<vector<int>, double>> patterns_tmp = test_interval_combinations_bfs(v1, v2); // only interval info
                for(auto& a : patterns_tmp){
                    (get<0>(a)).insert((get<0>(a)).end(), v2.begin(), v2.end()); // add gene info
                }
                testable_patterns.insert(testable_patterns.end(), patterns_tmp.begin(), patterns_tmp.end());
            }
        }
        //cout << cnt_process << endl;

    }
    //*/



    void process_genes(){
        gene_intervals = vector<vector<IntervalMap>>(data->n_genes);
        gene_intervals[0] = vector<IntervalMap>(data->n_genes);
        gene_intervals[1] = vector<IntervalMap>(data->n_genes);
        for(int i=0; i<data->n_genes; i++){
            gene_intervals[0][i] = make_gene_intervals_bfs(i, 0);
            gene_intervals[1][i] = make_gene_intervals_bfs(i, 1);
        }
        //cout << cnts << endl;
    }

    IntervalMap make_gene_intervals_bfs(int gene, int useSecond){ // this can return only closed patterns, dfs cannot
        pair<int, int> snps = make_pair( data->gene2snps[gene][0], data->gene2snps[gene].size() ); // start, length
        //cout << gene << " : " << snps.second << endl;
        IntervalMap intv_map = vector<vector<Bitset>>(snps.second);
        vector<vector<int>> non_closed = vector<vector<int>>(snps.second);
        for(int start = 0; start < snps.second; start++){
            intv_map[start] = vector<Bitset>(snps.second - start);
            non_closed[start] = vector<int>(snps.second - start);
            if(!useSecond)
                intv_map[start][0] = data->matrix[snps.first+start];
            else
                intv_map[start][0] = data->matrix2[snps.first+start];
        }
        for(int len=2; len<=snps.second; len++){
            if (maxlen > 0 && len > maxlen)
                break;
            for(int start = 0; start+len-1 < snps.second; start++){
                intv_map[start][len-1] = intv_map[start][len-2].bitwise_or(intv_map[start+1][len-2]);
                int cnt = intv_map[start][len-1].count();
                if(cnt == intv_map[start][len-2].count() || cnt == intv_map[start+1][len-2].count())
                    non_closed[start][len-1] = 1;
                //else cnts++;
            }
        }
        for(int len=2; len<=snps.second; len++){
            if (maxlen > 0 && len == maxlen)
                break;
            for(int start = 0; start+len-1 < snps.second; start++){
                if(non_closed[start][len-1]){
                    intv_map[start][len-1] = Bitset();
                }
            }
        }
        return intv_map;
    }

    /*
    int cnts = 0;

    IntervalMap make_gene_intervals_dfs(int gene){
        pair<int, int> snps = make_pair( data->gene2snps[gene][0], data->gene2snps[gene].size() ); // start, length
        IntervalMap intv_map = vector<vector<Bitset>>(snps.second);

        for(int start = 0; start < snps.second; start++){
            intv_map[start] = vector<Bitset>(snps.second - start);
            Bitset support = Bitset(data->n_trans);
            int cnt = 0, tmp; // count support

            for(int len = 1; len <= snps.second - start; len++){
                if(maxlen > 0 && len > maxlen)
                    break;

                support = support.bitwise_or(data->matrix[snps.first+start+len-1]);
                if((tmp = support.count()) == cnt){ // this is not enough to get only closed intervals
                    continue; // not closed
                }
                cnt = tmp;

                vector<int> pc_supports = cmh.per_class_supports(support);

                double env = cmh.compute_envelope(pc_supports);
                cnts++;
                //cout << gene << " " << start << " " << len << " : " << env << endl;
                if(!cmh.is_prunable(env)){
                    intv_map[start][len-1] = support;
                }
            }
        }
        return intv_map;
    }


    void test_interval_combinations_dfs_helper(vector<IntervalMap>& gene_itvl, vector<int>& gene_ids, ofstream& out_file){
        vector<pair<int, int>> intv(gene_itvl.size());
        Bitset supp;
        test_interval_combinations_dfs_rec(gene_itvl, gene_ids, intv, supp, out_file, 0, gene_itvl.size());
    }


    int test_interval_combinations_dfs_rec(vector<IntervalMap>& gene_itvl, vector<int>& gene_ids, vector<pair<int, int>>& intervals, Bitset& support,
                                            ofstream& out_file, int depth, int maxdepth)
    {
        if(depth == maxdepth){
            double pvalue = 1.0;
            double minpv = 1.0;
            double env = 1.0;
            cnt_process++;
            int cas = cmh.process_pattern_env(support, &pvalue, &minpv, &env);
            if(cas == 0){ // testable
                out_file << pvalue << "," << minpv;
                for(int i=0; i<maxdepth; i++){
                    out_file << "," << data->gene_names[gene_ids[i]] << "," << intervals[i].first << "_" << intervals[i].second;
                }
                out_file << endl;
            }
            if(cas == 2) return 1; // prunable
            return 0; // not prunable
        }

        Bitset old_supp = support;
        int n_snps = data->gene2snps[gene_ids[depth]].size();

        for(int start=0; start<n_snps; start++){
            for(int len=1; len <= n_snps-start; len++){
                if(gene_itvl[depth][start][len-1].n > 0){
                    if(depth == 0){
                        support = gene_itvl[depth][start][len-1];
                    }
                    else{
                        support = support.bitwise_or(gene_itvl[depth][start][len-1]);
                    }

                    if(depth < maxdepth-1){
                        vector<int> per_class_supps = cmh.per_class_supports(support);
                        double envelope = cmh.compute_envelope(per_class_supps);
                        if(envelope > cmh.deltas[cmh.delta_id]){
                            break;
                        }
                    }


                    intervals[depth] = make_pair(start, len);
                    int prunable = test_interval_combinations_dfs_rec(gene_itvl, gene_ids, intervals, support, out_file, depth+1, maxdepth);
                    support = old_supp;
                    if(prunable) // this works only on the last level of recursion
                        break;
                }
            }
        }
        return 0;
    }

    void test_interval_pairs(vector<IntervalMap>& gene_itvl, vector<int>& gene_ids, ofstream& out_file){
        int n_snps0 = data->gene2snps[gene_ids[0]].size();
        int n_snps1 = data->gene2snps[gene_ids[1]].size();

        for(int start0=0; start0<n_snps0; start0++){
            vector<int> prunable(n_snps1, n_snps1+1);
            for(int len0=1; len0 <= n_snps0-start0; len0++){
                if(gene_itvl[0][start0][len0-1].n == 0) continue;

                for(int start1=0; start1<n_snps1; start1++){
                    for(int len1=1; len1 <= n_snps1-start1; len1++){
                        if (len1 >= prunable[start1]) break;
                        if(gene_itvl[1][start1][len1-1].n == 0)continue;

                        Bitset support = gene_itvl[0][start0][len0-1].bitwise_or(gene_itvl[1][start1][len1-1]);
                        double pvalue = 1.0;
                        double minpv = 1.0;
                        double env = 1.0;
                        int cas = cmh.process_pattern_env(support, &pvalue, &minpv, &env);
                        cnt_process++;
                        if(cas == 0){ // testable
                            out_file << pvalue << "," << minpv;
                            out_file << "," << data->gene_names[gene_ids[0]] << "," << start0 << "_" << len0;
                            out_file << "," << data->gene_names[gene_ids[1]] << "," << start1 << "_" << len1;
                            out_file << endl;
                        }
                        if(cas == 2){
                            prunable[start1] = len1;
                            break;
                        }
                    }
                }
            }
        }
    }

    //*/






    void create_len1_intervals(FrontierMap& frontier, vector<vector<IntervalMap>>& gene_itvl, vector<int>& gene_ids,
                                vector<tuple<int, int, int>>& intervals, Bitset& support,int depth, int maxdepth)
    {
        if(depth == maxdepth){
            double pvalue = 1.0;
            double minpv = 0.0;
            double env = 0.0;
            cnt_process++;
            int cas = cmh.process_pattern_env(support, &pvalue, &minpv, &env);
            if(cas < 2){ // not prunable
                vector<int> key(maxdepth);
                for(int i=0; i<maxdepth; i++) key[i] = (get<2>(intervals[i])<<30) + (get<0>(intervals[i])<<14) + (get<1>(intervals[i])); // encode as isHomo<<30+start<<14+len
                frontier[key] = make_tuple(support, pvalue, env);
            }
            return;
        }
        Bitset old_supp = support;
        int n_snps = data->gene2snps[gene_ids[depth]].size();

        // heterozygote
        for(int start=0; start<n_snps; start++){
            if(depth == 0){
                support = gene_itvl[0][depth][start][0];
            }
            else{
                support = support.bitwise_or(gene_itvl[0][depth][start][0]);
            }
            intervals[depth] = make_tuple(start, 1, 0); // isHomo = 0
            create_len1_intervals(frontier, gene_itvl, gene_ids, intervals, support, depth+1, maxdepth);
            support = old_supp;
        }
        // homozygote
        for(int start=0; start<n_snps; start++){
            if(depth == 0){
                support = gene_itvl[1][depth][start][0];
            }
            else{
                support = support.bitwise_or(gene_itvl[1][depth][start][0]);
            }
            intervals[depth] = make_tuple(start, 1, 1); //isHomo = 1
            create_len1_intervals(frontier, gene_itvl, gene_ids, intervals, support, depth+1, maxdepth);
            support = old_supp;
        }
    }


    vector<tuple<vector<int>, double>> test_interval_combinations_bfs(vector<vector<IntervalMap>>& gene_itvl, vector<int>& gene_ids){
        int n_genes = gene_itvl[0].size();
        FrontierMap frontier;
        vector<tuple<vector<int>, double>> out;

        vector<tuple<int, int, int>> intv(n_genes);
        Bitset supp;
        create_len1_intervals(frontier, gene_itvl, gene_ids, intv, supp, 0, n_genes);
        for(auto& a : frontier){
            out.push_back(make_tuple(a.first, get<1>(a.second)));
        }

        int last_front = 0;
        while(1){
            FrontierMap new_frontier;
            for(auto& pattern : frontier){
                vector<int> key = pattern.first;
                if(get<2>(pattern.second) > cmh.deltas[cmh.delta_id]) continue; // env > delta, prunable

                int len0 = key[0]&((1<<14)- 1);

                for(int i=0; i<n_genes; i++){ // gene to enlarge interval
                    int keyi = key[i];

                    int starti = (keyi&((1<<30)-1))>>14; int leni = keyi&((1<<14)- 1); int isHomo = keyi>>30;
                    if(i > 0 && len0 > 1) continue; // this is to avoid generating the same pattern multiple times, but still generate all patterns
                    if( gene_itvl[isHomo][i][starti][leni].n == 0)continue; // check if [start, len+1] is in gene_itvl

                    key[i] += (1<<14); // increase by 1 starting position and keep len
                    auto it = frontier.find(key);
                    if(it != frontier.end()){ // we try to grow pattern
                        auto& pattern2 = it->second;

                        Bitset support = (get<0>(pattern.second)).bitwise_or(get<0>(pattern2));
                        int cnt = support.count();
                        if(cnt != (get<0>(pattern.second)).count() && cnt != (get<0>(pattern2)).count() ){ // pattern is closed
                            double pvalue = 1.0;
                            double minpv = 0.0;
                            double env = 0.0;
                            cnt_process++;
                            int cas = cmh.process_pattern_env(support, &pvalue, &minpv, &env);
                            if(cas < 2){ // not prunable
                                vector<int> key2 = key;
                                key2[i] = keyi; //restore key
                                key2[i] += 1; // [st, len+1]
                                new_frontier[key2] = make_tuple(support, pvalue, env);
                            }
                        }
                        else{ // pattern not closed, can be inserted in frontier to generate new patterns but does not count towards testable
                            vector<int> per_class_supps = cmh.per_class_supports(support);
                            double envelope = cmh.compute_envelope(per_class_supps);
                            if(envelope <= cmh.deltas[cmh.delta_id]){ // not prunable, we insert in frontier but not update testable
                                vector<int> key2 = key;
                                key2[i] = keyi; //restore key
                                key2[i] += 1; // [st, len+1]
                                new_frontier[key2] = make_tuple(support, 1.0, envelope);
                            }
                        }
                    }
                    key[i] = keyi; //restore key
                }
            }
            if(new_frontier.size() == 0){
                break;
            }
            else{
                for(auto& a : new_frontier){
                    out.push_back(make_tuple(a.first, get<1>(a.second)));
                }
                frontier = new_frontier;
                last_front++;
            }
        }
        return out;
    }



    void write_significant(string output_filename, bool print_pv){
        ofstream sign_file(output_filename+"_sign.txt");
        if(print_pv) {
            ofstream pv_file(output_filename+"_pv.txt");
            for(auto& a : testable_patterns){
                vector<int> pat = get<0>(a);
                double pv = get<1>(a);
                pv_file << pv;
                for(int i=0; i<pat.size()/2; i++){ // first half of pat holds intervals, secon half gene ids
                    int ii = pat[i];
                    pv_file << "," << data->gene_names[pat[pat.size()/2 + i]] << "," << ((ii&((1<<30)-1))>>14) << "_" << (ii&((1<<14)- 1)) << "#" << (ii>>30);
                }
                pv_file << endl;
            }
            pv_file.close();
        }
        double thr = cmh.corr_threshold();
        for(auto& a : testable_patterns){
            vector<int> pat = get<0>(a);
            double pv = get<1>(a);
            if(pv <= thr){
                sign_file << pv;
                for(int i=0; i<pat.size()/2; i++){
                    int ii = pat[i];
                    sign_file << "," << data->gene_names[pat[pat.size()/2 + i]] << "," << ((ii&((1<<30)-1))>>14) << "_" << (ii&((1<<14)- 1)) << "#" << (ii>>30);
                }
                cnt_signif++;
                sign_file << endl;
            }
        }
        sign_file.close();
    }



    /*
      Reads the significant itemsets from the pvalue.csv file of a Tarone run.
    int write_significant_edge_epi(double threshold,
                                   string input_filename,
                                   string output_filename)
    {

      ifstream in_file(input_filename);
      ofstream out_file(output_filename);

      out_file << "p-value,min_pval,gene0,start0_len0,gene1,start1_len1" ;
      out_file << endl;

      string line;
      string cell;
      getline(in_file,line);
      long long sig_count = 0;

      while(getline(in_file, line))
      {
        istringstream ss_line(line);
        while(getline(ss_line, cell, ','))
        {
          if (stod(cell) < threshold)
          {
            sig_count ++;
            out_file << line << endl;
          }
          break;
        }
        ss_line.clear();
      }
      return sig_count;
    }
    //*/

};



#endif
