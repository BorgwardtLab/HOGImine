#ifndef hogimine
#define hogimine

#include "cmhhandler.hpp"
#include "data.hpp"
#include "utils.hpp"


typedef vector<vector<Bitset>> IntervalMap;
typedef hash_map<vector<int>, tuple<Bitset, double, double>, VectorHash> FrontierMap;

class HOGImine{
public:
    Data* data;
    double alpha;
    int maxlen;
    string outfilename;

    CMHHandler cmh;
    vector<vector<IntervalMap>> gene_intervals;

    vector<tuple<vector<int>, double>> testable_patterns;

    int cnt_process = 0;
    int cnt_signif = 0;

    HOGImine(){}
    ~HOGImine(){}
    HOGImine(Data* data_, double alpha_, int maxlen_){
        data = data_;
        //map = map_;
        //alpha = alpha_;
        maxlen = maxlen_;
        //outfilename = outfilename_;

        cmh = CMHHandler(data_, alpha_);
        process_genes();
    }

    void process_edges(){
        //ofstream out_stream(outfilename);
        for(vector<int> edge : data->edgeList){
            vector<vector<IntervalMap>> intv_maps(2, vector<IntervalMap>(edge.size()));
            vector<int> genes(edge.size());
            for(int i=0; i<edge.size(); i++){
                intv_maps[0][i] = gene_intervals[0][edge[i]];
                intv_maps[1][i] = gene_intervals[1][edge[i]];
                genes[i] = edge[i];
            }
            vector<tuple<vector<int>, double>> patterns_tmp = test_interval_combinations_bfs(intv_maps, genes); // only interval info
            for(auto& a : patterns_tmp){
                (get<0>(a)).insert((get<0>(a)).end(), genes.begin(), genes.end()); // add gene info
            }
            testable_patterns.insert(testable_patterns.end(), patterns_tmp.begin(), patterns_tmp.end());
        }
    }




    void process_genes(){
        gene_intervals = vector<vector<IntervalMap>>(2);
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
        for(int start=0; (bool)(data->is_additive) && (start<n_snps); start++){
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
                    if( starti >= gene_itvl[isHomo][i].size() || leni >= gene_itvl[isHomo][i][starti].size() || gene_itvl[isHomo][i][starti][leni].n == 0)continue; // check if [start, len+1] is in gene_itvl

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



    void write_significant(string output_filename, bool print_pv, bool verbose){
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
        ofstream sign_file(output_filename+"_sign.txt");
        vector<tuple<vector<int>, double>> sign_patterns;
        for(auto& a : testable_patterns){
            vector<int> pat = get<0>(a);
            double pv = get<1>(a);
            if(pv <= thr){
                sign_patterns.push_back(a);
            }
        }
        std::sort(sign_patterns.begin(), sign_patterns.end(), 
            [](tuple<vector<int>, double> const &t1, tuple<vector<int>, double> const &t2) {
                return get<1>(t1) < get<1>(t2);
            }
        );

        if(!verbose){
            for(auto& a : sign_patterns){
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
        }
        else{
            for(auto& a : sign_patterns){
                vector<int> pat = get<0>(a);
                double pv = get<1>(a);
                if(pv <= thr){
                    sign_file << pv << "; ";
                    vector<string> snpsn;
                    for(int i=0; i<pat.size()/2; i++){
                        int ii = pat[i];
                        string genen = data->gene_names[pat[pat.size()/2 + i]];
                        int startsnp = ((ii&((1<<30)-1))>>14);
                        int lenghtsnp = (ii&((1<<14)- 1));
                        int encgene = (ii>>30);
                        if(data->is_additive) sign_file <<  genen << " (" << (encgene==0 ? "recessive" : "dominant") << ")";
                        else sign_file << genen;
                        if(i < (pat.size()/2 -1)) sign_file << ", ";
                        else sign_file << "; ";

                        for(int j=0; j<lenghtsnp; j++){
                            snpsn.push_back(data->snp_names[data->gene2snps[pat[pat.size()/2 + i]][startsnp+j]]);
                        }
                    }
                    cnt_signif++;
                    for(int j=0; j<snpsn.size(); j++){
                        sign_file << snpsn[j];
                        if(j < (snpsn.size()-1)) sign_file << ", ";
                    }
                    sign_file << endl;
                }
            }
        }
        sign_file.close();
    }





};



#endif
