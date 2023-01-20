#ifndef cmhhandler_h
#define cmhhandler_h

#include "utils.hpp"
#include "data.hpp"


#define N_DELTA 1000

class CMHHandler{
public:
    double alpha;
    Data* data;

    vector<long long> mask_start, mask_end;
    vector<double> gammat, gammabint;
    vector<int> hypercorner_bnd;

    vector<double> deltas;
    vector<int> testable_at_delta;
    int testable;
    int delta_id;



    CMHHandler(){}
    ~CMHHandler(){}
    CMHHandler(Data* data_, double alpha_){
        data = data_;
        alpha = alpha_;
        testable_at_delta = vector<int>(N_DELTA);
        deltas = vector<double>(N_DELTA);
        deltas[0] = 1.0;
        for(int i=1; i<N_DELTA; i++){
            deltas[i] = deltas[i-1]*0.95;
        }
        delta_id = 0; testable = 0;

        // precompute masks for per_class_supports
        mask_start = vector<long long>(data->num_covariates);
        mask_end = vector<long long>(data->num_covariates);
        for(int c=0; c<data->num_covariates; c++){
            int start = data->covariate_start[c];
            int end = data->covariate_start[c+1]-1; //inclusive
            if((start>>6) == (end>>6)){
                if((start&63) > 0)
                    mask_start[c] = ((1ull<<((end&63)+1))-1)^( (1ull<<(start&63))-1 );
                else
                    mask_start[c] = ((1ull<<((end&63)+1))-1);
                mask_end[c] = 0ll;
            }
            else{
                if((start&63) > 0)
                    mask_start[c] = (-1)^( (1ll<<(start&63))-1 );
                else
                    mask_start[c] = (-1);
                mask_end[c] = ((1ull<<((end&63)+1))-1);
            }
        }

        // precompute gamma
        for (int c=0; c<data->num_covariates; c++){
            gammat.push_back((double)data->positive_per_class[c] / data->trans_per_class[c]);
            gammabint.push_back(gammat[c]*(1-gammat[c]));

            int tmp_c = data->trans_per_class[c] - data->positive_per_class[c];  // number of negatives
            int maj = data->positive_per_class[c] > tmp_c ? data->positive_per_class[c] : tmp_c;
            hypercorner_bnd.push_back(maj);
          }

    }

    int process_pattern(Bitset& support, double* pval, double* minpv){
        vector<int> per_class_supps = per_class_supports(support);
        double min_pv = compute_minpval(per_class_supps);
        *minpv = min_pv;
        if( min_pv <= deltas[delta_id] ){ // is currently testable
            process_testable(min_pv);
        }
        if( min_pv <= deltas[delta_id] ){ // is still testable
            int a = compute_ap(support);
            double pvalue = compute_pval(a, per_class_supps);
            *pval = pvalue;
            return 0; // testable (and not prunable)
        }
        else{
            double envelope = compute_envelope(per_class_supps);
            if( envelope <= deltas[delta_id] )
                return 1; // not testable but not prunable
            else
                return 2; // not testable and prunable
        }
    }

    int process_pattern_env(Bitset& support, double* pval, double* minpv, double* env){
        vector<int> per_class_supps = per_class_supports(support);
        double min_pv = compute_minpval(per_class_supps);
        double envelope = compute_envelope(per_class_supps);
        *minpv = min_pv;
        *env = envelope;
        if( min_pv <= deltas[delta_id] ){ // is currently testable
            process_testable(min_pv);
        }
        if( min_pv <= deltas[delta_id] ){ // is still testable
            int a = compute_ap(support);
            double pvalue = compute_pval(a, per_class_supps);
            *pval = pvalue;
            return 0; // testable (and not prunable)
        }
        else{
            if( envelope <= deltas[delta_id] )
                return 1; // not testable but not prunable
            else
                return 2; // not testable and prunable
        }
    }

    void process_testable(double min_pv){
        int id;
        if(min_pv < 1e-40){
            id = N_DELTA-1;
        }
        else{
            id = min(N_DELTA-1, (int)floor(log(min_pv)/log(0.8)));
        }
        testable_at_delta[id]++;
        testable++;
        while( testable * deltas[delta_id] > alpha ){
            testable -= testable_at_delta[delta_id];
            delta_id++;
        }
    }

    bool is_testable(double min_pv){
        //return true;
        return ( min_pv <= deltas[delta_id] );
    }
    bool is_prunable(double env){
        //return false;
        return ( env > deltas[delta_id] );
    }

    double corr_threshold(){
        return alpha / (testable+0.0);
    }




    vector<int> per_class_supports_slow(Bitset& support){
        vector<int> out = vector<int>(data->num_covariates);
        for(int i=0; i<data->n_trans; i++){
            if(support.get(i)){
                out[data->covars[i]]++;
            }
        }
        return out;
    }

    vector<int> per_class_supports(Bitset& support){
        vector<int> out = vector<int>(data->num_covariates);
        for(int c=0; c<data->num_covariates; c++){
            int supp = 0;
            int i = (data->covariate_start[c])>>6, end = (data->covariate_start[c+1]-1)>>6;
            supp += __builtin_popcountll( support.bits[i] & mask_start[c]);
            i++;
            while( i < end ){
                supp += __builtin_popcountll(support.bits[i]);
                i++;
            }
            supp += __builtin_popcountll( support.bits[i] & mask_end[c]);
            out[c] = supp;
        }
        return out;
    }


    int compute_ap(Bitset& support){
        return (support.bitwise_and(data->labels)).count();
    }

    double compute_pval(int a, vector<int>& per_class_supps){
        double num = a;
        double den = 0.0;
        for(int c=0; c<per_class_supps.size(); c++){
            num -= per_class_supps[c] * gammat[c];
            den += per_class_supps[c] * (1.0 - ((double)per_class_supps[c])/data->trans_per_class[c]) * gammabint[c];
        }
        if(fabs(den) < 1e-10){
            return 1.0;
        }
        else{
            return Chi2_sf((num*num)/den, 1);
        }
    }

    double compute_minpval(vector<int>& per_class_supps){
        // Holds the numerator of T_amin
        double left_tail_num = 0;
        // Holds the numberator of T_amax
        double right_tail_num = 0;
        // Holds the denominator.
        double den = 0.0;
        // Auxiliaries for computation of T-scores.
        double aux1, aux2;

        for(int c=0; c<per_class_supps.size(); c++){
            aux1 = per_class_supps[c] - (data->trans_per_class[c] - data->positive_per_class[c]);
            aux2 = per_class_supps[c] * gammat[c];
            left_tail_num += ((aux1 > 0) ? aux1 : 0) - aux2;
            right_tail_num += ((per_class_supps[c] > data->positive_per_class[c]) ? data->positive_per_class[c] : per_class_supps[c]) - aux2;
            den += per_class_supps[c] * (1-((double)per_class_supps[c])/data->trans_per_class[c]) * gammabint[c];
        }

        if(fabs(den) < 1e-10){
            return 1.0;
        }
        else{
            left_tail_num *= left_tail_num;
            right_tail_num *= right_tail_num;
            double minpval = Chi2_sf(((left_tail_num > right_tail_num) ? left_tail_num : right_tail_num)/den, 1);
            return minpval;
        }
    }

    double compute_envelope(vector<int>& per_class_supps){
        std::vector<double> f_vals;
        std::vector<double> g_vals;
        std::vector<pair<double, int>> beta;
        std::vector<int> idx_beta_sorted;

        double f_sum = 0;
        double g_sum = 0;
        double Tcmh_aux_corner;
        double Tcmh_max_corner_r;
        double Tcmh_max_corner_l;


        // If for any of the tables, its margin is smaller than the maximum of n
        // and N-n, then we cannot prune the set (we are not in the "top-right"
        // hypercorner)
        for(int c=0; c<per_class_supps.size(); c++) {
            //cout << c << " :: " << per_class_supps[c] << " " << hypercorner_bnd[c] << endl;
            if(per_class_supps[c] < hypercorner_bnd[c]) return 0.0;
        }

        // compute the righthandside values.
        for(int c=0; c<per_class_supps.size(); c++) {
            if(per_class_supps[c] < data->trans_per_class[c]){
                double f = gammat[c] * (data->trans_per_class[c]-per_class_supps[c]);
                double g = gammabint[c] * per_class_supps[c] * (1-((double)per_class_supps[c])/data->trans_per_class[c]);
                f_vals.push_back(f);
                g_vals.push_back(g);
                beta.push_back({g/f, beta.size()});
            }
        }
        sort(beta.begin(), beta.end());
        // compute CMH test statistic.
        // Skip tables that only have one class, i.e. g_sum == 0.
        f_sum = 0;
        g_sum = 0;
        Tcmh_max_corner_r = 0;
        for(int k=0; k<beta.size(); k++){
            f_sum += f_vals[beta[k].second];
            g_sum += g_vals[beta[k].second];
            if (g_sum == 0){
                continue;
            }
            Tcmh_aux_corner = (f_sum*f_sum)/g_sum;
            Tcmh_max_corner_r = (Tcmh_max_corner_r >= Tcmh_aux_corner) ? Tcmh_max_corner_r : Tcmh_aux_corner;
        }

        g_vals.clear();
        f_vals.clear();
        beta.clear();


        // compute the lefthandside values.
        for(int c=0; c<per_class_supps.size(); c++) {
            if(per_class_supps[c] < data->trans_per_class[c]){
                double f = (1-gammat[c]) * (data->trans_per_class[c]-per_class_supps[c]);
                double g = gammabint[c] * per_class_supps[c] * (1-((double)per_class_supps[c])/data->trans_per_class[c]);
                f_vals.push_back(f);
                g_vals.push_back(g);
                beta.push_back({g/f, beta.size()});
            }
        }
        sort(beta.begin(), beta.end());

        f_sum = 0;
        g_sum = 0;
        Tcmh_max_corner_l = 0;
        for(int k=0; k<beta.size(); k++){
            f_sum += f_vals[beta[k].second];
            g_sum += g_vals[beta[k].second];
            if (g_sum == 0){
                continue;
            }
            Tcmh_aux_corner = (f_sum*f_sum)/g_sum;
            Tcmh_max_corner_l = (Tcmh_max_corner_l >= Tcmh_aux_corner) ? Tcmh_max_corner_l : Tcmh_aux_corner;
        }

        //cout << "Tcmh: " << ((Tcmh_max_corner_r >= Tcmh_max_corner_l) ? Tcmh_max_corner_r : Tcmh_max_corner_l) << endl;
        double min_pval = Chi2_sf(( (Tcmh_max_corner_r >= Tcmh_max_corner_l) ? Tcmh_max_corner_r : Tcmh_max_corner_l), 1);
        return min_pval;
    }

};









class CMHPermutationHandler{
public:
    double alpha;
    int nperm;
    Data* data;

    vector<long long> mask_start, mask_end;
    vector<double> gammat, gammabint;
    vector<int> hypercorner_bnd;

    vector<double> deltas;
    vector<int> testable_at_delta;
    int testable;
    int delta_id;

    vector<double> perm_minpv;
    vector<Bitset> perm_labels;

    CMHPermutationHandler(){}
    ~CMHPermutationHandler(){}
    CMHPermutationHandler(Data* data_, double alpha_, double _nperm){
        srand(0); // change this
        data = data_;
        alpha = alpha_;
        nperm = _nperm;
        testable_at_delta = vector<int>(N_DELTA);
        deltas = vector<double>(N_DELTA);
        deltas[0] = 1.0;
        for(int i=1; i<N_DELTA; i++){
            deltas[i] = deltas[i-1]*0.95;
        }
        delta_id = 0; testable = 0;
        // init minpv for each permutation
        perm_minpv = vector<double>(nperm, 1.0);
        // generate permuted labels
        perm_labels = vector<Bitset>(nperm, Bitset(data->n_trans));
        for(int i=0; i<nperm; i++){
            int id = 0;
            for(int c=0; c<data->num_covariates; c++){
                vector<int> indexes(data->trans_per_class[c]);
                for(int j=0; j<indexes.size(); j++) indexes[j] = j;
                random_shuffle(indexes.begin(), indexes.end());
                for(int j=0; j<indexes.size(); j++){
                    if(data->labels.get(id+j))
                        perm_labels[i].set(id+indexes[j]);
                }
                id += indexes.size();
            }

        }

        // precompute masks for per_class_supports
        mask_start = vector<long long>(data->num_covariates);
        mask_end = vector<long long>(data->num_covariates);
        for(int c=0; c<data->num_covariates; c++){
            int start = data->covariate_start[c];
            int end = data->covariate_start[c+1]-1; //inclusive
            if((start>>6) == (end>>6)){
                if((start&63) > 0)
                    mask_start[c] = ((1ull<<((end&63)+1))-1)^( (1ull<<(start&63))-1 );
                else
                    mask_start[c] = ((1ull<<((end&63)+1))-1);
                mask_end[c] = 0ll;
            }
            else{
                if((start&63) > 0)
                    mask_start[c] = (-1)^( (1ll<<(start&63))-1 );
                else
                    mask_start[c] = (-1);
                mask_end[c] = ((1ull<<((end&63)+1))-1);
            }
        }

        // precompute gamma
        for (int c=0; c<data->num_covariates; c++){
            gammat.push_back((double)data->positive_per_class[c] / data->trans_per_class[c]);
            gammabint.push_back(gammat[c]*(1-gammat[c]));

            int tmp_c = data->trans_per_class[c] - data->positive_per_class[c];  // number of negatives
            int maj = data->positive_per_class[c] > tmp_c ? data->positive_per_class[c] : tmp_c;
            hypercorner_bnd.push_back(maj);
          }

    }




    int process_pattern_env(Bitset& support, double* pval, double* minpv, double* env){
        vector<int> per_class_supps = per_class_supports(support);
        double min_pv = compute_minpval(per_class_supps);
        double envelope = compute_envelope(per_class_supps);
        *minpv = min_pv;
        *env = envelope;
        if( min_pv <= deltas[delta_id] ){ // is currently testable
            update_minpvals(support, per_class_supps);
            while(empirical_fwer() > alpha && delta_id < N_DELTA-1){
                delta_id++; // decrease significance threshold
            }
        }
        if( min_pv <= deltas[delta_id] ){ // is still testable
            int a = compute_ap(support);
            double pvalue = compute_pval(a, per_class_supps);
            *pval = pvalue;
            return 0; // testable (and not prunable)
        }
        else{
            if( envelope <= deltas[delta_id] )
                return 1; // not testable but not prunable
            else
                return 2; // not testable and prunable
        }
    }

    void update_minpvals(Bitset& support, vector<int>& per_class_supps){
        // compute possible values for a
        int a_min = 0;
        int a_max = 0;
        for (int i=0; i<data->num_covariates; i++){
            a_min += max(0, per_class_supps[i] - data->trans_per_class[i] + data->positive_per_class[i]);
            a_max += min(per_class_supps[i], data->positive_per_class[i]);
        }

        // precompute pvals
        vector<double> pre_pvals(a_max-a_min+1);
        for(int a=a_min; a<=a_max; a++)
            pre_pvals[a-a_min] = compute_pval(a, per_class_supps);
        // update minpvals
        for (int i=0; i<nperm; i++){
            int a = (support.bitwise_and(perm_labels[i])).count();
            perm_minpv[i] = min(perm_minpv[i], pre_pvals[a-a_min]);
        }

    }

    double empirical_fwer(){
        int count = 0;
        for (int i=0; i<nperm; i++){
            if (perm_minpv[i] <= deltas[delta_id]){
                count++;
            }
        }
        return ((double)count)/nperm;
    }




    double corr_threshold(){
        //return deltas[delta_id];
        vector<double> tmp(perm_minpv.begin(), perm_minpv.end());
        sort(tmp.begin(), tmp.end());
        return tmp[perm_minpv.size()*alpha];
    }

    bool is_prunable(double env){
        return ( env > deltas[delta_id] );
    }

    vector<int> per_class_supports(Bitset& support){
        vector<int> out = vector<int>(data->num_covariates);
        for(int c=0; c<data->num_covariates; c++){
            int supp = 0;
            int i = (data->covariate_start[c])>>6, end = (data->covariate_start[c+1]-1)>>6;
            supp += __builtin_popcountll( support.bits[i] & mask_start[c]);
            i++;
            while( i < end ){
                supp += __builtin_popcountll(support.bits[i]);
                i++;
            }
            supp += __builtin_popcountll( support.bits[i] & mask_end[c]);
            out[c] = supp;
        }
        return out;
    }

    int compute_ap(Bitset& support){
        return (support.bitwise_and(data->labels)).count();
    }

    double compute_pval(int a, vector<int>& per_class_supps){
        double num = a;
        double den = 0.0;
        for(int c=0; c<per_class_supps.size(); c++){
            num -= per_class_supps[c] * gammat[c];
            den += per_class_supps[c] * (1.0 - ((double)per_class_supps[c])/data->trans_per_class[c]) * gammabint[c];
        }
        if(fabs(den) < 1e-10){
            return 1.0;
        }
        else{
            return Chi2_sf((num*num)/den, 1);
        }
    }

    double compute_minpval(vector<int>& per_class_supps){
        // Holds the numerator of T_amin
        double left_tail_num = 0;
        // Holds the numberator of T_amax
        double right_tail_num = 0;
        // Holds the denominator.
        double den = 0.0;
        // Auxiliaries for computation of T-scores.
        double aux1, aux2;

        for(int c=0; c<per_class_supps.size(); c++){
            aux1 = per_class_supps[c] - (data->trans_per_class[c] - data->positive_per_class[c]);
            aux2 = per_class_supps[c] * gammat[c];
            left_tail_num += ((aux1 > 0) ? aux1 : 0) - aux2;
            right_tail_num += ((per_class_supps[c] > data->positive_per_class[c]) ? data->positive_per_class[c] : per_class_supps[c]) - aux2;
            den += per_class_supps[c] * (1-((double)per_class_supps[c])/data->trans_per_class[c]) * gammabint[c];
        }

        if(fabs(den) < 1e-10){
            return 1.0;
        }
        else{
            left_tail_num *= left_tail_num;
            right_tail_num *= right_tail_num;
            double minpval = Chi2_sf(((left_tail_num > right_tail_num) ? left_tail_num : right_tail_num)/den, 1);
            return minpval;
        }
    }

    double compute_envelope(vector<int>& per_class_supps){
        std::vector<double> f_vals;
        std::vector<double> g_vals;
        std::vector<pair<double, int>> beta;
        std::vector<int> idx_beta_sorted;

        double f_sum = 0;
        double g_sum = 0;
        double Tcmh_aux_corner;
        double Tcmh_max_corner_r;
        double Tcmh_max_corner_l;


        // If for any of the tables, its margin is smaller than the maximum of n
        // and N-n, then we cannot prune the set (we are not in the "top-right"
        // hypercorner)
        for(int c=0; c<per_class_supps.size(); c++) {
            //cout << c << " :: " << per_class_supps[c] << " " << hypercorner_bnd[c] << endl;
            if(per_class_supps[c] < hypercorner_bnd[c]) return 0.0;
        }

        // compute the righthandside values.
        for(int c=0; c<per_class_supps.size(); c++) {
            if(per_class_supps[c] < data->trans_per_class[c]){
                double f = gammat[c] * (data->trans_per_class[c]-per_class_supps[c]);
                double g = gammabint[c] * per_class_supps[c] * (1-((double)per_class_supps[c])/data->trans_per_class[c]);
                f_vals.push_back(f);
                g_vals.push_back(g);
                beta.push_back({g/f, beta.size()});
            }
        }
        sort(beta.begin(), beta.end());
        // compute CMH test statistic.
        // Skip tables that only have one class, i.e. g_sum == 0.
        f_sum = 0;
        g_sum = 0;
        Tcmh_max_corner_r = 0;
        for(int k=0; k<beta.size(); k++){
            f_sum += f_vals[beta[k].second];
            g_sum += g_vals[beta[k].second];
            if (g_sum == 0){
                continue;
            }
            Tcmh_aux_corner = (f_sum*f_sum)/g_sum;
            Tcmh_max_corner_r = (Tcmh_max_corner_r >= Tcmh_aux_corner) ? Tcmh_max_corner_r : Tcmh_aux_corner;
        }

        g_vals.clear();
        f_vals.clear();
        beta.clear();


        // compute the lefthandside values.
        for(int c=0; c<per_class_supps.size(); c++) {
            if(per_class_supps[c] < data->trans_per_class[c]){
                double f = (1-gammat[c]) * (data->trans_per_class[c]-per_class_supps[c]);
                double g = gammabint[c] * per_class_supps[c] * (1-((double)per_class_supps[c])/data->trans_per_class[c]);
                f_vals.push_back(f);
                g_vals.push_back(g);
                beta.push_back({g/f, beta.size()});
            }
        }
        sort(beta.begin(), beta.end());

        f_sum = 0;
        g_sum = 0;
        Tcmh_max_corner_l = 0;
        for(int k=0; k<beta.size(); k++){
            f_sum += f_vals[beta[k].second];
            g_sum += g_vals[beta[k].second];
            if (g_sum == 0){
                continue;
            }
            Tcmh_aux_corner = (f_sum*f_sum)/g_sum;
            Tcmh_max_corner_l = (Tcmh_max_corner_l >= Tcmh_aux_corner) ? Tcmh_max_corner_l : Tcmh_aux_corner;
        }

        //cout << "Tcmh: " << ((Tcmh_max_corner_r >= Tcmh_max_corner_l) ? Tcmh_max_corner_r : Tcmh_max_corner_l) << endl;
        double min_pval = Chi2_sf(( (Tcmh_max_corner_r >= Tcmh_max_corner_l) ? Tcmh_max_corner_r : Tcmh_max_corner_l), 1);
        return min_pval;
    }


};


#endif
