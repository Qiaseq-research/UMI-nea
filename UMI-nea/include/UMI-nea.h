#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <future>
#include <thread>
#include <mutex>
#include <algorithm>
#include <iostream>
#include <map>
#include <unordered_map>
#include <set>
#include <numeric>
#include <utility>
#include "edlib.h"
#include <boost/math/distributions/negative_binomial.hpp>

using namespace std;

typedef chrono::high_resolution_clock Clock;

#define PI 3.1415926535

typedef struct {
	int max_dist;
	int thread;
	int pool_size;
	int prod_size = 1000;  // UMIs handled by serial producer per batch
	int max_umi_len;
	int min_read_founder;
	bool greedy_mode;
	int verbose=1;
	string primer_id;
}UMI_clustering_parameters;

struct guardedvector {
   vector<string> myvector;
};

typedef struct {
	string UMI_seq;
	string padded_umi;  // pre-computed "NNN"+UMI_seq+"NNN", cached to avoid realloc in hot loops
	int founder_offset=0;
	string founder_temp;
	int founder_temp_dist;
	bool founder_temp_found=false;
	int input_pos=0;    // 1-based line number in input file; used to restore input order after drain merges
}UMI_item;

inline vector<int> repeat_n(int n) {
	vector<int> v;
	for (int i = 0; i < n; i++) v.push_back(n);
	return v;
}

template <typename T>
T median(const vector<T> v){
	if (v.size()%2==1)
                return v.at(v.size()/2);
        else
                return (v.at(v.size()/2)+v.at(v.size()/2-1))/2;

}

template <typename T>
vector<int> findItems(T const &v, int target) {
        vector<int> indices;
        //for (int i =0; i<v.size(); i++){
	int i = 0;
	for (auto u : v) {
                //if (v[i]==target)
		if (u == target)
                        indices.push_back(i);
		i++;
        }
        return indices;
}

// Convert any map-like container to std::multimap<V,K> sorted descending by value
template<typename Map>
std::multimap<typename Map::mapped_type, typename Map::key_type, std::greater<int>>
invertMap(Map const &map)
{
        std::multimap<typename Map::mapped_type, typename Map::key_type, std::greater<int>> multimap;
        for (auto const &pair: map)
                multimap.insert(std::make_pair(pair.second, pair.first));
        return multimap;
}

double mean(const vector<int> v);

double var(const vector<int> v);

double mad(const vector<int> v);

int count_umi(const string filename);

int calculate_dist_upper_bound(float error_rate, int max_umi_len);

int calculate_dist_upper_bound_old(float error_rate, int max_umi_len);

void do_nb(const string  updated_count_file, float nb_lowertail_p, int madfolds, int min_read_founder,  int nb_estimated_molecule, int median_rpu, ofstream & e_out_file);

void do_kp( const string updated_count_file, int  min_read_founder, int kp_estimated_molecule, int kp_angle, int median_rpu, ofstream & e_out_file  );

void fit_nb_model( const string  filename, float  p, int madfolds, int & min_read_founder, int &  nb_estimated_molecule, int & median_rpu);

void fit_nb_model_old( const string  filename, float  p, int madfolds, int & min_read_founder, int &  nb_estimated_molecule, int & median_rpu);

void fit_knee_plot ( const string  filename, int & min_read_founder, int & kp_estimated_molecule, int & kp_angle, int & median_rpu, int & after_rpucut_molecule );

bool align_umi(const string& umi, const string& f, int max_dist, int& endpos, int& dist);

string producer(const vector<UMI_item> &umi_pool_subset, const string& primer_id, const unsigned max_dist, const unsigned max_umi_len);

pair<string, vector<UMI_item>> consumer_pass(const vector<UMI_item>& umi_pool, int begin, int end, bool first_founder_mode, const string& primer_id, unsigned max_dist, unsigned max_umi_len, const array<vector<int>,4096>& kmer_idx, const vector<string>& founder_view, const string& padded_founders);

void parallel_processing( vector<UMI_item>& umi_pool,  ofstream& out, UMI_clustering_parameters parameters);

string founder_find ( const vector<UMI_item>& umi_pool, int begin, int end, const string primer_id, const unsigned  max_dist, const unsigned max_umi_len, const int pool, const array<vector<int>, 4096>& shared_kmer_idx, const vector<string>& shared_founder_view, const string& shared_padded_founders );

void parallel_founder_find ( vector<UMI_item> low_reads_umi_pool,  ofstream& out, const unsigned num_worker_threads, const string& curr_primer_id, const unsigned  max_dist, const unsigned max_umi_len  );

void clustering_umis(const string in_filename, const string out_filename,  UMI_clustering_parameters parameters );

void update_umi_reads_count(const string updated_count_filename, const string in_filename, const string out_filename);
