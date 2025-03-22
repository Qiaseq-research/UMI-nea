#include <condition_variable>
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
#include <iomanip>
#include <map>
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
	int max_umi_len;
	int min_read_founder;
	bool first_founder_mode;
	int verbose=1;
	string primer_id;
}UMI_clustering_parameters;

struct guardedvector {
   mutex guard;
   vector<string> myvector;
   int size_last_cycle=0;
};

typedef struct {
	string UMI_seq;
	int founder_offset=0;
	string founder_temp;
	int founder_temp_dist;
	bool founder_temp_found=false;
}UMI_item;

inline vector<int> repeat_n (int n){
        vector <int> repeats;
        for (int i=0; i<n; i++)
                repeats.push_back(n);
        return repeats;
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

// Function to convert a std::map<K,V> to std::multimap<V,K>
template<typename K, typename V>
        std::multimap<V,K,std::greater<int>> invertMap(std::map<K,V> const &map)
        {
        std::multimap<V,K,std::greater<int>> multimap;

        for (auto const &pair: map) {
                 multimap.insert(std::make_pair(pair.second, pair.first));
        }
        return multimap;
}

double mad  (const vector<int> v);

double mean(const vector<int> v);

double var(const vector<int> v);

void shared_writer(ofstream& out, const vector<string> lines);

int count_umi(const string filename);

int calculate_dist_upper_bound(float error_rate, int max_umi_len);

void do_nb(const string  updated_count_file, float nb_lowertail_p, int madfolds, int min_read_founder,  int nb_estimated_molecule, int median_rpu, ofstream & e_out_file);

void do_kp( const string updated_count_file, int  min_read_founder, int kp_estimated_molecule, int kp_angle, int median_rpu, ofstream & e_out_file  );

void fit_nb_model( const string  filename, float  p, int madfolds, int & min_read_founder, int &  nb_estimated_molecule, int & median_rpu);

void fit_knee_plot ( const string  filename, int & min_read_founder, int & kp_estimated_molecule, int & kp_angle, int & median_rpu, int & after_rpucut_molecule );

bool align_umi( string umi, string f,  int max_dist, int & endpos, int & dist);

bool producer(const vector<UMI_item> &umi_pool_subset, ofstream& out, const string primer_id, const unsigned max_dist, const unsigned max_umi_len);

vector<UMI_item> consumer(const int worker_index, ofstream& out, bool first_founder_mode, const string  primer_id, const unsigned  max_dist, const unsigned max_umi_len, map<int, vector<UMI_item>>  t_umi);

vector <int> split_umi_to_threads_on_founder(vector<UMI_item> umi_pool, int threads, int pool_size);

void parallel_processing( vector<UMI_item>& umi_pool,  ofstream& out, UMI_clustering_parameters parameters);

void founder_find ( vector<UMI_item>& umi_pool, ofstream& out,  string primer_id, const unsigned  max_dist, const unsigned max_umi_len, const int pool );

void parallel_founder_find ( vector<UMI_item> low_reads_umi_pool,  ofstream& out, const unsigned num_worker_threads, string curr_primer_id, const unsigned  max_dist, const unsigned max_umi_len  );

void clustering_umis(const string in_filename, const string out_filename,  UMI_clustering_parameters parameters );

void update_umi_reads_count(const string updated_count_filename, const string in_filename, const string out_filename);
