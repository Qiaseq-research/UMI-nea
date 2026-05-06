#pragma once

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <future>
#include <iostream>
#include <map>
#include <mutex>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

#include "edlib.h"
#include <boost/math/distributions/negative_binomial.hpp>

using Clock = std::chrono::high_resolution_clock;

static constexpr double PI = 3.1415926535;

struct UMI_clustering_parameters {
    int         max_dist;
    int         thread;
    int         pool_size;
    int         prod_size        = 1000;  // UMIs handled by serial producer per batch
    int         max_umi_len;
    int         min_read_founder;
    bool        greedy_mode;
    int         verbose          = 1;
    std::string primer_id;
};

struct guardedvector {
    std::vector<std::string> myvector;
};

struct UMI_item {
    std::string UMI_seq;
    std::string padded_umi;          // pre-computed "NNN"+UMI_seq+"NNN", cached to avoid realloc in hot loops
    int         founder_offset       = 0;
    std::string founder_temp;
    int         founder_temp_dist    = 0;
    bool        founder_temp_found   = false;
    int         input_pos            = 0;  // 1-based line number; used to restore input order after drain merges
};

template <typename T>
T median(const std::vector<T>& v) {
    if (v.size() % 2 == 1)
        return v.at(v.size() / 2);
    return (v.at(v.size() / 2) + v.at(v.size() / 2 - 1)) / 2;
}

template <typename T>
std::vector<int> findItems(const T& v, int target) {
    std::vector<int> indices;
    int i = 0;
    for (auto u : v) {
        if (u == target)
            indices.push_back(i);
        ++i;
    }
    return indices;
}

// Invert a map to multimap<V,K> sorted descending by value.
template <typename Map>
std::multimap<typename Map::mapped_type, typename Map::key_type,
              std::greater<typename Map::mapped_type>>
invertMap(const Map& m) {
    std::multimap<typename Map::mapped_type, typename Map::key_type,
                  std::greater<typename Map::mapped_type>> result;
    for (const auto& kv : m)
        result.insert({kv.second, kv.first});
    return result;
}

// --- Statistics ---
double mean(const std::vector<int>& v);
double var(const std::vector<int>& v);
double mad(const std::vector<int>& v);

// --- UMI utilities ---
int  count_umi(const std::string& filename);
int  calculate_dist_upper_bound(float error_rate, int max_umi_len);

// --- Molecule estimation ---
void fit_nb_model(const std::string& filename, float p, int madfolds,
                  int& min_read_founder, int& nb_estimated_molecule, int& median_rpu);

void fit_knee_plot(const std::string& filename, int& min_read_founder,
                   int& kp_estimated_molecule, int& kp_angle, int& median_rpu,
                   int& after_rpucut_molecule);

void do_nb(const std::string& updated_count_file, float nb_lowertail_p, int madfolds,
           int min_read_founder, int nb_estimated_molecule, int median_rpu,
           std::ofstream& e_out_file);

void do_kp(const std::string& updated_count_file, int min_read_founder,
           int kp_estimated_molecule, int kp_angle, int median_rpu,
           std::ofstream& e_out_file);

// --- Alignment ---
bool align_umi(const std::string& umi, const std::string& f,
               int max_dist, int& endpos, int& dist);

// --- Clustering pipeline ---
std::string producer(const std::vector<UMI_item>& umi_pool_subset,
                     const std::string& primer_id,
                     unsigned max_dist, unsigned max_umi_len);

std::pair<std::string, std::vector<UMI_item>> consumer_pass(
    const std::vector<UMI_item>& umi_pool, int begin, int end,
    bool first_founder_mode, const std::string& primer_id,
    unsigned max_dist, unsigned max_umi_len,
    const std::array<std::vector<int>, 4096>& kmer_idx,
    const std::vector<std::string>& founder_view,
    const std::string& padded_founders);

void parallel_processing(std::vector<UMI_item>& umi_pool, std::ofstream& out,
                         UMI_clustering_parameters parameters);

std::string founder_find(const std::vector<UMI_item>& umi_pool, int begin, int end,
                         const std::string& primer_id, unsigned max_dist,
                         unsigned max_umi_len, int pool,
                         const std::array<std::vector<int>, 4096>& shared_kmer_idx,
                         const std::vector<std::string>& shared_founder_view,
                         const std::string& shared_padded_founders);

void parallel_founder_find(std::vector<UMI_item> low_reads_umi_pool,
                           std::ofstream& out, unsigned num_worker_threads,
                           const std::string& curr_primer_id,
                           unsigned max_dist, unsigned max_umi_len);

void clustering_umis(const std::string& in_filename, const std::string& out_filename,
                     UMI_clustering_parameters parameters);

void update_umi_reads_count(const std::string& updated_count_filename,
                            const std::string& in_filename,
                            const std::string& out_filename);
