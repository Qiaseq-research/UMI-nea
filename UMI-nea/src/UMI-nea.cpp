#include "UMI-nea.h"

// Encode a 6-mer as a 12-bit int (2 bits/base). Returns -1 on N or unknown base.
static inline int encode_kmer(const char* s) {
    int code = 0;
    for (int i = 0; i < 6; ++i) {
        int b;
        switch (s[i]) {
            case 'A': case 'a': b = 0; break;
            case 'C': case 'c': b = 1; break;
            case 'G': case 'g': b = 2; break;
            case 'T': case 't': b = 3; break;
            default: return -1;
        }
        code = (code << 2) | b;
    }
    return code;
}

//global variable
atomic_int founder_added{0};
// Consumer poll loop uses founder_added atomic + 10µs sleep instead of condition_variable.
guardedvector founders;
atomic<bool> producer_done{false};
const int N_num=3;
const string padding(N_num, 'N');
bool verbose=false;

double mean(const vector<int> v)
{
      long double sum = 0;
      for (auto &each: v)
	    sum += each;
      return sum / v.size();
}

double var(const vector<int> v)
{
      long double square_sum_of_difference = 0;
      long double mean_var = mean(v);
      auto len = v.size();
      double tmp;
      for (auto &each: v) {
	    tmp = each - mean_var;
	    square_sum_of_difference += tmp * tmp;
      }
      return (square_sum_of_difference / (len -1 ) );
}

int calculate_dist_upper_bound(float error_rate, int max_umi_len){
      float z;
      int alpha_i;
      if (max_umi_len<=15)
            alpha_i=1; //alpha=0.95;
      else if (max_umi_len > 15 && max_umi_len <= 25)
            alpha_i=2; //alpha=0.99;
      else if (max_umi_len > 25 && max_umi_len <= 35)
            alpha_i=3; //alpha=0.999;
      else
            alpha_i=4; //alpha=0.9999;
      switch(alpha_i){
            case 1: //0.95:
                  z=1.96;
                  break;
            case 2: //0.99:
                  z=2.576;
                  break;
            case 3: //0.999:
                  z=3.29;
                  break;
            case 4: //0.9999:
                  z=3.89;
                  break;
            default:
                  z=1.96;
      }
      int min_dist=1;
      float add_dist_upper_bound=max_umi_len*error_rate+z*sqrt(error_rate*(1-error_rate)*max_umi_len ) ;
      int upper_ceil=ceil(add_dist_upper_bound);
      int upper_floor=floor(add_dist_upper_bound);
      int upper_round=round(add_dist_upper_bound);
      if (verbose)
	      cout<<"add_dist_upper_bound="<<add_dist_upper_bound<<" ceil="<<upper_ceil<<" floor="<<upper_floor<<" round="<<upper_round<<endl;
      int dist_upper_bound;
      if  (max_umi_len<=15)
	      dist_upper_bound=min_dist+upper_floor;
      else
	      dist_upper_bound=min_dist+upper_round;
      if (verbose)
	      cout<<"Dist="<<dist_upper_bound<<endl;
      return dist_upper_bound;
}

int count_umi(const string filename){
	int count=0;
	ifstream in_file(filename, ios::in )  ;
	if(in_file.fail()){
            cerr<<"input file for UMI counting: "<<filename <<" not readable"<<endl;
            exit(1);
	}
	string primer_id;
	string umi;
	int read_count;
	while (in_file>>primer_id>>umi>>read_count) {
		if (read_count > 0 )
			count++;
	}
	return count;
}

void read_umi_data(const string  filename, vector<int> & umi_data, vector<int> & read_replicated_data){
      ifstream in_file(filename, ios::in )  ;
      if(in_file.fail()){
            cerr<<"input file for data modeling: "<<filename <<" not readable"<<endl;
            exit(1);
      }
      string primer_id;
      string umi;
      int read_count;
      while (in_file>>primer_id>>umi>>read_count) {
            if (read_count > 0 ) {
                read_replicated_data.insert(read_replicated_data.end(), read_count, read_count);
                umi_data.push_back(read_count);
            }
      }
      if (umi_data.empty()){
	    cerr<<filename<< " : input data file has no record"<<endl;
            exit(1);
      }
}

void do_kp( const string updated_count_file, int  min_read_founder, int kp_estimated_molecule, int kp_angle, int median_rpu,  ofstream & e_out_file ){
	int after_rpucut_molecule=0;
        fit_knee_plot ( updated_count_file,  min_read_founder,  kp_estimated_molecule, kp_angle, median_rpu, after_rpucut_molecule);
	if (verbose)
		cout<<"After UMI clustering:"<<"\t"<<"rpu_cutoff using knee plot="<<min_read_founder<<"\t"<<"kp_estimate_molecules="<<kp_estimated_molecule<<endl;
	e_out_file<<"KP_estimate\tON"<<endl;
	e_out_file<<"median_rpu\t"<<median_rpu<<endl;
        e_out_file<<"rpu_cutoff\t"<<min_read_founder<<endl;
	e_out_file<<"estimated_molecules\t"<<after_rpucut_molecule<<endl;
}

double lengthSquare(pair<double,double> X, pair<double,double> Y)
{
    double xDiff = X.first - Y.first;
    double yDiff = X.second - Y.second;
    return xDiff*xDiff + yDiff*yDiff;
}

double get_angle( pair<double,double> A, pair<double,double> B , pair<double,double> C)
{
	// Square of lengths be a2, b2, c2
	double a2 = lengthSquare(B,C);
	double b2 = lengthSquare(A,C);
	double c2 = lengthSquare(A,B);

	// length of sides be a, b, c
	double a = sqrt(a2);
	double b = sqrt(b2);
	double c = sqrt(c2);

	double beta = acos((a2 + c2 - b2)/(2*a*c));
	beta = beta * 180 / PI;
	return (beta);
}

void fit_knee_plot (const string  filename, int & min_read_founder, int & kp_estimated_molecule, int & kp_angle, int & median_rpu, int & after_rpucut_molecule ){
	vector<int>  umi_data;
	set<int> umi_data_set;
	vector<int> read_replicated_data;
        read_umi_data(filename, umi_data, read_replicated_data);
	median_rpu=median(read_replicated_data);
	copy(umi_data.begin(), umi_data.end(), inserter(umi_data_set, umi_data_set.end()));
	if (umi_data.size() <=2){
		kp_estimated_molecule=umi_data.size();
		kp_angle=0;
		after_rpucut_molecule=kp_estimated_molecule;
		return;
	}
	vector<int> umi_cumsum(umi_data.size());
	std::partial_sum (umi_data.begin(), umi_data.end(), umi_cumsum.begin());
	int total_reads=umi_cumsum[umi_cumsum.size()-1] ;
	int total_umi=umi_cumsum.size() ;
	double rpu=total_reads*1.0/total_umi;
	pair<double,double> A = make_pair(1.0, 1.0);
	pair<double, double> C = make_pair(0.0, 0.0);
	vector<int> smooth_angles;
	for (int i=0; i<umi_cumsum.size()-2; i++){
		pair<double,double> B = make_pair(i*1.0/total_umi, umi_cumsum[i]*1.0/total_reads);
		double angle=get_angle( A, B, C);
		int smooth_angle=floor(angle+0.5);
		smooth_angles.push_back(smooth_angle);
	}
	int min = *min_element(smooth_angles.begin(), smooth_angles.end());
	kp_angle=min+1; //to get a bit more Molecule, so we select 1 degree above inflect point as cutoff
	std::vector<int> indices = findItems(smooth_angles, kp_angle);
	if (indices.size()==0)
		indices = findItems(smooth_angles, min);
	int last_index=indices[indices.size()-1];
	min_read_founder=umi_data[last_index];
	std::vector<int> min_read_founder_indices = findItems(umi_data, min_read_founder);
	kp_estimated_molecule=last_index+1;
	if (verbose)
		cout<<"min angle="<<min<<" min_read_founder="<<min_read_founder<<" min_read_founder_indices first="<<min_read_founder_indices[0]<<" min_read_founder_indices last="<<min_read_founder_indices[min_read_founder_indices.size()-1]<<" the index of smallest angle is "<<last_index << endl;
	int min_read_founder_indices_median=(min_read_founder_indices[0] + min_read_founder_indices[min_read_founder_indices.size()-1] ) / 2;
	if (last_index < min_read_founder_indices_median && min_read_founder==1 ){ //when the rpu cuttoff go to 1, we should check whether we should keep 1 as cutoff or move to higher rpu cutoff
		std::vector<int> min_read_founder_set_indices = findItems(umi_data_set, min_read_founder);
		if (min_read_founder_set_indices[0] < umi_data_set.size()){
			auto it = next(umi_data_set.begin(),min_read_founder_set_indices[0]+1 );
			min_read_founder=*it;
			if (verbose)
				cout<<"we should make min_read_founder ="<<*it<<endl;
		}
	}
	//when calcualte angle, the last UMI will not be included, so if the minimal is the second to last has same rpu as the last, we will assign the total mole as all UMIs
	if (last_index==total_umi-2 && umi_data[total_umi-2] == umi_data[total_umi-1])
		kp_estimated_molecule=total_umi;
	sort (umi_data.begin(), umi_data.end());
        auto lower_bound_it = lower_bound(umi_data.begin(), umi_data.end(), min_read_founder );
        after_rpucut_molecule=umi_data.size()-(lower_bound_it - umi_data.begin());

}

void do_nb(const string  updated_count_file, float nb_lowertail_p,  int min_read_founder,  int nb_estimated_molecule, int median_rpu, ofstream & e_out_file){
        fit_nb_model(updated_count_file, nb_lowertail_p,  min_read_founder,  nb_estimated_molecule, median_rpu);
	if (verbose)
        	cout<<"After UMI clustering:"<<"\t"<<"rpu_cutoff using NB model="<<min_read_founder<<"\t"<<"nb_estimate_molecules="<<nb_estimated_molecule<<"\t"<<"median_rpu="<<median_rpu<<endl;
	e_out_file<<"NB_estimate\tON"<<endl;
	e_out_file<<"median_rpu\t"<<median_rpu<<endl;
        e_out_file<<"rpu_cutoff\t"<<min_read_founder<<endl;
        e_out_file<<"estimated_molecules\t"<<nb_estimated_molecule<<endl;

}

void fit_nb_model( const string  filename, float  p,  int & min_read_founder, int &  nb_estimated_molecule, int & median_rpu){
      vector<int>  umi_data;
      vector<int> read_replicated_data;
      read_umi_data(filename, umi_data,read_replicated_data );
      if (umi_data.size() <= 1 ){
                min_read_founder=umi_data[0];
                nb_estimated_molecule=umi_data.size();
                cout<<"UMI input file only has one UMI "<<endl;
                return;
      }
      median_rpu=median(read_replicated_data);
      int median_estimated_molecule=read_replicated_data.size()/median_rpu;
      if (verbose)
	      cout<<"median_estimated_molecule="<<median_estimated_molecule<<endl;
      if (var(umi_data ) == 0){ //this is the case that only one rpu present
                min_read_founder=umi_data[umi_data.size()-1];
                nb_estimated_molecule=umi_data.size();
                if (verbose)
                        cout<<"Data to fit has var = 0"<<"\nMedian="<<median_rpu<<" min_read_founder="<<min_read_founder<<" nb_estimated_molecule="<<nb_estimated_molecule<<endl;
                return;
      }
      long double nb_p=mean(umi_data)/var(umi_data );
      long double nb_r=pow(mean(umi_data), 2)/(var(umi_data)-mean(umi_data));
      if (verbose)
                cout<<"\ndata_size="<<umi_data.size()<<"\nNB_p="<<nb_p<<"\nNB_r="<<nb_r<<"\nMedian="<<median_rpu<<"\nMean="<<mean(umi_data)<<"\nVar="<<var(umi_data )<<endl;

      if (nb_p<=0 || nb_p >=1 || nb_r <=0 ){ //in this case, nb fitting is bad and will cause esimate of nb_p or nb_r very inaccurate, (found in some high input samples!) so just give UMI clustering results
                min_read_founder=umi_data[umi_data.size()-1];
                nb_estimated_molecule=umi_data.size();
                return;
      }
      //int lower_nb = boost::math::quantile( boost::math::negative_binomial(nb_r, nb_p), p/2) ;
      int lower_nb = boost::math::quantile( boost::math::negative_binomial(nb_r, nb_p), p/2) ;
      if (lower_nb==0)
                lower_nb=1; //observed UMI has at least one read
      sort (umi_data.begin(), umi_data.end());
      auto lower_bound_it = lower_bound(umi_data.begin(), umi_data.end(), lower_nb);
      vector<int> umi_filtered_data(lower_bound_it, umi_data.end()) ;
      min_read_founder=lower_nb;
      nb_estimated_molecule=umi_filtered_data.size();
}

bool align_umi(const string& umi, const string& f, int max_dist, int& endpos, int& dist){
      EdlibAlignMode mode = EDLIB_MODE_HW;
      EdlibAlignTask task = EDLIB_TASK_DISTANCE;
      EdlibAlignConfig edlib_AlignConfig = {max_dist, mode, task, NULL, 0};
      EdlibAlignResult result = edlibAlign(umi.c_str(), umi.length(), f.c_str(), f.length(), edlib_AlignConfig);

      if (result.status == EDLIB_STATUS_OK ) {
	    if (result.editDistance <= max_dist && result.editDistance >=0) {
		  endpos=result.endLocations[0];
		  edlibFreeAlignResult(result);
		  dist=result.editDistance;
		  return true;
	    }
      }
      edlibFreeAlignResult(result);
      return false;
}

string producer(const vector<UMI_item> &umi_pool_subset, const string primer_id, const unsigned max_dist, const unsigned max_umi_len)
{
      producer_done.store(false);

      int founder_clustered=0;
      string buf;
      buf.reserve(umi_pool_subset.size() * 50);
      string padded_cached = padding;
      int cached_count = 0;
      const int K = 6;
      const int slot = (int)(max_umi_len + N_num);
      // Direct-indexed k-mer table (same as consumer): faster than unordered_map.
      array<vector<int>, 4096> prod_kmer_tbl;
      vector<bool> prod_in_cand;   // dedup bitset sized to cached_count
      vector<int> prod_cands;
      prod_cands.reserve(32);

      for (auto const & u: umi_pool_subset){
	    int offset = u.founder_offset;
	    const string& umi = u.UMI_seq;
	    string line;
	    line.reserve(primer_id.size() + 2 + umi.size() + 2 + umi.size() + 1);
	    line = primer_id; line += '\t'; line += umi; line += '\t'; line += umi; line += '\n';
	    if (founder_added>0){
		  int curr_size = founders.myvector.size();
		  if (curr_size > cached_count) {
			if ((int)prod_in_cand.size() < curr_size)
			      prod_in_cand.resize(curr_size, false);
			for (int i = cached_count; i < curr_size; ++i) {
			      const string& f = founders.myvector[i];
			      padded_cached.append(f);
			      padded_cached.append(max_umi_len + N_num - f.size(), 'N');
			      for (int p = 0; p + K <= (int)f.size(); ++p) {
				    int km = encode_kmer(f.data() + p);
				    if (km >= 0) prod_kmer_tbl[km].push_back(i);
			      }
			}
			cached_count = curr_size;
		  }
		  // K-mer pre-filter: find candidates in [offset, cached_count)
		  prod_cands.clear();
		  for (int p = 0; p + K <= (int)umi.size(); ++p) {
			int km = encode_kmer(umi.data() + p);
			if (km < 0) continue;
			for (int idx : prod_kmer_tbl[km]) {
			      if (idx >= offset && !prod_in_cand[idx]) {
				    prod_in_cand[idx] = true;
				    prod_cands.push_back(idx);
			      }
			}
		  }
		  for (int idx : prod_cands) prod_in_cand[idx] = false;
		  sort(prod_cands.begin(), prod_cands.end());

		  int endpos, dist;
		  bool matched = false;
		  int matched_ind = -1;
		  if (!prod_cands.empty()) {
			string small_target = padding;
			small_target.reserve(N_num + prod_cands.size() * slot);
			for (int idx : prod_cands) {
			      size_t start = N_num + (size_t)idx * slot;
			      small_target.append(padded_cached, start, slot);
			}
			if (align_umi(u.padded_umi, small_target, max_dist, endpos, dist)) {
			      matched = true;
			      matched_ind = prod_cands[(endpos - N_num) / slot];
			}
		  }
		  if (matched) {
			founder_clustered++;
			const string& clustered_founder = founders.myvector[matched_ind];
			const string& chosen = (u.founder_temp_found && dist >= u.founder_temp_dist) ? u.founder_temp : clustered_founder;
			line.clear();
			line.reserve(primer_id.size() + 2 + umi.size() + 2 + chosen.size() + 1);
			line = primer_id; line += '\t'; line += umi; line += '\t'; line += chosen; line += '\n';
		  }
		  else if (u.founder_temp_found ){
			line.clear();
			line.reserve(primer_id.size() + 2 + umi.size() + 2 + u.founder_temp.size() + 1);
			line = primer_id; line += '\t'; line += umi; line += '\t'; line += u.founder_temp; line += '\n';
		  }
		  else{
			founders.guard.lock();
			founders.myvector.push_back(umi);
			founders.guard.unlock();
			founder_added++;
		  }
	    }
	    else{
		  founders.guard.lock();
		  founders.myvector.push_back(umi);
		  founders.guard.unlock();
		  founder_added=1;
	    }
	    buf += line;
      }

      founders.guard.lock();
      founders.size_last_cycle=founders.myvector.size();
      founders.guard.unlock();
      producer_done.store(true);
      return buf;
}

pair<string, vector<UMI_item>> consumer(const int worker_index, bool first_founder_mode, const string  primer_id, const unsigned  max_dist, const unsigned max_umi_len, vector<UMI_item> umi_pool_subset, const int t_umi_size)
{
      const int K = 6;
      const int slot = (int)(max_umi_len + N_num);
      int stop_search_dist=1;
      if (first_founder_mode)
	    stop_search_dist=max_dist;
      if (umi_pool_subset.size()==0){
	    return {"", umi_pool_subset};
      }
      string buf;
      buf.reserve(umi_pool_subset.size() * 50);
      vector<UMI_item> updated_subset;
      int founders_start_offset=0;

      // Hoist per-cycle allocations outside the loop so memory is reused across cycles.
      // kmer_tbl: fixed 4096-slot direct array (12-bit encoded k-mer → founder indices).
      // Only touched slots are cleared via kmer_used, avoiding full array reset each cycle.
      array<vector<int>, 4096> kmer_tbl;
      vector<int> kmer_used;
      kmer_used.reserve(512);
      vector<bool> in_cand;    // dedup bitset; grown as needed, reset per UMI
      vector<int> candidates;
      candidates.reserve(32);
      string small_target;
      small_target.reserve(32 * slot);
      string padded_view;
      vector<string> founder_view;

      while (true) {
	    while ((int)founder_added.load() <= founders_start_offset && !producer_done.load())
		  this_thread::sleep_for(chrono::microseconds(10));
	    int curr_founders_size;
	    {
		  shared_lock<shared_mutex> lk(founders.guard);
		  curr_founders_size = (int)founders.myvector.size();
		  if (curr_founders_size > founders_start_offset)
			founder_view.assign(founders.myvector.begin() + founders_start_offset,
					   founders.myvector.end());
	    }

	    if (curr_founders_size <= founders_start_offset) {
		  return {buf, umi_pool_subset};
	    }

	    int step_size = (int)founder_view.size();

	    // Build k-mer index and pre-padded founder string for this cycle.
	    // Pigeonhole: any two 18-mers with Levenshtein dist ≤ 2 share ≥1 6-mer.
	    kmer_used.clear();
	    padded_view.clear();
	    padded_view = padding;
	    for (int i = 0; i < step_size; ++i) {
		  const string& f = founder_view[i];
		  padded_view.append(f);
		  padded_view.append(max_umi_len + N_num - f.size(), 'N');
		  for (int p = 0; p + K <= (int)f.size(); ++p) {
			int km = encode_kmer(f.data() + p);
			if (km < 0) continue;
			if (kmer_tbl[km].empty()) kmer_used.push_back(km);
			kmer_tbl[km].push_back(i);
		  }
	    }
	    // Grow in_cand if this cycle has more founders than any previous cycle.
	    if ((int)in_cand.size() < step_size) in_cand.resize(step_size, false);

	    for (int j = 0; j < (int)umi_pool_subset.size(); ++j) {
		  const UMI_item& one_umi_ref = umi_pool_subset[j];
		  const string& umi = one_umi_ref.UMI_seq;
		  int offset = one_umi_ref.founder_offset;
		  int relative_start = offset - founders_start_offset;
		  if (relative_start >= step_size) {
			updated_subset.push_back(one_umi_ref);
			continue;
		  }
		  if (relative_start < 0) relative_start = 0;

		  // K-mer pre-filter: collect candidates in [relative_start, step_size)
		  candidates.clear();
		  for (int p = 0; p + K <= (int)umi.size(); ++p) {
			int km = encode_kmer(umi.data() + p);
			if (km < 0) continue;
			for (int idx : kmer_tbl[km]) {
			      if (idx >= relative_start && !in_cand[idx]) {
				    in_cand[idx] = true;
				    candidates.push_back(idx);
			      }
			}
		  }
		  for (int idx : candidates) in_cand[idx] = false; // reset for next UMI
		  sort(candidates.begin(), candidates.end());

		  int endpos, dist;
		  bool matched = false;
		  int matched_ind = -1;
		  if (!candidates.empty()) {
			small_target.clear();
			small_target = padding;
			for (int idx : candidates)
			      small_target.append(padded_view, N_num + (size_t)idx * slot, slot);
			if (align_umi(one_umi_ref.padded_umi, small_target, max_dist, endpos, dist)) {
			      matched = true;
			      matched_ind = candidates[(endpos - N_num) / slot];
			}
		  }

		  if (matched) {
			const string& clustered_founder = founder_view[matched_ind];
			if (dist<=stop_search_dist){
			      buf += primer_id; buf += '\t'; buf += umi; buf += '\t'; buf += clustered_founder; buf += '\n';
			}
			else {
			      UMI_item one_umi = one_umi_ref;
			      if( (!one_umi.founder_temp_found ) || (one_umi.founder_temp_found  && dist<one_umi.founder_temp_dist )  ){
				    one_umi.founder_temp_found=true;
				    one_umi.founder_temp=clustered_founder;
				    one_umi.founder_temp_dist=dist;
			      }
			      one_umi.founder_offset = founders_start_offset + step_size;
			      updated_subset.push_back(std::move(one_umi));
			}
		  }
		  else{
			UMI_item one_umi = one_umi_ref;
			one_umi.founder_offset = founders_start_offset + step_size;
			updated_subset.push_back(std::move(one_umi));
		  }
		  if  (!first_founder_mode && producer_done.load() ){
			updated_subset.insert(updated_subset.end(), umi_pool_subset.begin()+j+1, umi_pool_subset.end());
			umi_pool_subset = updated_subset;
			return {buf, umi_pool_subset};
		  }
	    }
	    umi_pool_subset = updated_subset;
	    updated_subset.clear();
	    if (first_founder_mode)
		  return {buf, umi_pool_subset};
	    founders_start_offset += step_size;
	    // Clear only the k-mer table slots we touched this cycle.
	    for (int km : kmer_used) kmer_tbl[km].clear();
      }
}

vector <int> split_umi_to_threads_on_founder(vector<UMI_item> umi_pool, int threads, int pool_size){
      vector<int> thread_founder_ends;
      if (umi_pool.size()<1000){
	    for (int i=0; i<threads; i++){
		  thread_founder_ends.push_back(umi_pool.size());
	    }
	    return thread_founder_ends;
      }

      ulong total_founder_computation=0;
      int founders_size=founders.size_last_cycle;
      for (auto umi : umi_pool){
	    total_founder_computation+=(founders_size-umi.founder_offset);
      }
      //for the first round:
      if (total_founder_computation==0){
	    for (int i=0; i<threads; i++){
		  thread_founder_ends.push_back(pool_size*(i+1));
	    }
	    return thread_founder_ends;
      }

      ulong founder_computation_per_thread=ceil(total_founder_computation/threads);
      int umi_counter=0;
      int founder_counter=0;
      int ct=0;
      int umis_per_thread=0;
      for (auto umi : umi_pool){
	    founder_counter+=(founders.size_last_cycle-umi.founder_offset) ;
	    umi_counter++;
	    umis_per_thread++;
	    if ( umis_per_thread < 200 )
		 continue;
            if (founder_counter>=founder_computation_per_thread){
		  ct++;
                  thread_founder_ends.push_back(umi_counter);
                  founder_counter=0;
		  umis_per_thread=0;
            }
      }
      for (int i=thread_founder_ends.size(); i<threads; i++){
	    thread_founder_ends.push_back(umi_pool.size());
      }
      thread_founder_ends[threads-1]=umi_pool.size(); //make sure the last thread will include last UMI
      return thread_founder_ends;
}

void parallel_processing( vector<UMI_item>& umi_pool,  ofstream& out, UMI_clustering_parameters parameters)
{
      int max_dist = parameters.max_dist;
      int num_worker_threads=parameters.thread;
      int pool_size=parameters.pool_size;
      int max_umi_len=parameters.max_umi_len;
      int min_read_founder = parameters.min_read_founder;
      bool first_founder_mode=parameters.greedy_mode;
      string curr_primer_id= parameters.primer_id;

      vector<int> thread_founder_ends=split_umi_to_threads_on_founder(umi_pool, parameters.thread, pool_size);

      producer_done.store(false);
      vector<UMI_item> updated_umi_pool;
      map<int, vector<UMI_item>> t_umis;
      map<int, vector<UMI_item>> consumer_results;
      vector<future<pair<string, vector<UMI_item>>>> cl_consumers;
      const int t_umi_size = num_worker_threads - 1;

      vector<UMI_item> umi_pool_subset(umi_pool.begin(), umi_pool.begin()+min(static_cast<unsigned long>(thread_founder_ends[0]) , umi_pool.size()));
      future<string> fp_producer = async(launch::async, [umi_pool_subset, curr_primer_id, max_dist, max_umi_len]{
	    return producer(cref(umi_pool_subset), curr_primer_id, max_dist, max_umi_len);
      });
      for (int i = 1; i < num_worker_threads; i++) {
	    int first = min(static_cast<unsigned long>(thread_founder_ends[i-1]) , umi_pool.size());
	    int last = min(static_cast<unsigned long>(thread_founder_ends[i]) , umi_pool.size());
	    vector<UMI_item> one_d_umis;
	    if (first == umi_pool.size()) {
		  one_d_umis={};
	    }
	    else {
		  one_d_umis={umi_pool.begin()+first, umi_pool.begin()+last};
	    }
	    t_umis[i]=one_d_umis;
      }

      // Track which t_umis index each cl_consumers entry corresponds to so results
      // can be stored back under the correct key regardless of empty-slot skips.
      vector<int> consumer_thread_ids;
      for (int i = 1; i < num_worker_threads; i++) {
	    if (t_umis[i].size()>0) {
		  consumer_thread_ids.push_back(i);
		  // Capture only this thread's vector (not the whole map) to avoid N full-map copies.
		  cl_consumers.push_back(async(launch::async, [i, first_founder_mode, curr_primer_id, max_dist, max_umi_len, v=t_umis[i], t_umi_size]{
			return consumer(i, first_founder_mode, curr_primer_id, max_dist, max_umi_len, std::move(v), t_umi_size);
		  }));
	    }
      }
      // Wait for producer; it sets producer_done which unblocks consumers.
      // Write all output sequentially after join — no mutex needed.
      string producer_buf = fp_producer.get();
      if (!producer_buf.empty()) out << producer_buf;
      for (int j = 0; j < (int)cl_consumers.size(); j++) {
	    auto result = cl_consumers[j].get();
	    if (!result.first.empty()) out << result.first;
	    consumer_results[consumer_thread_ids[j]] = std::move(result.second);
      }
      for (int i = 1; i <= (int)t_umis.size(); i++) {
	    updated_umi_pool.insert(updated_umi_pool.end(), consumer_results[i].begin(), consumer_results[i].end());
      }
      umi_pool = updated_umi_pool;
}

// Accumulates output for umi_pool[begin..end) into a returned string; avoids per-line mutex.
string founder_find ( const vector<UMI_item>& umi_pool, int begin, int end, const string primer_id, const unsigned  max_dist, const unsigned max_umi_len, const int pool, const array<vector<int>, 4096>& kmer_idx, const vector<string>& founder_view, const string& padded_founders ){
      const int K = 6;
      const int slot = (int)(max_umi_len + N_num);
      string buf;
      buf.reserve((end - begin) * 50);

      int n_founders = (int)founder_view.size();
      vector<bool> in_cand(n_founders, false);
      vector<int> candidates;
      candidates.reserve(32);
      string small_target;
      small_target.reserve(32 * slot);

      for (int j = begin; j < end; ++j) {
	    const UMI_item& one_umi = umi_pool[j];
	    const string& umi = one_umi.UMI_seq;

	    candidates.clear();
	    for (int p = 0; p + K <= (int)umi.size(); ++p) {
		  int km = encode_kmer(umi.data() + p);
		  if (km < 0) continue;
		  for (int idx : kmer_idx[km]) {
			if (!in_cand[idx]) {
			      in_cand[idx] = true;
			      candidates.push_back(idx);
			}
		  }
	    }
	    for (int idx : candidates) in_cand[idx] = false;
	    sort(candidates.begin(), candidates.end());

	    if (!candidates.empty()) {
		  small_target.clear();
		  small_target = padding;
		  for (int idx : candidates)
			small_target.append(padded_founders, N_num + (size_t)idx * slot, slot);
		  int endpos, dist;
		  if (align_umi(one_umi.padded_umi, small_target, max_dist, endpos, dist)) {
			int ind = candidates[(endpos - N_num) / slot];
			const string& clustered_founder = founder_view[ind];
			buf += primer_id; buf += '\t'; buf += umi; buf += '\t'; buf += clustered_founder; buf += '\n';
		  } else {
			buf += primer_id; buf += '\t'; buf += umi; buf += '\t'; buf += umi; buf += '\n';
		  }
	    } else {
		  buf += primer_id; buf += '\t'; buf += umi; buf += '\t'; buf += umi; buf += '\n';
	    }
      }
      return buf;
}

void parallel_founder_find ( vector<UMI_item> low_reads_umi_pool,  ofstream& out, const unsigned num_worker_threads, const string& curr_primer_id, const unsigned  max_dist, const unsigned max_umi_len  ){
      // Build k-mer index, founder_view, and padded string once; share read-only across all threads.
      const int K = 6;
      const int slot = (int)(max_umi_len + N_num);
      vector<string> shared_founder_view(founders.myvector.begin(), founders.myvector.end());
      int n_founders = (int)shared_founder_view.size();
      array<vector<int>, 4096> shared_kmer_idx;
      string shared_padded_founders = padding;
      shared_padded_founders.reserve(N_num + (size_t)n_founders * slot);
      for (int i = 0; i < n_founders; ++i) {
	    const string& f = shared_founder_view[i];
	    shared_padded_founders.append(f);
	    shared_padded_founders.append(max_umi_len + N_num - f.size(), 'N');
	    for (int p = 0; p + K <= (int)f.size(); ++p) {
		  int km = encode_kmer(f.data() + p);
		  if (km >= 0) shared_kmer_idx[km].push_back(i);
	    }
      }

      // Divide by index range — avoids deep-copying 66k UMI_item objects (each with heap strings).
      int n = (int)low_reads_umi_pool.size();
      int chunk = (n + (int)num_worker_threads - 1) / (int)num_worker_threads;
      vector<future<string>> futures;
      for (int i = 0; i < (int)num_worker_threads; i++) {
	    int s = i * chunk;
	    int e = min(n, s + chunk);
	    futures.push_back(async(launch::async, [s, e, &low_reads_umi_pool, curr_primer_id, max_dist, max_umi_len, i,
	                                            &shared_kmer_idx, &shared_founder_view, &shared_padded_founders]() {
		  return founder_find(low_reads_umi_pool, s, e, curr_primer_id, max_dist, max_umi_len, i,
		                      shared_kmer_idx, shared_founder_view, shared_padded_founders);
	    }));
      }
      for (auto& f : futures) {
	    string buf = f.get();
	    if (!buf.empty()) out << buf;
      }
}

void clustering_umis(const string in_filename, const string out_filename,  UMI_clustering_parameters parameters )
{
      int max_dist = parameters.max_dist;
      int num_worker_threads=parameters.thread;
      int pool_size=parameters.pool_size;
      int max_umi_len=parameters.max_umi_len;
      int min_read_founder = parameters.min_read_founder;
      int total_umis_in_a_run = num_worker_threads * pool_size;
      vector<UMI_item> umi_pool;
      vector<UMI_item> low_reads_umi_pool;
      ifstream in_file(in_filename, ios::in )  ;
      if(in_file.fail()){
	    cerr<<"input file: "<<in_filename <<" not readable"<<endl;
	    exit(1);
      }
      ofstream out_file(out_filename, ios::out );
      if( out_file.fail()){
	    cerr<<"output file: "<<out_filename <<" not writable"<<endl;
	    exit(1);
      }

      string line;
      string primer_id;
      string umi;
      int read_count;
      UMI_item one_umi;
      string curr_primer_id;
      int round=0;
      int lines=0;

      while (getline(in_file, line)) {
	    istringstream ss(line);
	    lines++;
	    ss >> primer_id >> umi >> read_count;
	    if (umi.length()>max_umi_len)
		  umi.resize(max_umi_len);
	    one_umi.UMI_seq = umi;
	    one_umi.padded_umi = padding + umi + padding;
	    one_umi.founder_offset = 0;
	    if (read_count < min_read_founder ){ //reach the line where the read count for umi less than the minimal requirement of founder
		  while (!umi_pool.empty()) { //finish current umi_pool processing first
			parallel_processing( umi_pool,  out_file, parameters );
		  }
		  if (primer_id != curr_primer_id )   { //A new primer and this primer the first UMI < min_read_founder
			if (!low_reads_umi_pool.empty())
				parallel_founder_find(low_reads_umi_pool,  out_file, num_worker_threads, curr_primer_id,  max_dist,max_umi_len  );
			low_reads_umi_pool.clear();
		  }
		  low_reads_umi_pool.push_back(one_umi);
		  curr_primer_id=primer_id;
		  parameters.primer_id=curr_primer_id;
		  continue;
	    }

	    if (primer_id != curr_primer_id ) { //a new primer
		  if ( verbose)
			cout<<"A new group or primer:"<<primer_id<<" umi:"<<umi<<"\n";
		  while (!umi_pool.empty()) {
			parallel_processing( umi_pool,  out_file, parameters );
                  }
		  if (!low_reads_umi_pool.empty()){
			parallel_founder_find(low_reads_umi_pool,  out_file, num_worker_threads, curr_primer_id,  max_dist,max_umi_len  );
			low_reads_umi_pool.clear();
		  }
		  founders.myvector.clear();
		  founders.size_last_cycle=0;
		  round=0;
	    }
	    curr_primer_id=primer_id;
	    parameters.primer_id=curr_primer_id;
	    umi_pool.push_back(one_umi);
	    if ( umi_pool.size() == total_umis_in_a_run){
		  round++;
		  parallel_processing( umi_pool,  out_file, parameters);
	    }
      }
      //read last line of input file
      while (!umi_pool.empty()) {
	    parallel_processing(umi_pool, out_file, parameters);
      }
      if (!low_reads_umi_pool.empty()){
	    parallel_founder_find(low_reads_umi_pool,  out_file, num_worker_threads, curr_primer_id,  max_dist,max_umi_len  );
      }
      cout <<"All done!" << endl;
}

void update_umi_reads_count(const string updated_count_filename, const string in_filename, const string out_filename){
	unordered_map<string, int> umi_count;
	ifstream in_file(in_filename, ios::in )  ;
	if(in_file.fail()){
		cerr<<"input file: "<<in_filename <<" not readable"<<endl;
		exit(1);
	}
	string primer_id;
	string umi;
	int read_count;
	while (in_file>>primer_id>>umi>>read_count){
		umi_count[primer_id+"\t"+umi]=read_count;
	}
	ifstream cluster_file(out_filename, ios::in )  ;
	if(cluster_file.fail()){
                cerr<<"cluster file: "<<out_filename <<" not readable"<<endl;
                exit(1);
        }
	string child_umi;
	string founder_umi;
	string info;
	while (cluster_file>>primer_id>>child_umi>>founder_umi) {
		if  (child_umi.compare(founder_umi) != 0 ){
			umi_count[primer_id+"\t"+founder_umi]+=umi_count[primer_id+"\t"+child_umi];
			umi_count[primer_id+"\t"+child_umi]=0;
		}

        }
	multimap<int,string,greater<int>> count_umi = invertMap(umi_count);
	ofstream out_file(updated_count_filename, ios::out );
	for (auto it = count_umi.begin(); it != count_umi.end(); ++it)
		out_file<<it->second<<"\t"<<it->first<<endl;
}
