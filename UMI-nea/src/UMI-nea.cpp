#include "UMI-nea.h"

//global variable
mutex mut;
atomic_int founder_added{0};
// Removed global span; poll loop in parallel_processing now uses non-blocking wait_for(0) + 10µs sleep.
guardedvector founders;
bool producer_done=false;
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

void shared_writer(ofstream& out, const vector<string>& lines)
{
      lock_guard<mutex> lock(mut);
      for (const string& line : lines)
	    out << line;
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

bool producer(const vector<UMI_item> &umi_pool_subset, ofstream& out, const string primer_id, const unsigned max_dist, const unsigned max_umi_len)
{
      producer_done=false;
      int write_every=ceil(umi_pool_subset.size()/10)>10?ceil(umi_pool_subset.size()/10):10;

      int founder_clustered=0;
      vector <string> lines;
      // Incrementally extend cached padded founder string (starting from founder 0) to avoid
      // O(N*F) stringstream rebuild. For UMIs with non-zero founder_offset, a substring of
      // padded_cached starting at that offset is used so the index calculation stays correct.
      string padded_cached = padding;
      int cached_count = 0;
      string target_scratch;  // reused buffer for UMIs with founder_offset > 0
      for (auto const & u: umi_pool_subset){
	    int offset = u.founder_offset;
	    const string& umi = u.UMI_seq;
	    //new founder being found:
	    string line;
	    line.reserve(primer_id.size() + 2 + umi.size() + 2 + umi.size() + 1);
	    line = primer_id; line += '\t'; line += umi; line += '\t'; line += umi; line += '\n';
	    bool clustered = false;
	    if (founder_added>0){
		  int curr_size = founders.myvector.size();
		  if (curr_size > cached_count) {
			for (int i = cached_count; i < curr_size; ++i) {
			      const string& f = founders.myvector[i];
			      padded_cached.append(f);
			      padded_cached.append(max_umi_len + N_num - f.size(), 'N');
			}
			cached_count = curr_size;
		  }
		  // Build target starting from founder[offset]: "NNN" + padded_cached[N_num + offset*(slot)]
		  // This keeps the index formula (endpos-N_num)/(slot)+offset correct.
		  const string* target_ptr;
		  if (offset == 0) {
			target_ptr = &padded_cached;
		  } else {
			size_t byte_start = N_num + (size_t)offset * (max_umi_len + N_num);
			target_scratch.clear();
			if (byte_start < padded_cached.size()) {
			      target_scratch.reserve(N_num + padded_cached.size() - byte_start);
			      target_scratch = padding;
			      target_scratch.append(padded_cached, byte_start, string::npos);
			}
			target_ptr = &target_scratch;
		  }
		  int endpos;
		  int dist;
		  if (!target_ptr->empty() && align_umi(u.padded_umi, *target_ptr, max_dist, endpos, dist )) {
			founder_clustered++;
			int ind=(endpos-N_num)/(max_umi_len+N_num)+offset;
			const string& clustered_founder=founders.myvector[ind];
			const string& chosen = (u.founder_temp_found && dist >= u.founder_temp_dist) ? u.founder_temp : clustered_founder;
			line.clear();
			line.reserve(primer_id.size() + 2 + umi.size() + 2 + chosen.size() + 1);
			line = primer_id; line += '\t'; line += umi; line += '\t'; line += chosen; line += '\n';
		  }
		  else if (u.founder_temp_found ){
			//After checking all UMI in the front line, only option is to choose  temp founder
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
	    //shared_writer(out, line);
	    lines.push_back(line);
	    if (lines.size()>=write_every){
		  shared_writer(out, lines);
		  lines.clear();
	    }

      }
      if (lines.size()>=0){
	    shared_writer(out, lines);
      }

      producer_done=true;
      founders.guard.lock();
      founders.size_last_cycle=founders.myvector.size();
      founders.guard.unlock();
      return true;
}

vector<UMI_item> consumer(const int worker_index, ofstream& out, bool first_founder_mode, const string  primer_id, const unsigned  max_dist, const unsigned max_umi_len, vector<UMI_item> umi_pool_subset, const int t_umi_size)
{

      int stop_search_dist=1;
      if (first_founder_mode)
	    stop_search_dist=max_dist;
      if (umi_pool_subset.size()==0){
	    return umi_pool_subset;
      }
      int write_every=10;
      vector <string> lines;
      vector<UMI_item> updated_subset;
      int founders_start_offset=0;
      int step_size_small=2;
      int step_size_big;
      if (first_founder_mode)
	    step_size_big=0;
      else
	    step_size_big=t_umi_size/10;
      int step_size;
      int last_cycle_founders_size = founders.size_last_cycle;
      while (true) {
	    if (producer_done) {
		  if (lines.size()>0)
			shared_writer(out, lines);
		  return umi_pool_subset;
	    }
	    int curr_founders_size;
	    vector <string> founder_view;
	    founders.guard.lock();
	    curr_founders_size = founders.myvector.size();
	    if (curr_founders_size>founders_start_offset)
		  founder_view.assign(founders.myvector.begin()+founders_start_offset ,founders.myvector.end()  );
	    founders.guard.unlock();

	    if (founders_start_offset==0 ){
		  if (curr_founders_size == 0 && last_cycle_founders_size==0)
			continue;
		  else if (curr_founders_size != 0 && last_cycle_founders_size==0){
			step_size=step_size_big;
			if (curr_founders_size  < step_size ){
			      continue;
			}
		  }
		  else
			step_size=last_cycle_founders_size;
	    }
	    else{
		  if (last_cycle_founders_size==0)
			step_size=step_size_big;
		  else
			step_size=step_size_small;
		  if (curr_founders_size - founders_start_offset < step_size ){
			continue;
		  }
	    }
	    // Build padded founder string once per cycle; all UMIs with byte_start==0 share it.
	    string full_founder_string;
	    full_founder_string.reserve((size_t)step_size * (max_umi_len + N_num));
	    for (int i = 0; i < step_size; ++i) {
		  const string& fi = founder_view[i];
		  full_founder_string.append(fi);
		  full_founder_string.append(max_umi_len + N_num - fi.size(), 'N');
	    }
	    string padded_full;
	    padded_full.reserve(N_num + full_founder_string.size());
	    padded_full = padding;
	    padded_full += full_founder_string;
	    string padded_f_scratch;  // reused buffer for UMIs with non-zero byte_start
	    for (int j = 0; j < umi_pool_subset.size(); ++j) {
		  const UMI_item& one_umi_ref = umi_pool_subset[j];
		  const string& umi = one_umi_ref.UMI_seq;
		  int offset = one_umi_ref.founder_offset;
		  int relative_start = offset - founders_start_offset;
		  size_t byte_start = (size_t)relative_start * (max_umi_len + N_num);
		  // Common case: all UMIs in the first cycle have byte_start==0 and share padded_full.
		  // Minority case: reuse padded_f_scratch to avoid repeated malloc/free per UMI.
		  const string* target_ptr;
		  if (byte_start == 0) {
			target_ptr = &padded_full;
		  } else {
			padded_f_scratch.clear();
			padded_f_scratch.reserve(N_num + full_founder_string.size() - byte_start);
			padded_f_scratch = padding;
			padded_f_scratch.append(full_founder_string, byte_start, string::npos);
			target_ptr = &padded_f_scratch;
		  }
		  int endpos;
		  int dist;
		  if (align_umi(one_umi_ref.padded_umi, *target_ptr, max_dist, endpos, dist)) {
			int ind=(endpos-N_num)/(max_umi_len+N_num) + relative_start;
			const string& clustered_founder=founder_view[ind];

			if (dist<=stop_search_dist){
			      //founder being found in consumer pool as the distance to founder <= 1 or if first founder mode is on and dist is <= cutoff
			      string line;
			      line.reserve(primer_id.size() + 2 + umi.size() + 2 + clustered_founder.size() + 1);
			      line = primer_id; line += '\t'; line += umi; line += '\t'; line += clustered_founder; line += '\n';
			      lines.push_back(std::move(line));
			      if (lines.size()>=write_every){
				    shared_writer(out, lines);
				    lines.clear();
			      }
			}
			else {
			      UMI_item one_umi = one_umi_ref;
			      if( (!one_umi.founder_temp_found ) || (one_umi.founder_temp_found  && dist<one_umi.founder_temp_dist )  ){
				    one_umi.founder_temp_found=true;
				    one_umi.founder_temp=clustered_founder;
				    one_umi.founder_temp_dist=dist;
			      }
			      if (step_size!=last_cycle_founders_size)
				    one_umi.founder_offset += step_size;
			      else
				    one_umi.founder_offset = last_cycle_founders_size;
			      updated_subset.push_back(std::move(one_umi));
			}
		  }
		  else{
			UMI_item one_umi = one_umi_ref;
			if (step_size!=last_cycle_founders_size)
			      one_umi.founder_offset += step_size;
			else
			      one_umi.founder_offset = last_cycle_founders_size;
			updated_subset.push_back(std::move(one_umi));
		  }
		  if  (!first_founder_mode && producer_done ){ //when -d is not set and producer is done, we can exit loop now and will not affect resluts reproduciblity
			updated_subset.insert(updated_subset.end(), umi_pool_subset.begin()+j+1, umi_pool_subset.end());
			umi_pool_subset = updated_subset;
			if (lines.size()>0)
			      shared_writer(out, lines);
			return umi_pool_subset;
		  }
	    }
	    umi_pool_subset = updated_subset;
	    updated_subset.clear();
	    if (lines.size()>0){
		  shared_writer(out, lines);
		  lines.clear();
	    }
	    if (first_founder_mode  )
		  return umi_pool_subset;
	    founders_start_offset += step_size;
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

      producer_done=false;
      vector<UMI_item> updated_umi_pool;
      map<int, vector<UMI_item>> t_umis;
      map<int, vector<UMI_item>> consumer_results;
      vector<future<vector<UMI_item>>> cl_consumers;
      const int t_umi_size = num_worker_threads - 1;

      vector<UMI_item> umi_pool_subset(umi_pool.begin(), umi_pool.begin()+min(static_cast<unsigned long>(thread_founder_ends[0]) , umi_pool.size()));
      shared_future<bool>  fp_producer = async(launch::async, [umi_pool_subset, &out, curr_primer_id,max_dist, max_umi_len]{return producer( cref(umi_pool_subset), out, curr_primer_id, max_dist, max_umi_len);});
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

      for (int i = 1; i < num_worker_threads; i++) {
	    if (t_umis[i].size()>0)
		  // Capture only this thread's vector (not the whole map) to avoid N full-map copies.
		  cl_consumers.push_back(async(launch::async, [i,&out,first_founder_mode, curr_primer_id,max_dist,max_umi_len, v=t_umis[i], t_umi_size]{return consumer(i, out, first_founder_mode, curr_primer_id, max_dist, max_umi_len, std::move(v), t_umi_size); } ) ) ;

      }
      vector<bool> thread_done_bits(cl_consumers.size()+1, false);
      auto thread_done = count(thread_done_bits.begin(), thread_done_bits.end(), true);
      while(thread_done<cl_consumers.size()+1) {
	    if (producer_done){
		  thread_done_bits[0]=true;
	    }
	    for (int i = 1; i <= (int)cl_consumers.size(); i++) {
		  if ( thread_done_bits[i] ||  t_umis[i].size()==0){
			thread_done_bits[i]=true;
			continue;
		  }
		  // Non-blocking check; avoids blocking 1ms × N_threads per poll iteration.
		  else if (cl_consumers[i-1].wait_for(chrono::seconds(0)) == future_status::ready && cl_consumers[i-1].valid()) {
			thread_done_bits[i]=true;
			consumer_results[i]=cl_consumers[i-1].get();
		  }
	    }
	    thread_done = count(thread_done_bits.begin(), thread_done_bits.end(), true);
	    if (thread_done < (int)cl_consumers.size()+1)
		  this_thread::sleep_for(chrono::microseconds(10));
      }
      for (int i = 1; i <= (int)t_umis.size(); i++) {
	    updated_umi_pool.insert(updated_umi_pool.end(), consumer_results[i].begin(), consumer_results[i].end());
      }
      umi_pool = updated_umi_pool;
}

void founder_find ( vector<UMI_item>& umi_pool, ofstream& out,  string primer_id, const unsigned  max_dist, const unsigned max_umi_len, const int pool ){
      int write_every=100;
      vector <string> lines;

      vector <string> founder_view={founders.myvector.begin(), founders.myvector.end()};
      // Build founder string once — founders don't change during this function.
      string padded_founder = padding;
      padded_founder.reserve(N_num + founder_view.size() * (max_umi_len + N_num));
      for (const string& f : founder_view) {
	    padded_founder.append(f);
	    padded_founder.append(max_umi_len + N_num - f.size(), 'N');
      }
      for (int j = 0; j < umi_pool.size(); ++j) {
	    const UMI_item& one_umi = umi_pool[j];
	    const string& umi = one_umi.UMI_seq;
	    int endpos;
	    int dist;
	    string line;
	    if (align_umi(one_umi.padded_umi, padded_founder, max_dist, endpos, dist)) {
		  int ind=(endpos-N_num)/(max_umi_len+N_num)  ;
		  const string& clustered_founder=founder_view[ind];
		  //founder being found for UMIs not qualified for founder for other UMIs
		  line.reserve(primer_id.size() + 2 + umi.size() + 2 + clustered_founder.size() + 1);
		  line = primer_id; line += '\t'; line += umi; line += '\t'; line += clustered_founder; line += '\n';
	    }
	    else{
		  //no exising founder being found, this UMI is the founder for itself
		  line.reserve(primer_id.size() + 2 + umi.size() + 2 + umi.size() + 1);
		  line = primer_id; line += '\t'; line += umi; line += '\t'; line += umi; line += '\n';
	    }
	    lines.push_back(line);
	    if (lines.size()>=write_every){
		  shared_writer(out, lines);
		  lines.clear();
	    }
      }
      if (lines.size()>0){
	    shared_writer(out, lines);
	    lines.clear();
      }
}

void parallel_founder_find ( vector<UMI_item> low_reads_umi_pool,  ofstream& out, const unsigned num_worker_threads, const string& curr_primer_id, const unsigned  max_dist, const unsigned max_umi_len  ){
      int subpool_size=ceil(low_reads_umi_pool.size()*1.0/num_worker_threads);
      vector<vector<UMI_item>> t_umis(num_worker_threads);
      for (int i = 0; i < (int)num_worker_threads; i++) {
	    int first = min((size_t)(i*subpool_size), low_reads_umi_pool.size());
	    int last = min((size_t)((i+1)*subpool_size), low_reads_umi_pool.size());
	    t_umis[i].assign(low_reads_umi_pool.begin()+first, low_reads_umi_pool.begin()+last);
      }

      vector<thread> threads;
      for (int i = 0; i < (int)num_worker_threads; i++) {
	    threads.push_back(thread(founder_find, ref(t_umis[i]), ref(out), curr_primer_id, max_dist, max_umi_len, i));
      }
      for (auto &th : threads)
	    th.join();
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
	    int process_umi=umi_pool.size();
	    parallel_processing(umi_pool,  out_file, parameters);
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
