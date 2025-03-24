#include "UMI-nea.h"

//global variable
mutex mut;
atomic_int founder_added{0};
chrono::milliseconds span (1);
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

double mad(const vector<int> v){
      vector<double> diff;
      double median_v=median(v);
      for (auto &each: v){
	    diff.push_back(abs(each-median_v) );
      }
      sort(diff.begin(), diff.end());
      return median(diff);
}

int calculate_dist_upper_bound(float error_rate, int max_umi_len){
      int upper_ceil;
      int upper_floor;
      float z;
      int alpha_i;
      if (max_umi_len<=20)
	    alpha_i=1; //alpha=0.95;
      else if (max_umi_len > 20 && max_umi_len <= 30)
	    alpha_i=2; //alpha=0.99;
      else if (max_umi_len > 30 && max_umi_len <= 40)
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
      float error_rate_upper_bound=error_rate+z*sqrt(error_rate*(1-error_rate)/max_umi_len);
      upper_ceil=ceil(max_umi_len*error_rate_upper_bound);
      upper_floor=floor(max_umi_len*error_rate_upper_bound);
      if (max_umi_len<=20)
		if (upper_floor>0)
			return upper_floor;
		else
			return 1;
      else
		if (upper_ceil > 0 )
			return upper_ceil;
		else
			return 1;
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
                vector <int> repeatn = repeat_n (read_count);
                read_replicated_data.insert(read_replicated_data.end(), repeatn.begin(), repeatn.end()) ;
                umi_data.push_back(read_count);
            }
      }
      if (umi_data.empty()){
	    cerr<<filename<< " : input data file has no record"<<endl;
            exit(1);
      }
}

void do_kp( const string updated_count_file, int  min_read_founder, int kp_estimated_molecule, int kp_angle, int median_rpu,  ofstream & e_out_file ){
	cout<<"Before UMI clustering:"<<"\t"<<"rpu_cutoff using knee plot="<<min_read_founder<<"\t"<<"kp_estimate_molecules="<<kp_estimated_molecule<<endl;
	int after_rpucut_molecule=0;
        fit_knee_plot ( updated_count_file,  min_read_founder,  kp_estimated_molecule, kp_angle, median_rpu, after_rpucut_molecule);
	cout<<"After UMI clustering:"<<"\t"<<"rpu_cutoff using knee plot="<<min_read_founder<<"\t"<<"kp_estimate_molecules="<<kp_estimated_molecule<<endl;
	e_out_file<<"KP_estimate\tON"<<endl;
	e_out_file<<"median_rpu\t"<<median_rpu<<endl;
        e_out_file<<"rpu_cutoff\t"<<min_read_founder<<endl;
        e_out_file<<"estimated_molecules\t"<<kp_estimated_molecule<<endl;
	e_out_file<<"after_rpu-cutoff_molecules\t"<<after_rpucut_molecule<<endl;
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

void do_nb(const string  updated_count_file, float nb_lowertail_p, int madfolds, int min_read_founder,  int nb_estimated_molecule, int median_rpu, ofstream & e_out_file){
	 cout<<"Before UMI clustering:"<<"\t"<<"rpu_cutoff using NB model="<<min_read_founder<<"\t"<<"nb_estimate_molecules="<<nb_estimated_molecule<<"\t"<<"median_rpu="<<median_rpu<<endl;
        fit_nb_model(updated_count_file, nb_lowertail_p, madfolds, min_read_founder,  nb_estimated_molecule, median_rpu);
        cout<<"After UMI clustering:"<<"\t"<<"rpu_cutoff using NB model="<<min_read_founder<<"\t"<<"nb_estimate_molecules="<<nb_estimated_molecule<<"\t"<<"median_rpu="<<median_rpu<<endl;
	e_out_file<<"NB_estimate\tON"<<endl;
	e_out_file<<"median_rpu\t"<<median_rpu<<endl;
        e_out_file<<"rpu_cutoff\t"<<min_read_founder<<endl;
        e_out_file<<"estimated_molecules\t"<<nb_estimated_molecule<<endl;
	e_out_file<<"after_rpu-cutoff_molecules\t"<<nb_estimated_molecule<<endl;

}

void fit_nb_model( const string  filename, float  p, int madfolds, int & min_read_founder, int &  nb_estimated_molecule, int & median_rpu){
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
      int mad_rpu=mad(read_replicated_data);
      vector<int> read_replicated_filtered_data=read_replicated_data;
      int lower_bound_rpu;
      if (mad_rpu!=0){ //for high input samples, we saw rpu is below 6 for all UMIs, the MAD will be 0 in those cases.
	     lower_bound_rpu=median_rpu-3*mad_rpu;
             sort (read_replicated_data.begin(), read_replicated_data.end());
             auto lower_bound_it = lower_bound(read_replicated_data.begin(), read_replicated_data.end(), lower_bound_rpu);
	     vector <int> temp(lower_bound_it, read_replicated_data.end() );
	     read_replicated_filtered_data=temp;
      }
      if (var(read_replicated_filtered_data ) == 0){ //this is the case that after removing outliners or the in raw input, only one rpu present
		min_read_founder=umi_data[0];
		nb_estimated_molecule=umi_data.size();
		if (verbose)
			cout<<"Data to fit has var = 0"<<"\nMedian="<<median_rpu<<"\nMean="<<mean(read_replicated_filtered_data)<<"\nVar="<<var(read_replicated_filtered_data )<<"\nMAD="<<mad_rpu<<endl;
		return;
      }
      long double nb_p=mean(read_replicated_filtered_data)/var(read_replicated_filtered_data );
      long double nb_r=pow(mean(read_replicated_filtered_data), 2)/(var(read_replicated_filtered_data)-mean(read_replicated_filtered_data));
      if (verbose)
		cout<<"lower_bound_rpu="<<lower_bound_rpu<<"\ndata_size="<<read_replicated_filtered_data.size()<<"\nNB_p="<<nb_p<<"\nNB_r="<<nb_r<<"\nMedian="<<median_rpu<<"\nMean="<<mean(read_replicated_filtered_data)<<"\nVar="<<var(read_replicated_filtered_data )<<"\nMAD="<<mad_rpu<<endl;
      if (nb_p<=0 || nb_p >=1 || nb_r <=0 ){ //in this case, nb fitting is bad and will cause esimate of nb_p or nb_r not very accurate, (found in some high input samples!) so just give UMI clustering results
		min_read_founder=1;
		nb_estimated_molecule=umi_data.size();
		return;
      }
      int lower_nb = boost::math::quantile( boost::math::negative_binomial(nb_r, nb_p), p/2) ;
      if (lower_nb==0)
		lower_nb=1;
      sort (umi_data.begin(), umi_data.end());
      auto lower_bound_it2 = lower_bound(umi_data.begin(), umi_data.end(), lower_nb);
      vector<int> umi_filtered_data(lower_bound_it2, umi_data.end()) ;
      min_read_founder=lower_nb;
      nb_estimated_molecule=umi_filtered_data.size();
}

void shared_writer(ofstream& out, const vector<string> lines)
{
      lock_guard<mutex> lock(mut);
      for (string line : lines)
	    out << line;
}

bool align_umi( string umi, string f,  int max_dist, int & endpos, int & dist){
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
      for (auto const & u: umi_pool_subset){
	    int offset = u.founder_offset;
	    string umi = u.UMI_seq;
	    //new founder being found:
	    //string line = primer_id + "\t" + umi + "\t" + umi + "\t" + "founder-new\n";
	    string line = primer_id + "\t" + umi + "\t" + umi + "\n";
	    bool clustered = false;
	    if (founder_added>0){
		  stringstream ss;
		  for (int i = offset; i < founders.myvector.size(); ++i) {
			ss<<left<<setfill('N')<<setw(max_umi_len+3)<<founders.myvector[i];
		  }
		  string founder_string=ss.str();
		  int endpos;
		  int dist;
		  if (align_umi(padding+umi+padding, padding+founder_string, max_dist, endpos, dist )) {
			founder_clustered++;
			int ind=(endpos-N_num)/(max_umi_len+N_num)+offset;
			string clustered_founder=founders.myvector[ind];
			if (u.founder_temp_found  && dist<u.founder_temp_dist )
			      //line = primer_id + "\t" + umi + "\t" + clustered_founder + "\tf-ced\n";
			      //Real founder being found, replace temp founder
			      line = primer_id + "\t" + umi + "\t" + clustered_founder + "\n";
			else if (u.founder_temp_found  && dist>=u.founder_temp_dist ){
			      //line = primer_id + "\t" + umi + "\t" + u.founder_temp + "\tf-rec1\n";
			      //After checking all UMI in the front line, temp founder seems to be better
			      line = primer_id + "\t" + umi + "\t" + u.founder_temp + "\n";
			}
			else
			      //line = primer_id + "\t" + umi + "\t" + clustered_founder + "\tf-ced\n";
			      //No temp founder for UMI but founder being found now
			      line = primer_id + "\t" + umi + "\t" + clustered_founder + "\n";
		  }
		  else if (u.founder_temp_found ){
			//line = primer_id + "\t" + umi + "\t" + u.founder_temp + "\tf-rec2\n";
			//After checking all UMI in the front line, only option is to choose  temp founder
			line = primer_id + "\t" + umi + "\t" + u.founder_temp + "\n";
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

vector<UMI_item> consumer(const int worker_index, ofstream& out, bool first_founder_mode, const string  primer_id, const unsigned  max_dist, const unsigned max_umi_len, map<int, vector<UMI_item>>  t_umi)
{

      int stop_search_dist=1;
      if (first_founder_mode)
	    stop_search_dist=max_dist;
      vector<UMI_item> umi_pool_subset=t_umi[worker_index];
      if (umi_pool_subset.size()==0){
	    return umi_pool_subset;
      }
      int write_every=ceil(t_umi.size()/100)>10?ceil(t_umi.size()/100):10;
      vector <string> lines;
      vector<UMI_item> updated_subset;
      int founders_start_offset=0;
      int step_size_small=2;
      int step_size_big;
      if (first_founder_mode)
	    step_size_big=0;
      else
	    step_size_big=t_umi.size()/10;
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
	    for (int j = 0; j < umi_pool_subset.size(); ++j) {
		  UMI_item one_umi = umi_pool_subset[j];
		  string umi = one_umi.UMI_seq;
		  int offset = one_umi.founder_offset;
		  stringstream ss;
		  string founder_string;
		  for (int i = offset-founders_start_offset; i < step_size; ++i) {
			ss<<left<<setfill('N')<<setw(max_umi_len+3)<<founder_view[i];
		  }
		  founder_string=ss.str();
		  int endpos;
		  int dist;
		  if (align_umi( padding+umi+padding, padding+founder_string,  max_dist, endpos, dist)) {
			int ind=(endpos-N_num)/(max_umi_len+N_num) + offset-founders_start_offset  ;
			string clustered_founder=founder_view[ind];

			if (dist<=stop_search_dist){
			      string line = primer_id + "\t" + umi + "\t" + clustered_founder + "\t"+"customer"+to_string(worker_index)+"\n";
			      lines.push_back(line);
			      if (lines.size()>=write_every){
				    shared_writer(out, lines);
				    lines.clear();
			      }
			}
			else {
			      if( (!one_umi.founder_temp_found ) || (one_umi.founder_temp_found  && dist<one_umi.founder_temp_dist )  ){
				    one_umi.founder_temp_found=true;
				    one_umi.founder_temp=clustered_founder;
				    one_umi.founder_temp_dist=dist;
			      }
			      if (step_size!=last_cycle_founders_size)
				    one_umi.founder_offset += step_size;
			      else
				    one_umi.founder_offset = last_cycle_founders_size;
			      updated_subset.push_back(one_umi);
			}
		  }
		  else{
			if (step_size!=last_cycle_founders_size)
			      one_umi.founder_offset += step_size;
			else
			      one_umi.founder_offset = last_cycle_founders_size;
			updated_subset.push_back(one_umi);
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
      bool first_founder_mode=parameters.first_founder_mode;
      string curr_primer_id= parameters.primer_id;

      vector<int> thread_founder_ends=split_umi_to_threads_on_founder(umi_pool, parameters.thread, pool_size);

      producer_done=false;
      vector<UMI_item> updated_umi_pool;
      map<int, vector<UMI_item>> t_umis;
      vector<future<vector<UMI_item>>> cl_consumers;

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
		  cl_consumers.push_back(async(launch::async, [i,&out,first_founder_mode, curr_primer_id,max_dist,max_umi_len, t_umis]{return consumer(i,  out, first_founder_mode, curr_primer_id, max_dist, max_umi_len, t_umis ); } ) ) ;

      }
      vector<bool> thread_done_bits(cl_consumers.size()+1, false);
      auto thread_done = count(thread_done_bits.begin(), thread_done_bits.end(), true);
      while(thread_done<cl_consumers.size()+1) {
	    if (producer_done){
		  thread_done_bits[0]=true;
	    }
	    for (int i = 1; i <= cl_consumers.size(); i++) {
		  if ( thread_done_bits[i] ||  t_umis[i].size()==0){
			thread_done_bits[i]=true;
			continue;
		  }
		  else if (cl_consumers[i-1].wait_for(span) == future_status::ready && cl_consumers[i-1].valid()) {
			thread_done_bits[i]=true;
			t_umis[i]=cl_consumers[i-1].get();
		  }
	    }
	    thread_done = count(thread_done_bits.begin(), thread_done_bits.end(), true);
      }
      for (int i = 1; i <= t_umis.size(); i++) {
	    updated_umi_pool.insert(updated_umi_pool.end(), t_umis[i].begin(), t_umis[i].end());
      }
      umi_pool = updated_umi_pool;
}

void founder_find ( vector<UMI_item>& umi_pool, ofstream& out,  string primer_id, const unsigned  max_dist, const unsigned max_umi_len, const int pool ){
      int write_every=100;
      vector <string> lines;

      vector <string> founder_view={founders.myvector.begin(), founders.myvector.end()};
      for (int j = 0; j < umi_pool.size(); ++j) {
	    UMI_item one_umi = umi_pool[j];
	    string umi = one_umi.UMI_seq;
	    int offset = one_umi.founder_offset;
	    stringstream ss;
	    for (int i = 0; i < founder_view.size(); ++i) {
		  ss<<left<<setfill('N')<<setw(max_umi_len+3)<<founder_view[i];
	    }
	    string founder_string=ss.str();

	    int endpos;
	    int dist;
	    string line;
	    if (align_umi( padding+umi+padding, padding+founder_string,  max_dist, endpos, dist)) {
		  int ind=(endpos-N_num)/(max_umi_len+N_num)  ;
		  string clustered_founder=founder_view[ind];
		  line = primer_id + "\t" + umi + "\t" + clustered_founder + "\t"+"clustered"+to_string(pool)+"\n";
	    }
	    else{
		  line = primer_id + "\t" + umi + "\t" + umi + "\t"+"founder"+to_string(pool)+"\n";
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

void parallel_founder_find ( vector<UMI_item> low_reads_umi_pool,  ofstream& out, const unsigned num_worker_threads, string curr_primer_id, const unsigned  max_dist, const unsigned max_umi_len  ){
      map<int, vector<UMI_item>> t_umis;
      int subpool_size=ceil(low_reads_umi_pool.size()*1.0/num_worker_threads);
      for (int i = 0; i < num_worker_threads; i++) {
	    int first = i*subpool_size < low_reads_umi_pool.size() ? i*subpool_size : low_reads_umi_pool.size() ;
	    int last = (i+1)*subpool_size<low_reads_umi_pool.size()?(i+1)*subpool_size:low_reads_umi_pool.size() ;
	    vector<UMI_item> one_d_umis;
	    one_d_umis={low_reads_umi_pool.begin()+first, low_reads_umi_pool.begin()+last};
	    t_umis[i]=one_d_umis;
      }

      vector<thread> threads;

      for (int i = 0; i < num_worker_threads; i++) {
	    threads.push_back(thread(founder_find, ref(t_umis[i]) , ref(out), curr_primer_id, max_dist, max_umi_len,  i ));
      }
      for (auto &th : threads) {
	    th.join();
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
	map <string, int> umi_count;
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
	while (cluster_file>>primer_id>>child_umi>>founder_umi>>info) {
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
