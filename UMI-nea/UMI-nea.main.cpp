#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <algorithm>
#include <map>
#include <regex>
#include <stdlib.h>
#include <getopt.h>
#include "UMI-nea.h"

using namespace std;

string in_file;
string out_file;
string updated_count_file;
string estimate_out_file;
int max_dist=2;
int num_thread=10;
int pool_size=1000;
int min_read_founder_user=1;
int min_read_founder=1;
bool user_set_min_read_founder=false;
bool nb_estimate=false;
bool kp_estimate=false;
bool auto_estimate=true;
bool just_estimate=false;
int kp_fail_angle=120;
int kp_angle=0;
int nb_estimated_molecule=0;
int kp_estimated_molecule=0;
int auto_estimated_molecule=0;
int after_rpucut_molecule=0;
int median_rpu=0;
int max_umi_len=-1;
float error_rate=0.001;
bool set_error=false;
bool greedy_mode=false; //quick search mode, once founder found, immediately stop search,  no need to find the low distance founder, good enough for low sequencing error generated UMIs because -m will be small
float nb_lowertail_p=0.001;
int madfolds=3;

void PrintHelp_one_line(string p, string m){
	char para[p.length()+1];
	char message[m.length()+1];
	strcpy(para, p.c_str());
	strcpy(message, m.c_str());
	printf ("%-36s%-10s\n", para, message);
}

void PrintHelp(){
	PrintHelp_one_line("--maxL -l <int|default:-1>", "Max length of UMI, longer UMIs will be trimmed to the length set! Required to set for UMI clustering!");
	PrintHelp_one_line("--in -i <fname>:", "Input file, a tab delim file with three cols: GroupID, UMI, NumofReads, sorted in descending by NumofReads");
	PrintHelp_one_line("--out -o <fname>:", "Output file, output has three cols, GroupID, UMI, FounderUMI");
	PrintHelp_one_line("--minC -m <int|default:2>:", "Min levenshtein distance cutoff, overruled by -e");
	PrintHelp_one_line("--errorR -e <float (0,1]|default:0.001>:", "If set, -m will be calculated based on binomial CI, override -m");
	PrintHelp_one_line("--minF -f <int|default:1>", "Min reads count for a UMI to be founder, overruled by -a or -n or -k");
	PrintHelp_one_line("--nb -n <default:false>:", "Apply negative binomial model to decide min reads for a founder, override -f");
        PrintHelp_one_line("--kp -k <default:false>:", "Apply knee plot to decide min reads for a founder, override -f");
	PrintHelp_one_line("--auto -a <default:true>:", "Combine knee plot strategy and negative binomial model to decide min reads for a founder, override -f");
	PrintHelp_one_line("--just -j <default:false>:", "Just estimate molecule number and rpu cutoff");
	PrintHelp_one_line("--prob -q <float [0, 1]|default:0.001>:", "probability for nb lower tail quantile cutoff in quantification! Not reommended to change");
	PrintHelp_one_line("--greedy -d <default:false>:", "Greedy mode, first founder below cutoff will be selected once found, which speed up computation but affect reprouciblity. Default is false which enforce to find the best founder! Not reommended!");
	PrintHelp_one_line("--thread -t <int >1 |default:10>:", "Num of thread to use, minimal 2");
	PrintHelp_one_line("--pool -p <int >0 |default:1000>:", "Total UMIs to process in each thread at one time");
        PrintHelp_one_line("--help -h:", "Show help");
	exit(1);
}
void PrintOptions(){
	cout << endl<<"********************Input Parameters:************************"<<endl;
	cout << "inFile=" << in_file << endl;
	cout << "outFile=" << out_file << endl;
	cout << "maxlenUMI = " << max_umi_len <<endl;

	if (error_rate != 0)
                cout << "errorRate = "<< error_rate << endl;

	cout << "maxdist = " << max_dist << endl;

	if (nb_estimate){
                cout << "nb-Estimate = ON" << endl;
		cout << "nb_lowertail_p = "<< nb_lowertail_p <<endl;
	}

	if (kp_estimate)
                cout << "kp-Estimate = ON" << endl;

	if (auto_estimate)
                cout << "auto-Estimate = ON" << endl;

	cout << "min_ReadsPerUmi_founer  = " << min_read_founder << endl;

	if (greedy_mode )
                cout << "greedy_Mode = ON"<< endl;

	cout << "threads = " << num_thread << endl;
	cout << "poolSize = " << pool_size << endl;

	cout <<"********************************************"<<endl<<endl;
}
void ProcessArgs(int argc, char** argv)
{
      const char* const short_opts = "l:i:o:m:e:t:p:q:f:agjkngh";
      const option long_opts[] = {
            {"maxL", required_argument, nullptr, 'l'},
            {"in", required_argument, nullptr, 'i'},
            {"out", required_argument, nullptr, 'o'},
            {"minC", required_argument, nullptr, 'm'},
            {"errorR", required_argument, nullptr, 'e'},
	    {"minF", required_argument, nullptr, 'f'},
            {"nb", no_argument, nullptr, 'n'},
            {"kp", no_argument, nullptr, 'k'},
            {"auto", no_argument, nullptr, 'a'},
	    {"just", no_argument, nullptr, 'j'},
	    {"prob", required_argument, nullptr, 'q'},
	    {"greedy", no_argument, nullptr, 'g'},
	    {"thread", required_argument, nullptr, 't'},
	    {"pool", required_argument, nullptr, 'p'},
            {"help", no_argument, nullptr, 'h'},
            {nullptr, no_argument, nullptr, 0}
      };

      while (true)
      {
	    const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);
            if (-1 == opt)
                  break;
            switch (opt)
            {
		  case 'l':
                        max_umi_len = atoi(optarg);
                        break;
                  case 'i':
                        in_file = string(optarg);
                        break;
                  case 'o':
                        out_file = string(optarg);
                        break;
		  case 'm':
                        max_dist = atoi(optarg);
                        break;
                  case 'e':
			error_rate = atof(optarg);
			set_error=true;
			break;
		  case 'f':
			user_set_min_read_founder=true;
			min_read_founder_user = atoi(optarg);
			break;
		  case 'n':
                        nb_estimate = true;
			auto_estimate = false;
                        break;
		  case 'k':
			kp_estimate = true;
			auto_estimate= false;
			break;
		  case 'a':
			auto_estimate = true;
			break;
		  case 'j':
			just_estimate= true;
			max_umi_len=1; //max_umi_len does not matter when -j
			break;
		  case 'q':
			nb_lowertail_p = atof(optarg);
                        break;
		  case 'g':
                        greedy_mode = true;
                        break;
		  case 't':
                        num_thread  = atoi(optarg);
                        break;
                  case 'p':
                        pool_size  = atoi(optarg);
                        break;

		  case 'h': // -h or --help
                  case '?': // Unrecognized option
                  default:
                        PrintHelp();
                        break;
            }
      }
}

void capture_input_error(){
	ifstream my_file(in_file);
	if ( (auto_estimate && nb_estimate) || (auto_estimate && kp_estimate) || (nb_estimate && kp_estimate) ) {
                cerr<<"!!!!!!!!!!!!INPUT ERROR:\n"<<"You can only select one of the following -a, -n, -k"<<endl<<"!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
                PrintHelp();
        }
        if (! my_file){
            cerr<<"!!!!!!!!!!!!INPUT ERROR:\n"<<"input file "<<in_file<<" not readble"<<endl<<"!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
            PrintHelp();
        }
        if ( max_umi_len < 0){
                cerr<<"!!!!!!!!!!!!INPUT ERROR:\n"<<"max_umi_len="<<max_umi_len<<"; was Not set or Not set correctly"<<endl<<"!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
                PrintHelp();
        }
        if ( error_rate <= 0 || error_rate > 1 ){
		cerr<<"!!!!!!!!!!!!INPUT ERROR:\n"<<"error_rate="<<error_rate<<"; error_rate must set within (0,1]"<<endl<<"!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
                PrintHelp();
	}
	if (max_dist < 0) {
		cerr<<"!!!!!!!!!!!!INPUT ERROR:\n"<<"max_dist="<<max_dist<<"; Not 0 or positive number"<<endl<<"!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
                PrintHelp();

	}
	if ( pool_size < 1 ){
		cerr<<"!!!!!!!!!!!!INPUT ERROR:\n"<<"pool_size="<<pool_size <<", Not more than 0"<<endl<<"!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
                PrintHelp();
	}
	if (num_thread <2){
		cerr<<"!!!!!!!!!!!!INPUT ERROR:\n"<<"num_thread="<<num_thread<<", Not more than 1"<<endl<<"!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
                PrintHelp();
	}
	if (nb_lowertail_p <0 || nb_lowertail_p >1 ){
		 cerr<<"!!!!!!!!!!!!INPUT ERROR:\n"<<"nb_lowertail_p="<<nb_lowertail_p<<", must set within [0,1]"<<endl<<"!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
                PrintHelp();
	}
	if (min_read_founder_user<0){
		cerr<<"!!!!!!!!!!!!INPUT ERROR:\n"<<"min_read_founder="<<min_read_founder<<", must bigger than 0"<<endl<<"!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
                PrintHelp();
	}
}

int main(int argc, char **argv){
        ProcessArgs(argc, argv);
	capture_input_error();
	if (set_error){
                max_dist=calculate_dist_upper_bound (error_rate, max_umi_len );
        }
        if (out_file.empty())
		out_file=in_file + ".UMI.clustered" ;
	updated_count_file=out_file + ".umi.reads.count" ;
        estimate_out_file=out_file + ".estimate" ;
	if (nb_estimate){
		fit_nb_model(in_file, nb_lowertail_p, madfolds, min_read_founder,  nb_estimated_molecule, median_rpu);
		after_rpucut_molecule=nb_estimated_molecule;
		if (just_estimate)
			cout<<"NB_estimate\tON\nmedian_rpu\t"<<median_rpu<<endl<<"rpu_cutoff\t"<<min_read_founder<<endl<<"estimated_molecules\t"<<nb_estimated_molecule<<endl<<"after_rpu-cutoff_molecules\t"<<after_rpucut_molecule<<endl;
	}
	else if (kp_estimate){
		fit_knee_plot ( in_file,  min_read_founder,  kp_estimated_molecule, kp_angle, median_rpu, after_rpucut_molecule);
		if (just_estimate)
			cout<<"KP_estimate\tON\nknee_angle\t"<<kp_angle<<endl<<"median_rpu\t"<<median_rpu<<endl<<"rpu_cutoff\t"<<min_read_founder<<endl<<"estimated_molecules\t"<<kp_estimated_molecule<<endl<<"after_rpu-cutoff_molecules\t"<<after_rpucut_molecule<<endl;
	}
	else if (auto_estimate){
		fit_knee_plot ( in_file,  min_read_founder,  kp_estimated_molecule, kp_angle, median_rpu , after_rpucut_molecule);
		if (kp_angle > kp_fail_angle ){
			fit_nb_model(in_file, nb_lowertail_p, madfolds, min_read_founder,  nb_estimated_molecule, median_rpu);
			after_rpucut_molecule=nb_estimated_molecule;
			if (just_estimate)
				cout<<"NB_estimate\tON\nmedian_rpu\t"<<median_rpu<<endl<<"rpu_cutoff\t"<<min_read_founder<<endl<<"estimated_molecules\t"<<nb_estimated_molecule<<endl<<"after_rpu-cutoff_molecules\t"<<after_rpucut_molecule<<endl;
		}
		else{
			if (just_estimate)
				cout<<"KP_estimate\tON\nknee_angle\t"<<kp_angle<<endl<<"median_rpu\t"<<median_rpu<<endl<<"rpu_cutoff\t"<<min_read_founder<<endl<<"estimated_molecules\t"<<kp_estimated_molecule<<endl<<"after_rpu-cutoff_molecules\t"<<after_rpucut_molecule<<endl;
		}
	}
	if (user_set_min_read_founder)
                min_read_founder=min_read_founder_user;
	if (!just_estimate)
		PrintOptions();
	else
		exit(0);

	UMI_clustering_parameters parameters;
	parameters.max_umi_len=max_umi_len;
	parameters.max_dist=max_dist;
	parameters.min_read_founder=min_read_founder;
	parameters.greedy_mode=greedy_mode;
	parameters.thread=num_thread;
	parameters.pool_size=pool_size;
	clustering_umis(in_file, out_file,  parameters );
	update_umi_reads_count(updated_count_file, in_file, out_file);
	ofstream e_out_file(estimate_out_file, ios::out );
	if (nb_estimate)
	//negative binomial to estimate input
		do_nb(updated_count_file, nb_lowertail_p, madfolds, min_read_founder,  nb_estimated_molecule, median_rpu, e_out_file);
	else if (kp_estimate)
	//Knee plot to esitamte input
		do_kp(updated_count_file,  min_read_founder,  kp_estimated_molecule, kp_angle,  median_rpu,  e_out_file);
	else if (auto_estimate){
		fit_knee_plot ( updated_count_file,  min_read_founder,  kp_estimated_molecule, kp_angle, median_rpu, after_rpucut_molecule);
		if (kp_angle > kp_fail_angle)
			do_nb(updated_count_file, nb_lowertail_p, madfolds, min_read_founder,  nb_estimated_molecule, median_rpu, e_out_file);
		else
			do_kp(updated_count_file,  min_read_founder,  kp_estimated_molecule, kp_angle, median_rpu,  e_out_file);
	}
	else{
	// just do clustering and not estimating for cutoff(none of -n -k -a were selected), and estimated molecule will simply equal to # of clustered UMI groups
		int estimate_molecules=count_umi(in_file);
		estimate_molecules=count_umi(updated_count_file);
		cout<<"Just do UMI clustering:"<<"\t"<<"rpu_cutoff=1"<<"\t"<<"estimate_molecules="<<estimate_molecules<<endl;
		e_out_file<<"rpu_cutoff\t"<<1<<endl;
                e_out_file<<"estimated_molecules\t"<<estimate_molecules<<endl;
	}
}
