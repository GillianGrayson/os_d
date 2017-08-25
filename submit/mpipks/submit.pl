use File::Copy;
use Data::Dumper;
use Cwd;
use Math::Trig;
$dir = getcwd;

$data_path = "/data/biophys/yusipov/os_d";
$input_path = "qj_input";

$prefix = "qj_results/delta_0.1000/tt_1000";

for($curr_U = 0.01; $curr_U <= 0.750001; $curr_U += 0.01)
{
	print "curr_U = $curr_U\n";
	
	$N = 500;
	$U = $curr_U;
	
	$num_periods = 1000;
	$init_state_id = int($N/2);
	$propagation_type = 3;
	$num_dumps = 1001;
	$dump_type = 0;
	$num_periods_in_trans_proc = 1000;
	$num_omp_threads = 1;
	$num_trajectories = 2;
	$rnd_max = 2000000;
	$rnd_cur = 0;
	$calc_characteristics = 0;
	$dump_rho = 0;
	$mean_low_limit = 0.00;
	$mean_high_limit = 1.00;
	$avg_dump = 0;
	$dump_characteristics = 0;
	$btw_jump_times = 0;
	$borders_type = 0;
	$stationary = 0;
	$num_e_intervals = 1000;
	$num_mc_intervals = 1000;
	$energy_min = -2.5;
	$energy_max = 2.5;
	$after_dump = 0;
	$double_scale_dump = 0;
	$deep_characteristic = 0;
	$mc_specific = 0;
	$mc_type = 1;
	$num_att_trajectories = 0;
	$var_eps = 0.001;
	$delta_lim = 0.1;

	for($seed = 0; $seed < 1; $seed+=1)
	{
		$start = 0;
		$finish = 100;
		%exp = ();
		$i = 0;
		$i = $start;
	
		$U_str = sprintf("%.4f", $U);
		
		$input_file_name = sprintf('%s/main_data_N%d_U%0.4f.bin', $input_path, $N, $U);
		$aux_file_name 	 = sprintf('%s/main_data_N%d_U%0.4f.bin', $input_path, $N, $U);

		sub ForderName{
			$key_str = $_[0];
			
			return  "$data_path/$prefix/N_${N}/U_${U_str}/rnd_${key_str}";
		}

		mkdir "$data_path/$prefix";
		mkdir "$data_path/$prefix/N_${N}";
		mkdir "$data_path/$prefix/N_${N}/U_${U_str}";
		
		
		for($val = $start; $val < $finish; $val+=1)
		{
			$exp{ForderName($i)} = $val;
			$i++;
		}

		for($i = $start; $i < $finish; $i++)
		{
			$key = ForderName($i);    
			mkdir "$key";
			
			$rnd = $rnd_cur + $exp{$key} * $num_trajectories;
			
			print "$rnd \n";
			
			open( WF,">$key/config.txt"); 
			print WF "num_periods = $num_periods \n"; 
			print WF "init_state_id = $init_state_id \n";
			print WF "propagation_type = $propagation_type \n";
			print WF "num_dumps = $num_dumps \n";
			print WF "dump_type = $dump_type \n";
			print WF "num_periods_in_trans_proc = $num_periods_in_trans_proc \n";
			print WF "num_omp_threads = $num_omp_threads \n";
			print WF "num_trajectories = $num_trajectories \n"; 
			print WF "rnd_max = $rnd_max \n";
			print WF "rnd_cur = $rnd \n";
			print WF "calc_characteristics = $calc_characteristics \n";
			print WF "dump_rho = $dump_rho \n";
			print WF "mean_low_limit = $mean_low_limit \n";
			print WF "mean_high_limit = $mean_high_limit \n";
			print WF "avg_dump = $avg_dump \n";
			print WF "dump_characteristics = $dump_characteristics \n";
			print WF "btw_jump_times = $btw_jump_times \n";
			print WF "borders_type = $borders_type \n";
			print WF "stationary = $stationary \n";
			print WF "num_e_intervals = $num_e_intervals \n";
			print WF "num_mc_intervals = $num_mc_intervals \n";
			print WF "energy_min = $energy_min \n";
			print WF "energy_max = $energy_max \n";
			print WF "after_dump = $after_dump \n";
			print WF "double_scale_dump = $double_scale_dump \n";
			print WF "deep_characteristic = $deep_characteristic \n";
			print WF "mc_specific = $mc_specific \n";
			print WF "mc_type = $mc_type \n";
			print WF "num_att_trajectories = $num_att_trajectories \n";
			print WF "var_eps = $var_eps \n";
			print WF "delta_lim = $delta_lim \n";
			close WF;

			$test_file = "periods_evo.txt";
			
			unless (-e "$key/$test_file")
			{	
				print "qsub -wd $dir run_1_thread.sh $key $input_file_name $aux_file_name \n";
				system "qsub -wd $dir run_1_thread.sh $key $input_file_name $aux_file_name ";
			}
		}
	}

}
