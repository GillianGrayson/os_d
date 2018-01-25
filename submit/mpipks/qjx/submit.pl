use File::Copy;
use Data::Dumper;
use Cwd;
use Math::Trig;
$dir = getcwd;

$data_path = "/data/biophys/yusipov/os_d/qjx_results";

$PI = 3.1415926535897932384626433832795;

$num_runs = 100;

$eps_exp_shift = 0.1;
$eps_start = 1.0e-8;
$eps_num = 1;

for($curr_cd_eps_id = 0; $curr_cd_eps_id < $eps_num; $curr_cd_eps_id += 1)
{	
	$eps_mult = 10.0 ** ($eps_exp_shift * $curr_cd_eps_id);

	$curr_cd_eps = $eps_start * $eps_mult;

	for ($curr_cd_dim = 1; $curr_cd_dim <= 1; $curr_cd_dim += 1)
	{
		for($curr_U = 0.01; $curr_U <= 0.75000001; $curr_U += 0.01)
		{
			print "curr_U = $curr_U\n";
			print "curr_cd_eps = $curr_cd_eps\n";
			print "curr_cd_dim = $curr_cd_dim\n";
			
			$sys_id = 0;
			$task_id = 0;
			$is_debug = 0;
			$is_pp = 1;
			$init_fn = "";
			$path = "";
			$num_threads = 1;
			$qj_deep = 16;
			$qj_num_tp_periods = 1000;
			$qj_num_obs_periods = 1000;
			$qj_num_trajectories = 2;
			$qj_seed = 0;
			$qj_mns = 1000000;
			
			$type_lpn = 0;
			$eps_lpn = 1.0e-3;
			$delta_up_lpn = 1.0e-2;
			$delta_down_lpn = 1.0e-13;
			$is_obs_dump = 1;
			$is_adr_dump_sep = 0;
			$is_adr_dump_avg = 0;
			$is_evo_dump_sep = 1;
			$is_evo_dump_avg = 0;
			$dump_type = 0;
			$num_dumps = 1000;
			$N = 400;
			$diss_type = 1;
			$diss_gamma = 0.1;
			$diss_phase = 0.0;
			$drv_type = 0;
			$drv_ampl = 1.5;
			$drv_freq = 1.0;
			$drv_phase = 0.0;
			$prm_E = 1.0;
			$prm_U = $curr_U;
			$prm_J = 1.0;
			$start_type = 0;
			$start_state = 0;
			$cd_num_sub_steps = 128;
			$cd_dim = $curr_cd_dim;
			$cd_eps = $curr_cd_eps;
			$cd_dump_deep = 0;
			
			$diss_gamma_str = sprintf("%.4f", $diss_gamma);
			$diss_phase_str = sprintf("%.4f", $diss_phase);

			$drv_ampl_str = sprintf("%.4f", $drv_ampl);
			$drv_freq_str = sprintf("%.4f", $drv_freq);
			$drv_phase_str = sprintf("%.4f", $drv_phase);
			
			$prm_E_str = sprintf("%.4f", $prm_E);
			$prm_U_str = sprintf("%.4f", $prm_U);
			$prm_J_str = sprintf("%.4f", $prm_J);

			$cd_eps_str = sprintf("%.10f", $cd_eps);
			
			for($seed = 0; $seed < 1; $seed += 1)
			{
				$start = 0;
				$finish = $qj_num_trajectories * $num_runs;
				$step = $qj_num_trajectories;
				%exp = ();
				$i = $start;
			
				sub ForderName_task_2{
					$key_str = $_[0];
				
					return  "$data_path/main_${sys_id}_${task_id}/qj_${qj_deep}_${qj_num_tp_periods}_${qj_num_obs_periods}/ci_${cd_num_sub_steps}_${cd_dim}_${cd_eps_str}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}/drv_${drv_type}_${drv_ampl_str}_${drv_freq_str}_${drv_phase_str}/prm_${prm_E_str}_${prm_U_str}_${prm_J_str}/start_${start_type}_${start_state}/ss_${key_str}";
				}
				
				sub ForderName{
					$key_str = $_[0];
					
					return  "$data_path/main_${sys_id}_${task_id}/qj_${qj_deep}_${qj_num_tp_periods}_${qj_num_obs_periods}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}/drv_${drv_type}_${drv_ampl_str}_${drv_freq_str}_${drv_phase_str}/prm_${prm_E_str}_${prm_U_str}_${prm_J_str}/start_${start_type}_${start_state}/ss_${key_str}";
				}
				
				
				if ($task_id == 2)
				{
					mkdir "$data_path";
					mkdir "$data_path/main_${sys_id}_${task_id}";
					mkdir "$data_path/main_${sys_id}_${task_id}/qj_${qj_deep}_${qj_num_tp_periods}_${qj_num_obs_periods}";
					mkdir "$data_path/main_${sys_id}_${task_id}/qj_${qj_deep}_${qj_num_tp_periods}_${qj_num_obs_periods}/ci_${cd_num_sub_steps}_${cd_dim}_${cd_eps_str}";
					mkdir "$data_path/main_${sys_id}_${task_id}/qj_${qj_deep}_${qj_num_tp_periods}_${qj_num_obs_periods}/ci_${cd_num_sub_steps}_${cd_dim}_${cd_eps_str}/N_${N}";
					mkdir "$data_path/main_${sys_id}_${task_id}/qj_${qj_deep}_${qj_num_tp_periods}_${qj_num_obs_periods}/ci_${cd_num_sub_steps}_${cd_dim}_${cd_eps_str}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}";
					mkdir "$data_path/main_${sys_id}_${task_id}/qj_${qj_deep}_${qj_num_tp_periods}_${qj_num_obs_periods}/ci_${cd_num_sub_steps}_${cd_dim}_${cd_eps_str}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}/drv_${drv_type}_${drv_ampl_str}_${drv_freq_str}_${drv_phase_str}";
					mkdir "$data_path/main_${sys_id}_${task_id}/qj_${qj_deep}_${qj_num_tp_periods}_${qj_num_obs_periods}/ci_${cd_num_sub_steps}_${cd_dim}_${cd_eps_str}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}/drv_${drv_type}_${drv_ampl_str}_${drv_freq_str}_${drv_phase_str}/prm_${prm_E_str}_${prm_U_str}_${prm_J_str}";
					mkdir "$data_path/main_${sys_id}_${task_id}/qj_${qj_deep}_${qj_num_tp_periods}_${qj_num_obs_periods}/ci_${cd_num_sub_steps}_${cd_dim}_${cd_eps_str}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}/drv_${drv_type}_${drv_ampl_str}_${drv_freq_str}_${drv_phase_str}/prm_${prm_E_str}_${prm_U_str}_${prm_J_str}/start_${start_type}_${start_state}";
				}
				else
				{
					mkdir "$data_path";
					mkdir "$data_path/main_${sys_id}_${task_id}";
					mkdir "$data_path/main_${sys_id}_${task_id}/qj_${qj_deep}_${qj_num_tp_periods}_${qj_num_obs_periods}";
					mkdir "$data_path/main_${sys_id}_${task_id}/qj_${qj_deep}_${qj_num_tp_periods}_${qj_num_obs_periods}/N_${N}";
					mkdir "$data_path/main_${sys_id}_${task_id}/qj_${qj_deep}_${qj_num_tp_periods}_${qj_num_obs_periods}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}";
					mkdir "$data_path/main_${sys_id}_${task_id}/qj_${qj_deep}_${qj_num_tp_periods}_${qj_num_obs_periods}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}/drv_${drv_type}_${drv_ampl_str}_${drv_freq_str}_${drv_phase_str}";
					mkdir "$data_path/main_${sys_id}_${task_id}/qj_${qj_deep}_${qj_num_tp_periods}_${qj_num_obs_periods}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}/drv_${drv_type}_${drv_ampl_str}_${drv_freq_str}_${drv_phase_str}/prm_${prm_E_str}_${prm_U_str}_${prm_J_str}";
					mkdir "$data_path/main_${sys_id}_${task_id}/qj_${qj_deep}_${qj_num_tp_periods}_${qj_num_obs_periods}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}/drv_${drv_type}_${drv_ampl_str}_${drv_freq_str}_${drv_phase_str}/prm_${prm_E_str}_${prm_U_str}_${prm_J_str}/start_${start_type}_${start_state}";
				}
				
				for($val = $start; $val < $finish; $val += $step)
				{
					if ($task_id == 2)
					{
						$exp{ForderName_task_2($i)} = $val;
					}
					else
					{
						$exp{ForderName($i)} = $val;
					}
					$i++;
				}

				for($i = $start; $i < $finish; $i += $step)
				{
					if ($task_id == 2)
					{
						$key = ForderName_task_2($i); 
					}
					else
					{
						$key = ForderName($i); 
					}
				
					mkdir "$key";
					
					print "$i \n";
					
					open( WF,">$key/config.txt"); 
					print WF "sys_id $sys_id \n"; 
					print WF "task_id $task_id \n";
					print WF "is_debug $is_debug \n";
					print WF "is_pp $is_pp \n";
					print WF "init_fn $init_fn \n";
					print WF "path $path \n";
					print WF "num_threads $num_threads \n";
					print WF "qj_deep $qj_deep \n";
					print WF "qj_num_tp_periods $qj_num_tp_periods \n";
					print WF "qj_num_obs_periods $qj_num_obs_periods \n";
					print WF "qj_num_trajectories $qj_num_trajectories \n";
					print WF "qj_seed $i \n";
					print WF "qj_mns $qj_mns \n";
					close WF;
					
					open( WF,">$key/params.txt"); 
					print WF "type_lpn $type_lpn \n"; 
					print WF "eps_lpn $eps_lpn \n";
					print WF "delta_up_lpn $delta_up_lpn \n";
					print WF "delta_down_lpn $delta_down_lpn \n";
					print WF "is_obs_dump $is_obs_dump \n";
					print WF "is_adr_dump_sep $is_adr_dump_sep \n";
					print WF "is_adr_dump_avg $is_adr_dump_avg \n";
					print WF "is_evo_dump_sep $is_evo_dump_sep \n";
					print WF "is_evo_dump_avg $is_evo_dump_avg \n";
					print WF "dump_type $dump_type \n";
					print WF "num_dumps $num_dumps \n";
					print WF "N $N \n";
					print WF "diss_type $diss_type \n";
					print WF "diss_gamma $diss_gamma \n";
					print WF "diss_phase $diss_phase \n";
					print WF "drv_type $drv_type \n";
					print WF "drv_ampl $drv_ampl \n";
					print WF "drv_freq $drv_freq \n";
					print WF "drv_phase $drv_phase \n";
					print WF "prm_E $prm_E \n";
					print WF "prm_U $prm_U \n";
					print WF "prm_J $prm_J \n";
					print WF "start_type $start_type \n";
					print WF "start_state $start_state \n";
					print WF "cd_num_sub_steps $cd_num_sub_steps \n";
					print WF "cd_dim $cd_dim \n";
					print WF "cd_eps $cd_eps \n";
					print WF "cd_dump_deep $cd_dump_deep \n";
					close WF;
					
					if ($task_id == 2)
					{
						$test_file = sprintf('%s/ci_qjrnd(%d_%d)_N(%d)_diss(%d_%0.4f_%0.4f)_drv(%d_%0.4f_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d).txt', $key, $i, $qj_mns, $N, $diss_type, $diss_gamma, $diss_phase, $drv_type, $drv_ampl, $drv_freq, $drv_phase, $prm_E, $prm_U, $prm_J, $start_type, $start_state);

					}
					elsif($task_id == 0)
					{
						$test_file = sprintf('%s/mean_qjrnd(%d_%d)_N(%d)_diss(%d_%0.4f_%0.4f)_drv(%d_%0.4f_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d).txt', $key, $i, $qj_mns, $N, $diss_type, $diss_gamma, $diss_phase, $drv_type, $drv_ampl, $drv_freq, $drv_phase, $prm_E, $prm_U, $prm_J, $start_type, $start_state);

					}
					else
					{
						$test_file = sprintf('%s/adr_avg_qjrnd(%d_%d)_N(%d)_diss(%d_%0.4f_%0.4f)_drv(%d_%0.4f_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d).txt', $key, $i, $qj_mns, $N, $diss_type, $diss_gamma, $diss_phase, $drv_type, $drv_ampl, $drv_freq, $drv_phase, $prm_E, $prm_U, $prm_J, $start_type, $start_state);

					}
					
					#print "$test_file \n";
					
					unless (-e "$test_file")
					{	
						print "qsub -wd $dir run_1_thread.sh $key \n";
						system "qsub -wd $dir run_1_thread.sh $key";
					}
				}
			}
		}
	}
}
