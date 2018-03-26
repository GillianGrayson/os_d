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
		for($curr_U = 0.01; $curr_U <= 0.010001; $curr_U += 0.01)
		{
			print "curr_U = $curr_U\n";
			print "curr_cd_eps = $curr_cd_eps\n";
			print "curr_cd_dim = $curr_cd_dim\n";
			
			$sys_id = 0;
			$task_id = 0;
			$prop_id = 0;
			$is_debug = 0;
			$is_pp = 1;
			$init_fn = "";
			$path = "";
			$num_threads = 1;
			$seed = 0;
			$mns = 1000000;
			$num_trajectories = 2;
			$num_tp_periods = 1000;
			$num_obs_periods = 100;
			$ex_deep = 16;
			$rk_ns = 10000;
			
			$lpn_type = 0;
			$lpn_eps = 1.0e-1;
			$lpn_delta_up = 1.0e-2;
			$lpn_delta_down = 1.0e-13;
			$dump_obs = 1;
			$dump_adr_sep = 0;
			$dump_adr_avg = 0;
			$dump_evo_sep = 1;
			$dump_evo_avg = 0;
			$dump_type = 0;
			$dump_num = 100;
			$N = 100;
			$diss_type = 0;
			$diss_gamma = 0.1;
			$diss_phase = 0.0;
			$drv_dimer_type = 0;
			$drv_dimer_ampl = 1.5;
			$drv_dimer_freq = 1.0;
			$drv_dimer_phase = 0.0;
			$prm_dimer_E = 1.0;
			$prm_dimer_U = $curr_U;
			$prm_dimer_J = 1.0;
			$start_type = 0;
			$start_state = 0;
			$cd_dim = $curr_cd_dim;
			$cd_eps = $curr_cd_eps;
			$deep_num_steps = 128;
			$deep_dump = 0;
			$jump = 0;
			
			$diss_gamma_str = sprintf("%.4f", $diss_gamma);
			$diss_phase_str = sprintf("%.4f", $diss_phase);

			$drv_dimer_ampl_str = sprintf("%.4f", $drv_dimer_ampl);
			$drv_dimer_freq_str = sprintf("%.4f", $drv_dimer_freq);
			$drv_dimer_phase_str = sprintf("%.4f", $drv_dimer_phase);
			
			$prm_dimer_E_str = sprintf("%.4f", $prm_dimer_E);
			$prm_dimer_U_str = sprintf("%.4f", $prm_dimer_U);
			$prm_dimer_J_str = sprintf("%.4f", $prm_dimer_J);
			
			$lpn_eps_str = sprintf("%.4f", $lpn_eps);
			$lpn_delta_up_str = sprintf("%.4f", $lpn_delta_up);

			$cd_eps_str = sprintf("%.10f", $cd_eps);
			
			for($seed = 0; $seed < 1; $seed += 1)
			{
				$start = 0;
				$finish = $num_trajectories * $num_runs;
				$step = $num_trajectories;
				%exp = ();
				$i = $start;
				
				sub ForderName_task_0{
					$key_str = $_[0];
				
					return  "$data_path/main_${sys_id}_${prop_id}_${task_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/lpn_${lpn_eps_str}_${lpn_delta_up_str}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}/drv_${drv_dimer_type}_${drv_dimer_ampl_str}_${drv_dimer_freq_str}_${drv_dimer_phase_str}/prm_${prm_dimer_E_str}_${prm_dimer_U_str}_${prm_dimer_J_str}/start_${start_type}_${start_state}/ss_${key_str}";
				}
			
				sub ForderName_task_2{
					$key_str = $_[0];
				
					return  "$data_path/main_${sys_id}_${prop_id}_${task_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/ci_${deep_num_steps}_${cd_dim}_${cd_eps_str}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}/drv_${drv_dimer_type}_${drv_dimer_ampl_str}_${drv_dimer_freq_str}_${drv_dimer_phase_str}/prm_${prm_dimer_E_str}_${prm_dimer_U_str}_${prm_dimer_J_str}/start_${start_type}_${start_state}/ss_${key_str}";
				}
				
				sub ForderName{
					$key_str = $_[0];
					
					return  "$data_path/main_${sys_id}_${prop_id}_${task_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}/drv_${drv_dimer_type}_${drv_dimer_ampl_str}_${drv_dimer_freq_str}_${drv_dimer_phase_str}/prm_${prm_dimer_E_str}_${prm_dimer_U_str}_${prm_dimer_J_str}/start_${start_type}_${start_state}/ss_${key_str}";
				}
				
				if ($task_id == 0)
				{
					mkdir "$data_path";
					mkdir "$data_path/main_${sys_id}_${prop_id}_${task_id}";
					mkdir "$data_path/main_${sys_id}_${prop_id}_${task_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}";
					mkdir "$data_path/main_${sys_id}_${prop_id}_${task_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/lpn_${lpn_eps_str}_${lpn_delta_up_str}";
					mkdir "$data_path/main_${sys_id}_${prop_id}_${task_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/lpn_${lpn_eps_str}_${lpn_delta_up_str}/N_${N}";
					mkdir "$data_path/main_${sys_id}_${prop_id}_${task_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/lpn_${lpn_eps_str}_${lpn_delta_up_str}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}";
					mkdir "$data_path/main_${sys_id}_${prop_id}_${task_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/lpn_${lpn_eps_str}_${lpn_delta_up_str}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}/drv_${drv_dimer_type}_${drv_dimer_ampl_str}_${drv_dimer_freq_str}_${drv_dimer_phase_str}";
					mkdir "$data_path/main_${sys_id}_${prop_id}_${task_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/lpn_${lpn_eps_str}_${lpn_delta_up_str}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}/drv_${drv_dimer_type}_${drv_dimer_ampl_str}_${drv_dimer_freq_str}_${drv_dimer_phase_str}/prm_${prm_dimer_E_str}_${prm_dimer_U_str}_${prm_dimer_J_str}";
					mkdir "$data_path/main_${sys_id}_${prop_id}_${task_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/lpn_${lpn_eps_str}_${lpn_delta_up_str}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}/drv_${drv_dimer_type}_${drv_dimer_ampl_str}_${drv_dimer_freq_str}_${drv_dimer_phase_str}/prm_${prm_dimer_E_str}_${prm_dimer_U_str}_${prm_dimer_J_str}/start_${start_type}_${start_state}";
				}
				elsif ($task_id == 2)
				{
					mkdir "$data_path";
					mkdir "$data_path/main_${sys_id}_${prop_id}_${task_id}";
					mkdir "$data_path/main_${sys_id}_${prop_id}_${task_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}";
					mkdir "$data_path/main_${sys_id}_${prop_id}_${task_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/ci_${deep_num_steps}_${cd_dim}_${cd_eps_str}";
					mkdir "$data_path/main_${sys_id}_${prop_id}_${task_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/ci_${deep_num_steps}_${cd_dim}_${cd_eps_str}/N_${N}";
					mkdir "$data_path/main_${sys_id}_${prop_id}_${task_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/ci_${deep_num_steps}_${cd_dim}_${cd_eps_str}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}";
					mkdir "$data_path/main_${sys_id}_${prop_id}_${task_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/ci_${deep_num_steps}_${cd_dim}_${cd_eps_str}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}/drv_${drv_dimer_type}_${drv_dimer_ampl_str}_${drv_dimer_freq_str}_${drv_dimer_phase_str}";
					mkdir "$data_path/main_${sys_id}_${prop_id}_${task_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/ci_${deep_num_steps}_${cd_dim}_${cd_eps_str}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}/drv_${drv_dimer_type}_${drv_dimer_ampl_str}_${drv_dimer_freq_str}_${drv_dimer_phase_str}/prm_${prm_dimer_E_str}_${prm_dimer_U_str}_${prm_dimer_J_str}";
					mkdir "$data_path/main_${sys_id}_${prop_id}_${task_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/ci_${deep_num_steps}_${cd_dim}_${cd_eps_str}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}/drv_${drv_dimer_type}_${drv_dimer_ampl_str}_${drv_dimer_freq_str}_${drv_dimer_phase_str}/prm_${prm_dimer_E_str}_${prm_dimer_U_str}_${prm_dimer_J_str}/start_${start_type}_${start_state}";
				}
				else
				{
					mkdir "$data_path";
					mkdir "$data_path/main_${sys_id}_${prop_id}_${task_id}";
					mkdir "$data_path/main_${sys_id}_${prop_id}_${task_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}";
					mkdir "$data_path/main_${sys_id}_${prop_id}_${task_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/N_${N}";
					mkdir "$data_path/main_${sys_id}_${prop_id}_${task_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}";
					mkdir "$data_path/main_${sys_id}_${prop_id}_${task_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}/drv_${drv_dimer_type}_${drv_dimer_ampl_str}_${drv_dimer_freq_str}_${drv_dimer_phase_str}";
					mkdir "$data_path/main_${sys_id}_${prop_id}_${task_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}/drv_${drv_dimer_type}_${drv_dimer_ampl_str}_${drv_dimer_freq_str}_${drv_dimer_phase_str}/prm_${prm_dimer_E_str}_${prm_dimer_U_str}_${prm_dimer_J_str}";
					mkdir "$data_path/main_${sys_id}_${prop_id}_${task_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}/drv_${drv_dimer_type}_${drv_dimer_ampl_str}_${drv_dimer_freq_str}_${drv_dimer_phase_str}/prm_${prm_dimer_E_str}_${prm_dimer_U_str}_${prm_dimer_J_str}/start_${start_type}_${start_state}";
				}
				
				for($val = $start; $val < $finish; $val += $step)
				{
					if ($task_id == 0)
					{
						$exp{ForderName_task_0($i)} = $val;
					}
					elsif ($task_id == 2)
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
					if ($task_id == 0)
					{
						$key = ForderName_task_0($i);
					}
					elsif ($task_id == 2)
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
					print WF "prop_id $prop_id \n";
					print WF "is_debug $is_debug \n";
					print WF "is_pp $is_pp \n";
					print WF "init_fn $init_fn \n";
					print WF "path $path \n";
					print WF "num_threads $num_threads \n";
					print WF "num_tp_periods $num_tp_periods \n";
					print WF "num_obs_periods $num_obs_periods \n";
					print WF "num_trajectories $num_trajectories \n";
					print WF "seed $i \n";
					print WF "mns $mns \n";
					print WF "ex_deep $ex_deep \n";
					print WF "rk_ns $rk_ns \n";
					close WF;
					
					open( WF,">$key/params.txt"); 
					print WF "lpn_type $lpn_type \n"; 
					print WF "lpn_eps $lpn_eps \n";
					print WF "lpn_delta_up $lpn_delta_up \n";
					print WF "lpn_delta_down $lpn_delta_down \n";
					print WF "dump_obs $dump_obs \n";
					print WF "dump_adr_sep $dump_adr_sep \n";
					print WF "dump_adr_avg $dump_adr_avg \n";
					print WF "dump_evo_sep $dump_evo_sep \n";
					print WF "dump_evo_avg $dump_evo_avg \n";
					print WF "dump_type $dump_type \n";
					print WF "dump_num $dump_num \n";
					print WF "N $N \n";
					print WF "diss_type $diss_type \n";
					print WF "diss_gamma $diss_gamma \n";
					print WF "diss_phase $diss_phase \n";
					print WF "drv_dimer_type $drv_dimer_type \n";
					print WF "drv_dimer_ampl $drv_dimer_ampl \n";
					print WF "drv_dimer_freq $drv_dimer_freq \n";
					print WF "drv_dimer_phase $drv_dimer_phase \n";
					print WF "prm_dimer_E $prm_dimer_E \n";
					print WF "prm_dimer_U $prm_dimer_U \n";
					print WF "prm_dimer_J $prm_dimer_J \n";
					print WF "start_type $start_type \n";
					print WF "start_state $start_state \n";
					print WF "cd_dim $cd_dim \n";
					print WF "cd_eps $cd_eps \n";
					print WF "deep_num_steps $deep_num_steps \n";
					print WF "deep_dump $deep_dump \n";
					print WF "jump $jump \n";
					close WF;
					
					if ($task_id == 2)
					{
						$test_file = sprintf('%s/ci_config(%d_%d_%d)_rnd(%d_%d)_N(%d)_diss(%d_%0.4f_%0.4f)_drv(%d_%0.4f_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d).txt', $key, $sys_id, $task_id, $prop_id, $i, $mns, $N, $diss_type, $diss_gamma, $diss_phase, $drv_dimer_type, $drv_dimer_ampl, $drv_dimer_freq, $drv_dimer_phase, $prm_dimer_E, $prm_dimer_U, $prm_dimer_J, $start_type, $start_state);

					}
					elsif($task_id == 0)
					{
						$test_file = sprintf('%s/mean_config(%d_%d_%d)_rnd(%d_%d)_N(%d)_diss(%d_%0.4f_%0.4f)_drv(%d_%0.4f_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d).txt', $key, $sys_id, $task_id, $prop_id, $i, $mns, $N, $diss_type, $diss_gamma, $diss_phase, $drv_dimer_type, $drv_dimer_ampl, $drv_dimer_freq, $drv_dimer_phase, $prm_dimer_E, $prm_dimer_U, $prm_dimer_J, $start_type, $start_state);

					}
					else
					{
						$test_file = sprintf('%s/adr_avg_config(%d_%d_%d)_rnd(%d_%d)_N(%d)_diss(%d_%0.4f_%0.4f)_drv(%d_%0.4f_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d).txt', $key, $sys_id, $task_id, $prop_id, $i, $mns, $N, $diss_type, $diss_gamma, $diss_phase, $drv_dimer_type, $drv_dimer_ampl, $drv_dimer_freq, $drv_dimer_phase, $prm_dimer_E, $prm_dimer_U, $prm_dimer_J, $start_type, $start_state);
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
