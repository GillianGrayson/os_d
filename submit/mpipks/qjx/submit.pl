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
		for($curr_U = 0.01; $curr_U <= 0.750001; $curr_U += 0.01)
		{
			for($curr_jcs_drv_ampl = 0.01; $curr_jcs_drv_ampl <= 0.0100001; $curr_jcs_drv_ampl += 0.01)
			{
				print "curr_U = $curr_U\n";
				print "curr_cd_eps = $curr_cd_eps\n";
				print "curr_cd_dim = $curr_cd_dim\n";
				
				print "curr_jcs_drv_ampl = $curr_jcs_drv_ampl\n";
				
				$sys_id = 0;
				$task_id = 0;
				$prop_id = 0;
				$is_debug = 0;
				$is_pp = 1;
				$init_fn = "";
				$path = "";
				$seed = 0;
				$mns = 1000000;
				$num_threads = 1;
				$num_trajectories = 2;
				$num_tp_periods = 100;
				$num_obs_periods = 100;
				$ex_deep = 16;
				$rk_ns = 10000;
				
				$lpn_type = 0;
				$lpn_eps = 1.0e-2;
				$lpn_eps_change = 1.25;
				$lpn_delta_s_high = 1.0e-3;
				$lpn_delta_s_low = 1.0e-4;
				$lpn_delta_f_high = 1.0e-2;
				$lpn_delta_f_low = 1.0e-10;
				$dump_obs = 1;
				$dump_adr_sep = 0;
				$dump_adr_avg = 0;
				$dump_evo_sep = 1;
				$dump_evo_avg = 0;
				$dump_type = 0;
				$dump_num = 100;
				$N = 500;
				$diss_type = 0;
				$diss_gamma = 0.1;
				$diss_phase = 0.0;
				$dimer_drv_type = 0;
				$dimer_drv_ampl = 1.5;
				$dimer_drv_freq = 1.0;
				$dimer_drv_phase = 0.0;
				$dimer_prm_E = 1.0;
				$dimer_prm_U = $curr_U;
				$dimer_prm_J = 1.0;
				$jcs_drv_part_1 = 0.98;
				$jcs_drv_part_2 = 1.0;
				$jcs_drv_ampl = $curr_jcs_drv_ampl;
				$jcs_prm_alpha = 5.0;
				$start_type = 0;
				$start_state = 0;
				$cd_dim = $curr_cd_dim;
				$cd_eps = $curr_cd_eps;
				$deep_num_steps = 128;
				$deep_dump = 0;
				$jump = 0;
				
				$diss_gamma_str = sprintf("%.4f", $diss_gamma);
				$diss_phase_str = sprintf("%.4f", $diss_phase);

				$dimer_drv_ampl_str = sprintf("%.4f", $dimer_drv_ampl);
				$dimer_drv_freq_str = sprintf("%.4f", $dimer_drv_freq);
				$dimer_drv_phase_str = sprintf("%.4f", $dimer_drv_phase);
				$dimer_prm_E_str = sprintf("%.4f", $dimer_prm_E);
				$dimer_prm_U_str = sprintf("%.4f", $dimer_prm_U);
				$dimer_prm_J_str = sprintf("%.4f", $dimer_prm_J);
				
				$jcs_drv_part_1_str = sprintf("%.4f", $jcs_drv_part_1);
				$jcs_drv_part_2_str = sprintf("%.4f", $jcs_drv_part_2);
				$jcs_drv_ampl_str = sprintf("%.4f", $jcs_drv_ampl);
				$jcs_prm_alpha_str = sprintf("%.4f", $jcs_prm_alpha);
				
				$lpn_eps_str = sprintf("%.4f", $lpn_eps);
				$lpn_delta_f_high_str = sprintf("%.4f", $lpn_delta_f_high);

				$cd_eps_str = sprintf("%.10f", $cd_eps);
				
				for($seed = 0; $seed < 1; $seed += 1)
				{
					$start = 0;
					$finish = $num_trajectories * $num_runs;
					$step = $num_trajectories;
					%exp = ();
					$i = $start;
					
					sub ForderName_sys_0_sys_0_task_0{
						$key_str = $_[0];
					
						return  "$data_path/main_${sys_id}_${prop_id}_${task_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/lpn_${lpn_eps_str}_${lpn_delta_f_high_str}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}/drv_${dimer_drv_type}_${dimer_drv_ampl_str}_${dimer_drv_freq_str}_${dimer_drv_phase_str}/prm_${dimer_prm_E_str}_${dimer_prm_U_str}_${dimer_prm_J_str}/start_${start_type}_${start_state}/ss_${key_str}";
					}
				
					sub ForderName_sys_0_sys_0_task_2{
						$key_str = $_[0];
					
						return  "$data_path/main_${sys_id}_${prop_id}_${task_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/ci_${deep_num_steps}_${cd_dim}_${cd_eps_str}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}/drv_${dimer_drv_type}_${dimer_drv_ampl_str}_${dimer_drv_freq_str}_${dimer_drv_phase_str}/prm_${dimer_prm_E_str}_${dimer_prm_U_str}_${dimer_prm_J_str}/start_${start_type}_${start_state}/ss_${key_str}";
					}
					
					sub ForderName_sys_0{
						$key_str = $_[0];
						
						return  "$data_path/main_${sys_id}_${prop_id}_${task_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}/drv_${dimer_drv_type}_${dimer_drv_ampl_str}_${dimer_drv_freq_str}_${dimer_drv_phase_str}/prm_${dimer_prm_E_str}_${dimer_prm_U_str}_${dimer_prm_J_str}/start_${start_type}_${start_state}/ss_${key_str}";
					}
					
					sub ForderName_sys_1{
						$key_str = $_[0];
						
						return  "$data_path/main_${sys_id}_${task_id}_${prop_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}/drv_${jcs_drv_part_1_str}_${jcs_drv_part_2_str}_${jcs_drv_ampl_str}/prm_${jcs_prm_alpha_str}/start_${start_type}_${start_state}/ss_${key_str}";
					}
					
					
					if($sys_id == 0)
					{
						if ($task_id == 0)
						{
							mkdir "$data_path";
							mkdir "$data_path/main_${sys_id}_${prop_id}_${task_id}";
							mkdir "$data_path/main_${sys_id}_${prop_id}_${task_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}";
							mkdir "$data_path/main_${sys_id}_${prop_id}_${task_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/lpn_${lpn_eps_str}_${lpn_delta_f_high_str}";
							mkdir "$data_path/main_${sys_id}_${prop_id}_${task_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/lpn_${lpn_eps_str}_${lpn_delta_f_high_str}/N_${N}";
							mkdir "$data_path/main_${sys_id}_${prop_id}_${task_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/lpn_${lpn_eps_str}_${lpn_delta_f_high_str}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}";
							mkdir "$data_path/main_${sys_id}_${prop_id}_${task_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/lpn_${lpn_eps_str}_${lpn_delta_f_high_str}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}/drv_${dimer_drv_type}_${dimer_drv_ampl_str}_${dimer_drv_freq_str}_${dimer_drv_phase_str}";
							mkdir "$data_path/main_${sys_id}_${prop_id}_${task_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/lpn_${lpn_eps_str}_${lpn_delta_f_high_str}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}/drv_${dimer_drv_type}_${dimer_drv_ampl_str}_${dimer_drv_freq_str}_${dimer_drv_phase_str}/prm_${dimer_prm_E_str}_${dimer_prm_U_str}_${dimer_prm_J_str}";
							mkdir "$data_path/main_${sys_id}_${prop_id}_${task_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/lpn_${lpn_eps_str}_${lpn_delta_f_high_str}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}/drv_${dimer_drv_type}_${dimer_drv_ampl_str}_${dimer_drv_freq_str}_${dimer_drv_phase_str}/prm_${dimer_prm_E_str}_${dimer_prm_U_str}_${dimer_prm_J_str}/start_${start_type}_${start_state}";
						}
						elsif ($task_id == 2)
						{
							mkdir "$data_path";
							mkdir "$data_path/main_${sys_id}_${prop_id}_${task_id}";
							mkdir "$data_path/main_${sys_id}_${prop_id}_${task_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}";
							mkdir "$data_path/main_${sys_id}_${prop_id}_${task_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/ci_${deep_num_steps}_${cd_dim}_${cd_eps_str}";
							mkdir "$data_path/main_${sys_id}_${prop_id}_${task_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/ci_${deep_num_steps}_${cd_dim}_${cd_eps_str}/N_${N}";
							mkdir "$data_path/main_${sys_id}_${prop_id}_${task_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/ci_${deep_num_steps}_${cd_dim}_${cd_eps_str}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}";
							mkdir "$data_path/main_${sys_id}_${prop_id}_${task_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/ci_${deep_num_steps}_${cd_dim}_${cd_eps_str}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}/drv_${dimer_drv_type}_${dimer_drv_ampl_str}_${dimer_drv_freq_str}_${dimer_drv_phase_str}";
							mkdir "$data_path/main_${sys_id}_${prop_id}_${task_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/ci_${deep_num_steps}_${cd_dim}_${cd_eps_str}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}/drv_${dimer_drv_type}_${dimer_drv_ampl_str}_${dimer_drv_freq_str}_${dimer_drv_phase_str}/prm_${dimer_prm_E_str}_${dimer_prm_U_str}_${dimer_prm_J_str}";
							mkdir "$data_path/main_${sys_id}_${prop_id}_${task_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/ci_${deep_num_steps}_${cd_dim}_${cd_eps_str}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}/drv_${dimer_drv_type}_${dimer_drv_ampl_str}_${dimer_drv_freq_str}_${dimer_drv_phase_str}/prm_${dimer_prm_E_str}_${dimer_prm_U_str}_${dimer_prm_J_str}/start_${start_type}_${start_state}";
						}
						else
						{
							mkdir "$data_path";
							mkdir "$data_path/main_${sys_id}_${prop_id}_${task_id}";
							mkdir "$data_path/main_${sys_id}_${prop_id}_${task_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}";
							mkdir "$data_path/main_${sys_id}_${prop_id}_${task_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/N_${N}";
							mkdir "$data_path/main_${sys_id}_${prop_id}_${task_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}";
							mkdir "$data_path/main_${sys_id}_${prop_id}_${task_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}/drv_${dimer_drv_type}_${dimer_drv_ampl_str}_${dimer_drv_freq_str}_${dimer_drv_phase_str}";
							mkdir "$data_path/main_${sys_id}_${prop_id}_${task_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}/drv_${dimer_drv_type}_${dimer_drv_ampl_str}_${dimer_drv_freq_str}_${dimer_drv_phase_str}/prm_${dimer_prm_E_str}_${dimer_prm_U_str}_${dimer_prm_J_str}";
							mkdir "$data_path/main_${sys_id}_${prop_id}_${task_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}/drv_${dimer_drv_type}_${dimer_drv_ampl_str}_${dimer_drv_freq_str}_${dimer_drv_phase_str}/prm_${dimer_prm_E_str}_${dimer_prm_U_str}_${dimer_prm_J_str}/start_${start_type}_${start_state}";
						}
					}
					elsif ($sys_id == 1)
					{
						if ($task_id == 1)
						{
							mkdir "$data_path";
							mkdir "$data_path/main_${sys_id}_${task_id}_${prop_id}";
							mkdir "$data_path/main_${sys_id}_${task_id}_${prop_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}";
							mkdir "$data_path/main_${sys_id}_${task_id}_${prop_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/N_${N}";
							mkdir "$data_path/main_${sys_id}_${task_id}_${prop_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}";
							mkdir "$data_path/main_${sys_id}_${task_id}_${prop_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}/drv_${jcs_drv_part_1_str}_${jcs_drv_part_2_str}_${jcs_drv_ampl_str}";
							mkdir "$data_path/main_${sys_id}_${task_id}_${prop_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}/drv_${jcs_drv_part_1_str}_${jcs_drv_part_2_str}_${jcs_drv_ampl_str}/prm_${jcs_prm_alpha_str}";
							mkdir "$data_path/main_${sys_id}_${task_id}_${prop_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}/drv_${jcs_drv_part_1_str}_${jcs_drv_part_2_str}_${jcs_drv_ampl_str}/prm_${jcs_prm_alpha_str}/start_${start_type}_${start_state}";
							
						}
						else
						{
							mkdir "$data_path";
							mkdir "$data_path/main_${sys_id}_${task_id}_${prop_id}";
							mkdir "$data_path/main_${sys_id}_${task_id}_${prop_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}";
							mkdir "$data_path/main_${sys_id}_${task_id}_${prop_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/N_${N}";
							mkdir "$data_path/main_${sys_id}_${task_id}_${prop_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}";
							mkdir "$data_path/main_${sys_id}_${task_id}_${prop_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}/drv_${jcs_drv_part_1_str}_${jcs_drv_part_2_str}_${jcs_drv_ampl_str}";
							mkdir "$data_path/main_${sys_id}_${task_id}_${prop_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}/drv_${jcs_drv_part_1_str}_${jcs_drv_part_2_str}_${jcs_drv_ampl_str}/prm_${jcs_prm_alpha_str}";
							mkdir "$data_path/main_${sys_id}_${task_id}_${prop_id}/run_${ex_deep}_${rk_ns}_${num_tp_periods}_${num_obs_periods}/N_${N}/diss_${diss_type}_${diss_gamma_str}_${diss_phase_str}/drv_${jcs_drv_part_1_str}_${jcs_drv_part_2_str}_${jcs_drv_ampl_str}/prm_${jcs_prm_alpha_str}/start_${start_type}_${start_state}";
						}
					}
					
					for($val = $start; $val < $finish; $val += $step)
					{
						if ($sys_id == 0)
						{
							if ($task_id == 0)
							{
								$exp{ForderName_sys_0_sys_0_task_0($i)} = $val;
							}
							elsif ($task_id == 2)
							{
								$exp{ForderName_sys_0_sys_0_task_2($i)} = $val;
							}
							else
							{
								$exp{ForderName_sys_0($i)} = $val;
							}
						}
						elsif ($sys_id == 1)
						{
							if ($task_id == 1)
							{
								$exp{ForderName_sys_1($i)} = $val;
							}
							else
							{
								$exp{ForderName_sys_1($i)} = $val;
							}	
						}
						
						$i++;
					}

					for($i = $start; $i < $finish; $i += $step)
					{
						if ($sys_id == 0)
						{
							if ($task_id == 0)
							{
								$key = ForderName_sys_0_sys_0_task_0($i);
							}
							elsif ($task_id == 2)
							{
								$key = ForderName_sys_0_sys_0_task_2($i); 
							}
							else
							{
								$key = ForderName_sys_0($i); 
							}
						}
						elsif ($sys_id == 1)
						{
							if ($task_id == 1)
							{
								$key = ForderName_sys_1($i);
							}
							else
							{
								$key = ForderName_sys_1($i); 
							}
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
						print WF "seed $i \n";
						print WF "mns $mns \n";
						print WF "num_threads $num_threads \n";
						print WF "num_trajectories $num_trajectories \n";
						print WF "num_tp_periods $num_tp_periods \n";
						print WF "num_obs_periods $num_obs_periods \n";
						print WF "ex_deep $ex_deep \n";
						print WF "rk_ns $rk_ns \n";
						close WF;
						
						open( WF,">$key/params.txt"); 
						print WF "lpn_type $lpn_type \n"; 
						print WF "lpn_eps $lpn_eps \n";
						print WF "lpn_eps_change $lpn_eps_change \n";
						print WF "lpn_delta_s_high $lpn_delta_s_high \n";
						print WF "lpn_delta_s_low $lpn_delta_s_low \n";
						print WF "lpn_delta_f_high $lpn_delta_f_high \n";
						print WF "lpn_delta_f_low $lpn_delta_f_low \n";
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
						print WF "dimer_drv_type $dimer_drv_type \n";
						print WF "dimer_drv_ampl $dimer_drv_ampl \n";
						print WF "dimer_drv_freq $dimer_drv_freq \n";
						print WF "dimer_drv_phase $dimer_drv_phase \n";
						print WF "dimer_prm_E $dimer_prm_E \n";
						print WF "dimer_prm_U $dimer_prm_U \n";
						print WF "dimer_prm_J $dimer_prm_J \n";
						print WF "jcs_drv_part_1 $jcs_drv_part_1 \n";
						print WF "jcs_drv_part_2 $jcs_drv_part_2 \n";
						print WF "jcs_drv_ampl $jcs_drv_ampl \n";
						print WF "jcs_prm_alpha $jcs_prm_alpha \n";
						print WF "start_type $start_type \n";
						print WF "start_state $start_state \n";
						print WF "cd_dim $cd_dim \n";
						print WF "cd_eps $cd_eps \n";
						print WF "deep_num_steps $deep_num_steps \n";
						print WF "deep_dump $deep_dump \n";
						print WF "jump $jump \n";
						close WF;
						
						
						if ($sys_id == 0)
						{
							if ($task_id == 0)
							{
								$test_file = sprintf('%s/mean_config(%d_%d_%d)_rnd(%d_%d)_N(%d)_diss(%d_%0.4f_%0.4f)_drv(%d_%0.4f_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d).txt', $key, $sys_id, $task_id, $prop_id, $i, $mns, $N, $diss_type, $diss_gamma, $diss_phase, $dimer_drv_type, $dimer_drv_ampl, $dimer_drv_freq, $dimer_drv_phase, $dimer_prm_E, $dimer_prm_U, $dimer_prm_J, $start_type, $start_state);
							}
							elsif ($task_id == 2)
							{
								$test_file = sprintf('%s/ci_config(%d_%d_%d)_rnd(%d_%d)_N(%d)_diss(%d_%0.4f_%0.4f)_drv(%d_%0.4f_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d).txt', $key, $sys_id, $task_id, $prop_id, $i, $mns, $N, $diss_type, $diss_gamma, $diss_phase, $dimer_drv_type, $dimer_drv_ampl, $dimer_drv_freq, $dimer_drv_phase, $dimer_prm_E, $dimer_prm_U, $dimer_prm_J, $start_type, $start_state);
							}
							else
							{
								$test_file = sprintf('%s/adr_avg_config(%d_%d_%d)_rnd(%d_%d)_N(%d)_diss(%d_%0.4f_%0.4f)_drv(%d_%0.4f_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d).txt', $key, $sys_id, $task_id, $prop_id, $i, $mns, $N, $diss_type, $diss_gamma, $diss_phase, $dimer_drv_type, $dimer_drv_ampl, $dimer_drv_freq, $dimer_drv_phase, $dimer_prm_E, $dimer_prm_U, $dimer_prm_J, $start_type, $start_state);
							}
						}
						elsif ($sys_id == 1)
						{
							if ($task_id == 1)
							{
								$test_file = sprintf('%s/spec_config(%d_%d_%d)_rnd(%d_%d)_N(%d)_diss(%d_%0.4f_%0.4f)_drv(%0.4f_%0.4f_%0.4f)_prm(%0.4f)_start(%d_%d).txt', $key, $sys_id, $task_id, $prop_id, $i, $mns, $N, $diss_type, $diss_gamma, $diss_phase, $jcs_drv_part_1, $jcs_drv_part_2, $jcs_drv_ampl, $jcs_prm_alpha, $start_type, $start_state);
							}
							else
							{
								$test_file = sprintf('%s/spec_config(%d_%d_%d)_rnd(%d_%d)_N(%d)_diss(%d_%0.4f_%0.4f)_drv(%0.4f_%0.4f_%0.4f)_prm(%0.4f)_start(%d_%d).txt', $key, $sys_id, $task_id, $prop_id, $i, $mns, $N, $diss_type, $diss_gamma, $diss_phase, $jcs_drv_part_1, $jcs_drv_part_2, $jcs_drv_ampl, $jcs_prm_alpha, $start_type, $start_state);
							}
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
}
