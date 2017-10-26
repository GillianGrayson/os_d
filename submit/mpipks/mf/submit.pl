use File::Copy;
use Data::Dumper;
use Cwd;
use Math::Trig;
use POSIX;
$dir = getcwd;

$data_path = "/data/biophys/yusipov/os_d/mf_results";

$PI = 3.1415926535897932384626433832795;

for($curr_U = 2.8; $curr_U <= 2.8000001; $curr_U += 0.03)
{
	print "curr_U = $curr_U\n";
	
	$eps_exp_shift = 0.1;
	$eps_start = 1.0e-8;
	$eps_num = 81;
	
	for($curr_cd_eps_id = 0; $curr_cd_eps_id < $eps_num; $curr_cd_eps_id += 1)
	{	
		$eps_mult = 10.0 ** ($eps_exp_shift * $curr_cd_eps_id);
	
		$curr_cd_eps = $eps_start * $eps_mult;
	
		print "curr_cd_eps = $curr_cd_eps\n";
		
		$task = 2;
		$path = "";
		$U_start = $curr_U;
		$U_shift = 0.03;
		$U_num = 1;
		$cd_eps = $curr_cd_eps;
		$cd_eps_num = 1;
		$cd_eps_ndpd = 10;
		$seed_start = 0;
		$seed_num = 1;
		$max_num_seeds = 1000000;
		$mt = 0;
		$num_steps = 10000;
		$npt = 2000;
		$np = 100000;
		$E = 1.0;
		$A = 1.5; 
		$omega = 1.0;
		$phase = 0.0;
		$gamma = 0.1;
		$J = 1.0;
		$cd_dim = 4;
		
		$omega_str = sprintf("%.4f", $omega);
		$phase_str = sprintf("%.4f", $phase);
		$gamma_str = sprintf("%.4f", $gamma);
		$J_str = sprintf("%.4f", $J);
		$E_str = sprintf("%.4f", $E);
		$A_str = sprintf("%.4f", $A);
		$U_str = sprintf("%.4f", $U_start);
		
		$eps_str = sprintf("%.10f", $cd_eps);
		
		$start = 1;
		$finish = 2;
		%exp = ();
		$i = 0;
		$i = $start;
		
		sub ForderName{
			$key_str = $_[0];
			
			return  "$data_path/task_${task}/mt_${mt}/omega_${omega_str}/phase_${phase_str}/g_${gamma_str}/J_${J_str}/E_${E_str}/A_${A_str}/U_${U_str}/eps_${eps_str}/m_${cd_dim}/seed_${key_str}";
		}
		
		mkdir "$data_path/task_${task}";
		mkdir "$data_path/task_${task}/mt_${mt}";
		mkdir "$data_path/task_${task}/mt_${mt}/omega_${omega_str}";
		mkdir "$data_path/task_${task}/mt_${mt}/omega_${omega_str}/phase_${phase_str}";
		mkdir "$data_path/task_${task}/mt_${mt}/omega_${omega_str}/phase_${phase_str}/g_${gamma_str}/";
		mkdir "$data_path/task_${task}/mt_${mt}/omega_${omega_str}/phase_${phase_str}/g_${gamma_str}/J_${J_str}";
		mkdir "$data_path/task_${task}/mt_${mt}/omega_${omega_str}/phase_${phase_str}/g_${gamma_str}/J_${J_str}/E_${E_str}";
		mkdir "$data_path/task_${task}/mt_${mt}/omega_${omega_str}/phase_${phase_str}/g_${gamma_str}/J_${J_str}/E_${E_str}/A_${A_str}";
		mkdir "$data_path/task_${task}/mt_${mt}/omega_${omega_str}/phase_${phase_str}/g_${gamma_str}/J_${J_str}/E_${E_str}/A_${A_str}/U_${U_str}";
		mkdir "$data_path/task_${task}/mt_${mt}/omega_${omega_str}/phase_${phase_str}/g_${gamma_str}/J_${J_str}/E_${E_str}/A_${A_str}/U_${U_str}/eps_${eps_str}";
		mkdir "$data_path/task_${task}/mt_${mt}/omega_${omega_str}/phase_${phase_str}/g_${gamma_str}/J_${J_str}/E_${E_str}/A_${A_str}/U_${U_str}/eps_${eps_str}/m_${cd_dim}";

		for($val = $start; $val < $finish; $val+=1)
		{
			$exp{ForderName($i)} = $val;
			$i++;
		}

		for($i = $start; $i < $finish; $i++)
		{
			$key = ForderName($i);    
			mkdir "$key";
			
			$seed_start = $i;
				
			print "$seed_start \n";
				
			open( WF,">$key/config.txt"); 
			print WF "task $task \n"; 
			print WF "path $path \n";
			print WF "U_start $U_start \n";
			print WF "U_shift $U_shift \n";
			print WF "U_num $U_num \n";
			print WF "cd_eps $cd_eps \n";
			print WF "cd_eps_num $cd_eps_num \n";
			print WF "cd_eps_ndpd $cd_eps_ndpd \n";
			print WF "seed_start $seed_start \n";
			print WF "seed_num $seed_num \n";
			print WF "max_num_seeds $max_num_seeds \n";
			print WF "mt $mt \n";
			print WF "num_steps $num_steps \n";
			print WF "npt $npt \n";
			print WF "np $np \n";
			print WF "E $E \n";
			print WF "A $A \n";
			print WF "omega $omega \n";
			print WF "phase $phase \n";
			print WF "gamma $gamma \n";
			print WF "J $J \n";
			print WF "cd_dim $cd_dim \n";
			close WF;

			$test_file = sprintf('%s/exps_lpn_mt(%d)_omega(%0.4f)_phase(%0.4f)_g(%0.4f)_J(%0.4f)_E(%0.4f)_A(%0.4f)_U(%0.4f)_seed(%d).txt', $key, $mt, $omega, $phase, $gamma, $J, $E, $A, $U_start, $i);
			
			if ($task == 2)
			{
				$test_file = sprintf('%s/ci_eps(%0.10f)_m(%d)_mt(%d)_omega(%0.4f)_phase(%0.4f)_g(%0.4f)_J(%0.4f)_E(%0.4f)_A(%0.4f)_U(%0.4f)_seed(%d).txt', $key, $cd_eps, $cd_dim, $mt, $omega, $phase, $gamma, $J, $E, $A, $U_start, $i);
			}
				
			unless (-e "$test_file")
			{	
				print "qsub -wd $dir run_1_thread.sh $key \n";
				system "qsub -wd $dir run_1_thread.sh $key";
			}
		}
	}
}
