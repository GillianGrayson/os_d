use File::Copy;
use Data::Dumper;
use Cwd;
use Math::Trig;
$dir = getcwd;

$root = "/home/yusipov_i/Work/os_d";

$PI = 3.14159265358979323846;

for($curr_U = 0.04; $curr_U <= 4.00000001; $curr_U += 0.04)
{
	print "curr_U=$curr_U\n";
	
	$N = 100;
	$J = -1; 
	$E0 = 0.0; 
	$U = $curr_U;
	$g = 0.1; 
	$CalcEig = 0; 
	$hasDriving = 1;
	$A0 = -3.4; 
	$w = 1;
	$N_T = 21;
	$NSTEP = 20000;
	
	$seed_begin = 1;
	$seed_end = 2;
	
	%exp = ();
	$i = $seed_begin;

	$J_str = sprintf("%.4f", $J);
	$E0_str = sprintf("%.4f", $E0);
	$U_str = sprintf("%.4f", $U);
	$g_str = sprintf("%.4f", $g);
	$A0_str = sprintf("%.4f", $A0);
	$omega_str = sprintf("%.4f", $w);

	sub ForderName{
		$key_str = $_[0];
		return  "$root/data/N_${N}/E0_${E0_str}/J_${J_str}/U_${U_str}/g_${g_str}/A0_${A0_str}/omega_${omega_str}/seed_${key_str}";
	}

	mkdir "$root/data";
	mkdir "$root/data/N_${N}";
	mkdir "$root/data/N_${N}/E0_${E0_str}";
	mkdir "$root/data/N_${N}/E0_${E0_str}/J_${J_str}/";
	mkdir "$root/data/N_${N}/E0_${E0_str}/J_${J_str}/U_${U_str}";
	mkdir "$root/data/N_${N}/E0_${E0_str}/J_${J_str}/U_${U_str}/g_${g_str}/";
	mkdir "$root/data/N_${N}/E0_${E0_str}/J_${J_str}/U_${U_str}/g_${g_str}/A0_${A0_str}";
	mkdir "$root/data/N_${N}/E0_${E0_str}/J_${J_str}/U_${U_str}/g_${g_str}/A0_${A0_str}/omega_${omega_str}";
	
	for($val = $seed_begin; ($val < $seed_end ); $val+=1)
	{
		$exp{ForderName($i)} = $val;
		$i++;
	}

	for($i = $seed_begin; $i < $seed_end; $i++)
	{
		$key = ForderName($i);    
		mkdir "$key";
		
		$expId = $exp{$key};
		
		open( WF,">$key/config.txt");
		print WF "N, $N \n";  
		print WF "J,  $J \n"; 
		print WF "E0, $E0 \n"; 
		print WF "U,  $U \n";
		print WF "g,  $g \n";
		print WF "CalcEig, $CalcEig \n";
		print WF "hasDriving, $hasDriving\n";
		print WF "A0, $A0 \n";
		print WF "w, $w\n";
		print WF "N_T, $N_T\n";
		print WF "NSTEP, $NSTEP\n";
		close WF;
				
		$test_file = sprintf('%s/purity_avg.txt', $key);
		
		unless (-e "$test_file") 
		{
			print "sbatch run.sh $key \n";
			system "sbatch run.sh $key ";
		}
	}
}
