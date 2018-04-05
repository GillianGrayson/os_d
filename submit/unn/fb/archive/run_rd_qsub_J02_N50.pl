use File::Copy;
use Data::Dumper;
use Cwd;
$dir = getcwd;
	
$pef_url = "/scratch/denysov/Qdiss";

$N = 50;
$J = 0.02;
$J_str = 'J02';
$g = 0.1;
$hasDriving = 0;
$CalcEig = 1;
$E0 = 0;

$U_min = 0.0;#/4;
$U_max = 0.5;#/4;
$step  = 0.01;#/4;
$N_T = 100;
$NSTEP = 1000;
$A0 = 0;

%exp = ();
$i = 0;
$Nstart = 0;
$Nfinish = 50;
$i = $Nstart; 

#1 - run
#2 - merge
$task = 2;

sub ForderName{
   $key_str = $_[0];
   return  "N${N}/J${J_str}/t${key_str}N${N}_J${J_str}";
}

#mkdir -p "$pef_url/N${N}";
#print "$pef_url/N${N}/J${J_str}\n\n";
#mkdir "$pef_url";
#mkdir "$pef_url/N${N}";
#mkdir "$pef_url/N${N}/J${J_str}";
mkdir "N${N}";
mkdir "N${N}/J${J_str}";

for($val = $U_min + $i * $step; 
   ($val <= $U_max) &&  ($val < $U_min + $Nfinish * $step); 
    $val+=$step)
{
  $exp{ForderName($i)} = $val;
  $i++;
}

if($task == 1)
{
  for($i = $Nstart; $i < $Nfinish; $i++)
  {
    $key = ForderName($i);    
    mkdir "$key";
    $U = $exp{$key};
    
    open( W,">$key/config.txt");
    print W "N, $N \n";  
    print W "J,  $J \n"; 
    print W "E0, $E0 \n"; 
    print W "U,  $U \n";
    print W "g,  $g \n";
    print W "CalcEig, $CalcEig \n";
    print W "hasDriving, $hasDriving\n";
    print W "A0, $A0 \n";
    print W "w, 1\n";
    print W "N_T, $N_T\n";
    print W "NSTEP, $NSTEP\n";
    close W;

    print $key."\n";
    
    print "qsub -wd $dir script.sh $key $dir\n";
    system "qsub -wd $dir script.sh $key $dir";
  }
}
if($task == 2)
{
  @file_m = ("absRho.txt", "angleRho.txt", "traseRho2.txt");

  foreach $f_w (@file_m)
  {
    open( W,">$f_w");
  
    for($i = $Nstart; $i < $Nfinish; $i++)
    #for($i = $Nfinish - 1; $i >= $Nstart; $i--)
    {
      $key = ForderName($i);
      $U = $exp{$key};
      
      open (R,"./$key/$f_w");
      
      $j = 0;
      while (<R>)
      {
        @arr = split(' ', $_);
        $t = $arr[$j];
        print W "$t ";
        $j++;	
      }
      print W "\n";
      close(R);
    }
    close (W);
  }
  $f_w = "eigs.txt";
  open(W,">$f_w");
  for($i = $Nstart; $i < $Nfinish; $i++)
  #for($i = $Nfinish - 1; $i >= $Nstart; $i--)
  {
    $key = ForderName($i);
    $U = $exp{$key};
    
    open (R,"./$key/$f_w");
    $val = <R>;
    ($re_m, $im_m) = split(' ', $val);
    while (<R>)
    {
      ($re, $im) = split(' ', $_);
      if($re > $re_m) {$re_m = $re;}
      if($im > $im_m) {$im_m = $im;}      
    }
    print W "$re_m\n";
    
    close(R);
  }
  close(W);
}