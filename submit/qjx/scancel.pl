use File::Copy;
use Data::Dumper;
use Cwd;
$dir = getcwd;

 

for($val = 5255488; ($val <= 5256130); $val+=1)
{
	system "scancel $val";
}