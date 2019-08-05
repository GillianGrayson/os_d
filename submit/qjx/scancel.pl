use File::Copy;
use Data::Dumper;
use Cwd;
$dir = getcwd;

 

for($val = 12616178; ($val <= 12616279); $val+=1)
{
	system "scancel $val";
}