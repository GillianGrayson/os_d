use File::Copy;
use Data::Dumper;
use Cwd;
$dir = getcwd;

 
for($val = 34199409; ($val <= 34199443); $val+=1)
{
	system "scancel $val";
}