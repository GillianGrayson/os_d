use File::Copy;
use Data::Dumper;
use Cwd;
$dir = getcwd;

 

for($val = 10314600; ($val <= 10314697); $val+=1)
{
	system "scancel $val";
}