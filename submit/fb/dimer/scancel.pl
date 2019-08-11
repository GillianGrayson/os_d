use File::Copy;
use Data::Dumper;
use Cwd;
$dir = getcwd;

 

for($val = 12879856; ($val <= 12880000); $val+=1)
{
	system "scancel $val";
}