use File::Copy;
use Data::Dumper;
use Cwd;
$dir = getcwd;

 
for($val = 34035119; ($val <= 34035173); $val+=1)
{
	system "scancel $val";
}