use File::Copy;
use Data::Dumper;
use Cwd;
$dir = getcwd;

for($val = 2128211; ($val <= 2129110); $val+=1)
{
	system "qdel $val";
}