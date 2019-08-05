use File::Copy;
use Data::Dumper;
use Cwd;
$dir = getcwd;

for($val = 12574243; ($val <= 12577072); $val+=1)
{
	system "qdel $val";
}