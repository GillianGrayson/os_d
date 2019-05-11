use File::Copy;
use Data::Dumper;
use Cwd;
$dir = getcwd;

 

for($val = 2125628; ($val <= 2126526); $val+=1)
{
	system "qdel $val";
}