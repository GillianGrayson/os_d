use File::Copy;
use Data::Dumper;
use Cwd;
$dir = getcwd;


$id_start = 2115810;
$num_id = 1000;
 

for($val = 2116510; ($val <= 2117020); $val+=1)
{
	system "qdel $val";
}