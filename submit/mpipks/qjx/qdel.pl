use File::Copy;
use Data::Dumper;
use Cwd;
$dir = getcwd;


$id_start = 5017824;
$num_id = 2000;
 

for($val = $id_start; ($val < $id_start + $num_id); $val+=1)
{
	system "qdel $val";
}