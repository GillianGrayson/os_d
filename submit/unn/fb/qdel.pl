use File::Copy;
use Data::Dumper;
use Cwd;
$dir = getcwd;


$id_start = 835955;
$num_id = 200;
 

for($val = $id_start; ($val < $id_start + $num_id); $val+=1)
{
	system "qdel $val";
}