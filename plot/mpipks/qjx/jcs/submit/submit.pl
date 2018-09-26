use File::Copy;
use Data::Dumper;
use Cwd;
use Math::Trig;
$dir = getcwd;

print "qsub -wd $dir run.sh\n";
system "qsub -wd $dir run.sh";
