if(@ARGV!=2) {die "Usage: <Images> <scalingfactor> \n"}
#print OUT @ARGV, "\n";

@files = glob($ARGV[0]);
$scale = $ARGV[1];
#print "Here \n";

$dirout = "/Volumes/DANIEL/RFOB/RESULTS/";

#sleep(60);
my $sttime = time;

foreach $file (@files)
{	
	print $file, "\n";
	@ss1 = split(/\//, $file);
	$ss2 = $ss1[-1];
	@ss = split(/\./, $ss2);
	$image2save = $dirout.$ss[0];
	$table2save = $dirout.$ss[0].".txt";
	
	$command = "/Applications/MATLAB_R2015a.app/bin/matlab -nodisplay -nojvm -r 'image_processing $file $image2save $table2save $scale; exit;'";
	system($command);
}

print "Total elapsed time... ", (time - $sttime)/60.0, " minuts \n";

