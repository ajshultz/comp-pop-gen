use strict ; 
use warnings ; 

foreach my $species ( `ls /n/holylfs/LABS/informatics/ashultz/CompPopGen/SPECIES_DATASETS/` ) { 
	chomp $species ; 
	system("mkdir $species") ; 
	system("mkdir ${species}/log") ; 
	system("perl -pi -e 's/SPECIESS/$species/g' < run_qc.sbatch > ${species}/run_qc.sbatch") ; 
	system("cp qc.sh $species/qc.sh") ; 
	chdir("$species/") ; 
	system("sbatch run_qc.sbatch") ; 
	chdir("..") ;
}
