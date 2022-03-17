#  the input bam without readgroup information #	
### updated 1/9/2018 ###
#!/usr/bin/perl
use strict;
use warnings;
#use POSIX;
use Getopt::Long;
my $version = 0.1;
#color code
my $red = "\e[31m";
my $gray = "\e[37m";
my $yellow = "\e[33m";
my $green = "\e[32m";
my $purple = "\e[35m";
my $cyan = "\e[36m";
my $normal = "\e[0m";
#usage information

(my $usage = <<OUT) =~ s/\t+//g;
This script will do readcount for tumor and normal bam based on an input vcf file for hotspot mutations. 
Pipeline version: $version
$yellow     Usage: perl $0  --rdir --vcf --ref --log --groupname --users --q --step

$normal

<rdir> = full path of the folder holding files for this sequence run (user must provide)
<groupname> = job group name
<users> = user name for job group
<log> = full path of the folder for saving log file; usually upper folder of rdir 
<step> run this pipeline step by step. (user must provide)
<ref> the human reference: 
<q> which queue for submitting job; research-hpc, ding-lab, long (default)

GDC HG38: /storage1/fs1/songcao/Active/Database/hg38_database/GRCh38.d1.vd1/GRCh38.d1.vd1.fa 

<run_folder> = full path of the folder holding files for this sequence run
<step_number> run this pipeline step by step. (running the whole pipeline if step number is 0)

$red      	 [0]  Do bam index if index file is not existing  
$red         	 [1]  Run bamreadcount

$normal
OUT

#die $usage unless @ARGV == 2;
#my ( $run_dir, $step_number ) = @ARGV;
#if ($run_dir =~/(.+)\/$/) {
#    $run_dir = $1;
#}
#die $usage unless ($step_number >=0)&&(($step_number <= 10));
#GENOMEVIP_SCRIPTS=/gscmnt/gc2525/dinglab/rmashl/Software/bin/genomevip
# obtain script path
#my $run_script_path = `dirname $0`;
#__DEFAULT NUMBER OF BINS IE (MUST BE INTEGER)
my $step_number = -1;
my $status_rg = 1;
my $status_rerun=0;
#__HELP (BOOLEAN, DEFAULTS TO NO-HELP)
my $help = 0;

#__FILE NAME (STRING, NO DEFAULT)
my $run_dir="";
my $log_dir="";
my $h38_REF="";
my $q_name="";
my $chr_status=0;
my $f_vcf; 
my $f_vcf_cut;
my $compute_username="";
my $group_name="";

my $status = &GetOptions (
      "step=i" => \$step_number,
      "rdir=s" => \$run_dir,
      "vcf=s" => \$f_vcf,
      "groupname=s" => \$group_name,
      "users=s" => \$compute_username,
      "ref=s"  => \$h38_REF,
      "log=s"  => \$log_dir,
      "q=s" => \$q_name,
      "help" => \$help,
    );

print $group_name,"\n"; 
print $compute_username, "\n"; 

if ($help || $run_dir eq "" || $log_dir eq ""  || $f_vcf eq "" || $group_name eq "" || $compute_username eq "" || $step_number<0 || $step_number>8) {
      print $usage;
      exit;
   }

print "run dir=",$run_dir,"\n";
print "step num=",$step_number,"\n";
print "queue name=",$q_name,"\n";
print "job group=",$group_name,"\n";
print "user group=",$compute_username,"\n";

if($q_name eq "")
{
    $q_name="general";
}

if ($run_dir =~/(.+)\/$/) {
    $run_dir = $1;
}

my $email = "scao\@wustl\.edu";
my $HOME = $ENV{HOME};
my $working_name= (split(/\//,$run_dir))[-2];
my $HOME1=$log_dir;

if (! -d $HOME1)
{
`mkdir $HOME1`;
}

if (! -d $HOME1."/temprc") {
    `mkdir $HOME1"/temprc"`;
}

my $job_files_dir = $HOME1."/temprc";

if (! -d $HOME1."/LSF_DIR_RC") {
    `mkdir $HOME1"/LSF_DIR_RC"`;
}
my $lsf_file_dir = $HOME1."/LSF_DIR_RC";

my $run_script_path =`echo \$PWD`;
chomp $run_script_path;
my $script_dir=$run_script_path;

$run_script_path = "/usr/bin/perl ".$run_script_path."/";
print $run_script_path,"\n";
my $hold_RM_job = "norm";
my $current_job_file = "";#cannot be empty
my $hold_job_file = "";
my $bsub_com = "";
my $sample_full_path = "";
my $sample_name = "";
my $h38_REF_bai=$h38_REF.".fai";
my $samtools="/storage1/fs1/songcao/Active/Software/samtools/1.2/bin";
#my $STRELKA_DIR="/gscmnt/gc2525/dinglab/rmashl/Software/bin/strelka/1.0.14/bin";

my $bamrc="/storage1/fs1/songcao/Active/Software/bam-readcount/0.7.4/bam-readcount";
#my $f_vcf = $script_dir."/hotspot.tcga.driver.chr.pos.tsv";
#my $f_mut_tcga_hotspot = $script_dir."/liftover_hg38.tcga.driver.gene.mut";

$f_vcf_cut = $f_vcf.".cut"; 
open(OUT,">$f_vcf_cut"); 
my %chrpos; 
my $chrpos2; 
foreach my $l (`cat $f_vcf`) 
{
my $ltr=$l;
my @t=split("\t",$ltr);  
$chrpos2=$t[0]."-".$t[1]."-".$t[2];
if(!defined $chrpos{$chrpos2})
{
print OUT $t[0],"\t",$t[1],"\t",$t[2],"\n"; 
}

}

close OUT; 

my $first_line=`head -n 1 $h38_REF`; 

if($first_line=~/^\>chr/) { $chr_status=1; }

opendir(DH, $run_dir) or die "Cannot open dir $run_dir: $!\n";
my @sample_dir_list = readdir DH;
close DH;

if ($step_number < 8) {
    for (my $i=0;$i<@sample_dir_list;$i++) {#use the for loop instead. the foreach loop has some problem to pass the global variable $sample_name to the sub functions
        $sample_name = $sample_dir_list[$i];
        if (!($sample_name =~ /\./ || $sample_name=~/worklog/)) {
            $sample_full_path = $run_dir."/".$sample_name;
            if (-d $sample_full_path) { # is a full path directory containing a sample
                print $yellow, "\nSubmitting jobs for the sample ",$sample_name, "...",$normal, "\n";
                $current_job_file="";
                if($step_number==0)
                {
                   &bsub_bam();
                }elsif ($step_number == 1) {
                    &bsub_rc(1);
		}
		}
	}
}
}   

sub bsub_bam{

    $current_job_file = "j0_bam_".$sample_name.".sh";
    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    if(-e $lsf_out)
    {
    `rm $lsf_out`;
    `rm $lsf_err`;
    `rm $current_job_file`;
    }
    my $IN_bam_T = $sample_full_path."/".$sample_name.".T.bam";
    my $IN_bam_N = $sample_full_path."/".$sample_name.".N.bam";

#    if (! -e $IN_bam_N) {#make sure there is a input fasta file 
#        print $red,  "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
#        print "Warning: Died because there is no input bam file for bwa:\n";
#        print "File $IN_bam_N does not exist!\n";
#        die "Please check command line argument!", $normal, "\n\n";

 #   }
 #   if (! -s $IN_bam_N) {#make sure input fasta file is not empty
 #       print $red, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
 #       die "Warning: Died because $IN_bam_N is empty!", $normal, "\n\n";
 #   }

    open(BAM, ">$job_files_dir/$current_job_file") or die $!;
    print BAM "#!/bin/bash\n";
    print BAM "TBAM=".$sample_full_path."/".$sample_name.".T.bam\n";
    print BAM "NBAM=".$sample_full_path."/".$sample_name.".N.bam\n";
    print BAM "TBAM_bai=".$sample_full_path."/".$sample_name.".T.bam.bai\n";
    print BAM "NBAM_bai=".$sample_full_path."/".$sample_name.".N.bam.bai\n";
    print BAM "export SAMTOOLS_DIR=$samtools\n";
    print BAM "if [  -e \${TBAM} ]\n";
    print BAM "then\n";
    print BAM "if [ ! -e \${TBAM_bai} ]\n";
    print BAM "then\n";
    print BAM "\${SAMTOOLS_DIR}/samtools index \${TBAM}\n";
    print BAM "fi\n";
    print BAM "fi\n";
    print BAM "if [  -e \${NBAM} ]\n";
    print BAM "then\n";
    print BAM "if [ ! -e \${NBAM_bai} ]\n";
    print BAM "then\n";
    print BAM "\${SAMTOOLS_DIR}/samtools index \${NBAM}\n";
    print BAM "fi\n";
    print BAM "fi\n";
    close BAM;

    my $sh_file=$job_files_dir."/".$current_job_file;


    $bsub_com = "LSF_DOCKER_ENTRYPOINT=/bin/bash LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -g /$compute_username/$group_name -q $q_name -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(scao/dailybox)\' -o $lsf_out -e $lsf_err bash $sh_file\n";
    print $bsub_com;
    system ($bsub_com);

}

sub bsub_rc{

    my ($step_by_step) = @_;
    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }

    $current_job_file = "j1_rc_".$sample_name.".sh";
    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    if(-e $lsf_out)
    {
    `rm $lsf_out`;
    `rm $lsf_err`;
    `rm $current_job_file`;
    }
    my $IN_bam_T = $sample_full_path."/".$sample_name.".T.bam";
    my $IN_bam_N = $sample_full_path."/".$sample_name.".N.bam";
    my $f_rc_t_out = $sample_full_path."/".$sample_name.".T.rc.tsv";
    my $f_rc_n_out = $sample_full_path."/".$sample_name.".N.rc.tsv";
    my $f_vaf_t_out = $sample_full_path."/".$sample_name.".T.rc.vaf";
    my $f_vaf_n_out = $sample_full_path."/".$sample_name.".N.rc.vaf";
    open(RC, ">$job_files_dir/$current_job_file") or die $!;
    print RC "#!/bin/bash\n";
    print RC "TBAM=".$sample_full_path."/".$sample_name.".T.bam\n";
    print RC "NBAM=".$sample_full_path."/".$sample_name.".N.bam\n";
    print BAM "if [ -e \${TBAM} ]\n";
    print BAM "then\n";   
    print RC "$bamrc -q 10 -b 10 \${TBAM} -f $h38_REF -l $f_vcf_cut > $f_rc_t_out","\n";  
    print RC "     ".$run_script_path."bamReadcount2vaf.pl -s $sample_name -l $f_vcf $f_rc_t_out > $f_vaf_t_out","\n"; 
    print RC "fi\n";
    print RC "if [ -e \${NBAM} ]\n";	 
    print RC "$bamrc -q 10 -b 10 \${NBAM} -f $h38_REF -l $f_vcf_cut > $f_rc_n_out","\n";
    print RC "     ".$run_script_path."bamReadcount2vaf.pl -s $sample_name -l $f_vcf $f_rc_n_out > $f_vaf_n_out","\n"; 
    print RC "fi\n";  
    close RC;
    my $sh_file=$job_files_dir."/".$current_job_file;


    $bsub_com = "LSF_DOCKER_ENTRYPOINT=/bin/bash LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -g /$compute_username/$group_name -q $q_name -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(scao/dailybox)\' -o $lsf_out -e $lsf_err bash $sh_file\n";
    print $bsub_com;
    system ($bsub_com);

} 
