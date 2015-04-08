#!/user/bin/perl
use warnings;
use strict;
use Getopt::Long;
sub Usage{
    die<< 'EOT';
    Usage: command
           --help        -h       help
           --input       -i       input file name

EOT
}

#main(){
my ($opt_h, $opt_i, $opt_o);

GetOptions("h|help" => \$opt_h,
	   "i|input=s" => \$opt_i,
	   )||Usage();

Usage() if $opt_h;

if(! defined $opt_i){
    Usage();
    exit;
}

for (my $i = 1; $i <23; $i++)
{
  system("plink --bfile $opt_i --chr $i --recode --out dataChr$i");
  system("plink --file dataChr$i --freq --out Chr$i")
}

