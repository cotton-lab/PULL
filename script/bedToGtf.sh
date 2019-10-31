perl -n -e '
chomp; 
my @f=split("\t",$_);
my @starts=split(",",$f[11]);
my @blocks=split(",",$f[10]);
my $count=$f[9]; 
my $group="gene_id \"$f[3]\"; transcript_id \"$f[3]\";  exon_number $count;";
my $xxx=$f[1]+1;
print "$f[0]\tCufflinks\ttranscript\t$xxx\t$f[2]\t$f[4]\t$f[5]\t.\t$group\n";
for(my $i=0; $i<$count; $i++){
my $start=$f[1]+$starts[$i];
my $end=$start+$blocks[$i];
my $xxx=$start+1;
print "$f[0]\tCufflinks\texon\t$xxx\t$end\t$f[4]\t$f[5]\t.\t$group\n";
}'  $1 

