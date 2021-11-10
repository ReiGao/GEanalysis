#V3
#perl $0 

my ($outdir,$barfile,@file)=@ARGV;

#my ($outdir)=$barfile=~/(\S+)\.csv/;
if(!$barfile){
        print "USAGE:\nperl $0 ATCG CG bar.csv /public/source/fastq/R1.fq.gz /public/source/fastq/R1.fq.gz\n";
        exit;
}
print "Outdir:$outdir\n";
#chomp $file;
my $hashout;
open(BAR,"$barfile");
if(-e "$outdir"){
	print "Exists $outdir file\n";
	my $cc=`rm -R $outdir`;
	my $cmd=`mkdir $outdir`;
}else{
	my $cmd=`mkdir $outdir`;
}
my %barcode,%name;
while(<BAR>){
   chomp $_;
   my ($name,$bar,$barR2)=split("\t",$_);
   $bar = uc($bar);
   $barR2 = uc($barR2);
   chomp $bar;
   chomp $barR2;
   $bar=~s/\s//g;
   $barR2=~s/\s//g;
#   $bar.="T";
   #print "$bar\t$name\n";
   $barcode{$bar}{$barR2}=$name;
   $name{$name}=$bar.$barR2;
   
   my $cmd=`mkdir $outdir/$name`;
   	my $R1=$bar.$barR2."R1";
	my $R2=$bar.$barR2."R2";
	$hashout{$R1}=$R1;
	$hashout{$R2}=$R2;
	open($hashout{$R1},">./$outdir/$name/R1.fastq");
	open($hashout{$R2},">./$outdir/$name/R2.fastq");
   
   
}
close BAR;
my $all=0;
my %hash;
foreach my $file(@file){
chomp $file;
my $file2=$file;
$file2=~s/r1/r2/;
print "Deal with $file\n$file2\n.....\n";
open(IN,"gzip -dc $file|") or die "$file not exists!\n";
open(IN2,"gzip -dc $file2|") or die "$file2 not exists!\n";
open(UN1,">>$outdir/R1_unBarcode.fastq");
open(UN2,">>$outdir/R2_unBarcode.fastq");
while(defined(my $a=<IN>) and defined ($b=<IN2>)){
	my $name1=$a;
	my $seq=<IN>;
    my $s1=<IN>;
	my $qu=<IN>;
	#$seq=~s/\.//;
	#$s1=~s/^\.//;
	#$qu=~s/^\.//;
	my $name2=$b;
	my $seq2=<IN2>;
	my $s2=<IN2>;
	my $qu2=<IN2>;
	
	#my $s2=$b.<IN2>.<IN2>.<IN2>;
	$all+=1;
   my $flag=0;
	#print "$.===$seq\n";
      foreach my $bar(keys %barcode){
			my @barR2=keys %{$barcode{$bar}};
			foreach my $bar2(@barR2){
		 #print "$seq$bar\t$barcode{$bar}\n\n";
	       if(($seq =~ /^$bar/ and $seq2=~/^$bar2/)  or  ($seq2=~/^$bar/ and $seq=~/^$bar2/)){
			$flag=1;
			#print "now====$bar\n";
			$hash{$bar.$bar2}+=1;
#			$seq=~s/^$bar//;
#			$qu=substr($qu,(length $bar));
			#open(OUT1,">>$outdir/$barcode{$bar}/R1.fastq");
			#print OUT1 $name1.$seq.$s1.$qu;
			#close OUT1;
			my $out1=$bar.$bar2."R1";
			#printf {$hashout{$out1}} %s,"$name1$seq$s1$qu";
			print {$hashout{$out1}} "$name1$seq$s1$qu";
			my $out2=$bar.$bar2."R2";
#			$seq2=~s/^$bar2//;
#			$qu2=substr($qu2,(length $bar2));
			print {$hashout{$out2}} "$name2$seq2$s2$qu2";
			#open(OUT2,">>$outdir/$barcode{$bar}/R2.fastq");
			#print OUT2 $s2;
			#close OUT2;
			
			last;		
		}
		}
       }
       if($flag==0){
		$hash{UnBarcode}+=1;
		print UN1 $name1.$seq.$s1.$qu;
		print UN2 $name2.$seq2.$s2.$qu2;
	}
	#
}
close IN;
close IN2;
close UN1;
close UN2;
}
open(RE,">$outdir/$barfile.Result.txt");
foreach(sort keys %name){
	print RE "$_\t$name{$_}\t$hash{$name{$_}}\t".(int($hash{$name{$_}}/$all*10000)/100)."%\n";
}
 print RE "UnBarcode\tUnbarcode\t$hash{UnBarcode}\t".(int($hash{UnBarcode}/$all*10000)/100)."%\n";
print RE "All READS\t$all\n";

close RE;
