#########################################################################
#	File Name: BaseCall.pl
#	> Author: QiangGao
#	> Mail: qgao@genetics.ac.cn 
#	Created Time: Sun 15 Oct 2017 04:41:41 PM CST
#########################################################################

#!/usr/bin/perl -w
use strict;
my $QUALCUTOFF=30;
my ($R1,$gene,$out)=@ARGV;
my $name;
if($R1=~/extendedFrags/){
 ($name)=$R1=~/.*\/(.*?)\.extendedFrags/;
}else{
  ($name)=$R1=~/.*\/(.*?)\/R/;
}
if($out){
	$name=$name;
}
my $database="/public-supool/home/gaolab/gene/gene2.fa";
my %database;
$gene=uc($gene);
open(IN,"$database") or die "$database is not exists\n";
while(<IN>){
	if($_=~/^>/){
		chomp $_;
		my $name=$_;
		$name=~s/>//;
		$name=~s/\-//;
		$name=~s/\s+//g;
		$name=uc($name);
		my $L=<IN>;
		$L=~s/\s+//g;
		chomp $L;
		my $target=<IN>;
		$target=~s/\s+//g;
		chomp $target;
		my $R=<IN>;
		$R=~s/\s+//g;
		chomp $R;
		$database{$name}{L}=uc($L);
		$database{$name}{TARGET}=uc($target);
		$database{$name}{R}=uc($R);
		#print "$name\tR=$R\tt=$target\tL=$L\n";
	}
}
close IN;
if(!(-e "$R1") or !(exists $database{$gene})){
	print "$R1 is not exists or $gene is not in database;\n";
	exit;
}
my $R2=$R1;
$R2=~s/R1\.fastq/R2\.fastq/ if($R2=~/R1.fastq/);
open(IN,"$R1");
my $count=0;
my $passqual=0;
my %seq;
while($a=<IN>){
	$a=~s/\s+//g;
	my $seq=<IN>;
	my $b=<IN>;
	my $qual=<IN>;
	my ($flag,$seq)=qual($seq,$qual);
	if ($flag==1){
		$passqual+=1;
		$seq{$a}=$seq;
	}
	$count+=1;
}
close IN;
#if(-e "$R2" and $R2=~/R2/){
#print "open R2=$R2\n";
#open(IN,"$R2");
#while(my $a=<IN>){
#	$a=~s/\s+//g;
#	my $seq=<IN>;
#	my $b=<IN>;
#	my $qual=<IN>;
#	my $seq2=reverse $seq;
#	$seq2 =~tr/ACGTacgt/TGCAtgca/;
#	my ($flag,$seq)=qual($seq2,$qual);
#	if ($flag==1){
#		$passqual+=1;
#		$seq{$a}=$seq;
#	}
#    $count+=1;
#}
#close IN;
#}
my @name=sort keys %seq;
my $big=@name;
my $in=0;
my $del=0;
my $edit=0;
my $noedit=0;
my $all=0;
my %change;
my %site;
my %site2;
my %all;
my @base=("A","T","C","G");
my %indel;
my %indelT;
foreach(@name){
	my $seq=$seq{$_};
	my @R=split($database{$gene}{L},$seq);
	if(@R!=2){
		my $seq2=reverse $seq;
		$seq2 =~tr/ACGTacgt/TGCAtgca/;
		@R=split($database{$gene}{L},$seq2);
	}
	next if(@R!=2);
	my @tar=split($database{$gene}{R},$R[1]);
	next if(@tar!=2);
	$all+=1;
	if(length $tar[0]==length $database{$gene}{TARGET}){
		$all{$tar[0]}+=1;
		if($tar[0]=~/$database{$gene}{TARGET}/){
			$noedit+=1;
			my ($e)=compare($tar[0],$database{$gene}{TARGET});
		}else{
			my ($e)=compare($tar[0],$database{$gene}{TARGET});
			$edit+=1 if($e>0);
		}
	}elsif(length $tar[0] > length $database{$gene}{TARGET}){
		$indel{$tar[0]}+=1;
		$indelT{$tar[0]}="INS";
		$in+=1;
	}elsif(length $tar[0]<length $database{$gene}{TARGET}){
		$indel{$tar[0]}+=1;
		$indelT{$tar[0]}="DEL";
		$del+=1;
	}else{
		print "ERRO=$tar[0]\n";
	}
}
print "$out\t$count\t$passqual\t$all\t$noedit\t$edit\t$in\t$del\n";

open(OUT,">$out.detil.csv");
print OUT "Key,Num\n";
foreach(sort {$all{$b}<=>$all{$a}} keys %all){
	print OUT "$_,$all{$_}\n";
}
close OUT;
open(OUT,">$out.Indeldetil.csv");
print OUT "Key,Type,Num\n";
foreach(sort {$indel{$b}<=>$indel{$a}} keys %indel){
	print OUT "$_,$indelT{$_},$indel{$_}\n";
}
close OUT;
###################20180905####################
my %allhash;
my %allhashtype;
foreach my $key (keys %all){
	$allhash{$key}=$all{$key};
	$allhashtype{$key}="SNP";
}
foreach my $key (keys %indel){
	$allhash{$key}=$all{$key};
	$allhashtype{$key}=$indelT{$key};
	
}
open(OUT,">$out.AllDetail.csv");
print OUT "Name,Key,Type,Num,AllCount,Ratio\n";
my $lastcount;
my $CUTOFF=0.01;
foreach(sort {$allhash{$b}<=>$allhash{$a}} keys %allhash){
	if($allhash{$_}/$passqual>$CUTOFF and $_!~/N/){
		$lastcount+=$allhash{$_};
	}
}
foreach(sort {$allhash{$b}<=>$allhash{$a}} keys %allhash){
	if($allhash{$_}/$passqual>$CUTOFF and $_!~/N/){
		my $p=int(($allhash{$_}/$lastcount)*10000)/100;
		#my $p=$allhash{$_}/$lastcount;
		#$count
		#$passqual
		print OUT "$out,$_,$allhashtype{$_},$allhash{$_},$lastcount,$p\n";
	}
}
close OUT;

##############################################
open(OUT,">$out.indel.csv");
print OUT "$out\n";
print OUT "Insert,$in\n";
print OUT "Deletion,$del\n";
print OUT "All,".($in+$del)."\n";
close OUT;


my @base=("A","T","C","G","N");
open(OUT,">$out.Count.csv");
print OUT "$name\n";
print OUT "Base";
my @ss=split("",$database{$gene}{TARGET});
for(my $i=0;$i<(length $database{$gene}{TARGET});$i++){
	print OUT ",".($i+1);
}
print OUT "\n";
for(my $i=0;$i<(length $database{$gene}{TARGET});$i++){
	print OUT ",$ss[$i]";
}
print OUT "\n";
foreach my $base (@base){
	print OUT "$base";
	for(my $i=0;$i<(length $database{$gene}{TARGET});$i++){
		if(exists $site2{$base}{$i}){
			print OUT ",$site2{$base}{$i}";
		}else{
			print OUT ",0";
		}
	}
	print OUT "\n";
}
close OUT;

my @base=("A","T","C","G");
open(OUT,">$out.Percent.csv");
print OUT "$name\n";
print OUT "Base";
for(my $i=0;$i<(length $database{$gene}{TARGET});$i++){
	print OUT ",".($i+1);
}
print OUT "\n";
for(my $i=0;$i<(length $database{$gene}{TARGET});$i++){
	print OUT ",$ss[$i]";
}
print OUT "\n";
foreach my $base (@base){
	print OUT "$base";
	for(my $i=0;$i<(length $database{$gene}{TARGET});$i++){
		if(exists $site2{$base}{$i}){
			print OUT ",".int($site2{$base}{$i}/$site2{$i}{all}*10000)/100;
		}else{
			print OUT ",0";
		}
	}
	print OUT "\n";
}
close OUT;



sub qual(){
	my ($seq,$qual)=@_;
	chomp $qual;
	chomp $seq;
	my @tmp=split("",$qual);
	my $count=0;
	my @seq=split("",$seq);
	for(my $i=0;$i<@tmp;$i++){
		my $q=ord($tmp[$i])-33;
		if($q<$QUALCUTOFF){
			#print "$seq[$i]\t";
			$seq[$i]="N";
			#print "$seq[$i]\n";
			$count+=1;
		}
		
	}
	my $leng=scalar @tmp;
	my $nseq=join("",@seq);
	if($count/$leng>=0.5){
		return (0,$nseq);
	}else{
		return (1,$nseq)
	}
}
sub compare(){
	my ($tar,$org)=@_;
	my $count=0;
	my @tar=split("",$tar);
	my @org=split("",$org);
	for(my $i=0;$i<@tar;$i++){
		#next if($tar[$i]=~/$org[$i]/);
		$site{$i}{$org[$i]}{$tar[$i]}+=1;
		$site2{$i}{'all'}+=1 if($tar[$i]!~/N/);
		$site2{$tar[$i]}{$i}+=1;
		if($org[$i]=~/C/ and $tar[$i]=~/T/){
			$count+=1;
		}
		$change{$org[$i]}{$tar[$i]}+=1;
	}
	return $count;

}
