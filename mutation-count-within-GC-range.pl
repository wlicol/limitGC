#!/usr/bin/perl

# purpose: count number of de novo mutations in each category (e.g. intergenic), 
# from either control or all samples, within a given GC% range, and from which local
# neighborhood a GC% attached to a location/mutation is used (1kb, 2kb, 5kb, 10kb, 20kb). 
#
# the input file: denovo-db.ALL.SHORT.SORTED.SNP.COM1.COM2.much.expanded.tsv.gz
#
# choice for "function class": intron,intergenic,missense,stream,synonymous,UTR,non-coding,splice-,stop
# (note: "splice-" needs a "-")
#
# choice for sample phenotype: control, all
#
# choice for local neighborhood size where the GC% is calculated: 1k, 2k, 5k, 10k, 20k
#
# bracket for GC% range: lower bound and upper bound.
#



if($#ARGV !=4){
 print "[.pl] intron/intergenic/missense/stream/synonymous/UTR/non-coding/splice-/stop control/all GC0 GC1 [1k/2k/5k/10k/20k]\n";
 exit;
}

$f="denovo-db.ALL.SHORT.SORTED.SNP.COM1.COM2.much.expanded.tsv.gz";
if($f =~ m/\.gz/){ open(IN,"gunzip -c $f | "); }else{ open(IN,"< $f"); }

 
$category=$ARGV[0];
$sample =$ARGV[1];
if($sample eq "all"){ $sample="[a-z]";}

$GC0=$ARGV[2];
$GC1=$ARGV[3];

# specific for this denovo-db.ALL.SHORT.SORTED.SNP.COM1.COM2.expanded.tsv.gz file
$wchoice=$ARGV[4];
if($wchoice eq "1k"){	# 7 columns from (0-init) 13-19
 $I=13;
 $wsize=1000;
}elsif($wchoice eq "2k"){	# 20-26
 $I=20;
 $wsize=2000;
}elsif($wchoice eq "5k"){	# 27-33
 $I=27;
 $wsize=5000;
}elsif($wchoice eq "10k"){	# 34-40
 $I=34;
 $wsize=10000;
}elsif($wchoice eq "20k"){	# 41-47 (last column)
 $I=41;
 $wsize=20000;
}

$nww=$nws=$nws2=0;
$nsw=$nss2=$nss=0;
$ns2w=$ns2s=$ns2s2=0;
$ntotal=0;
$nat=$ncg=$ncg2=0;

$line=0;
while($_=<IN>){ if($line >0){	# there is a header
 chop;
 @a=split(/\t/, $_);

# 0 SampleID, 1 StudyName, 2 PrimaryPhenotype, <-
#  3 Chr, 4 Position, 5 Variant ++, 6 Gene
# 7 FunctionClass <-
# 8 CaddScore, 9 SSC, 10 Count
# 11 tri <-
# 12 target tri
# 13-19 1k, 20-26 2k, 27-33 5k, 34-40 10k, 41-47 20k

 # num_GC / num_typed
 $tmpgc= $a[$I+1]/$a[$I];

# if(($a[7] =~ m/$category/ ) & ($a[2] =~ m/$sample/) & ($tmpgc >= $GC0) & ($tmpgc < $GC1) ){
# one more requirement? require the number of typed should > (%) $wsize (e.g. 50%)
if(($a[7] =~ m/$category/ ) & ($a[2] =~ m/$sample/) & ($tmpgc >= $GC0) & ($tmpgc < $GC1) & ($a[$I] > 0.5*$wsize) ) {

# print "TEST $a[7] $a[2] $a[11] -> $a[12] $a[10]\n";
 $ntotal += $a[10];

 # print "$a[7] $a[2] $a[5] $a[11] $a[12]\n";
 @tri0=split(//, $a[11]);
 $b0 = lc($tri0[1]);
 @tri1=split(//, $a[12]);
 $b1 = lc($tri1[1]);

# w->
 if($b0 eq "a" | $b0 eq "t"){
  # not plus 1, but add a[10]!
  $nat +=$a[10];
  # w->w
  if($b1 eq "a" | $b1 eq "t"){ $nww += $a[10]; }
  # w->s2
  elsif($a[12] =~ m/cg/){ $nws2 += $a[10];}
  # w->s
  elsif($b1 eq "c" | $b1 eq "g"){ $nws += $a[10];}
  else{ print "sth wrong $a[11] -> $a[12] \n"; exit;}

#s2 ->
 }elsif($a[11] =~ m/cg/){
  # not plus 1, but add a[10]!
  $ncg2 +=$a[10];
  # s2->w
  if($b1 eq "a" | $b1 eq "t"){ $ns2w += $a[10]; }
  # s2->s2
  elsif($a[12] =~ m/cg/){ $ns2s2 += $a[10];}
  # s2->s
  elsif($b1 eq "c" | $b1 eq "g"){ $ns2s += $a[10];}
  else{ print "sth wrong $a[11]-> $a[12]\n"; exit;}
# s ->
 }elsif($b0 eq "c" | $b0 eq "g"){
  # not plus 1, but add a[10]!
  $ncg +=$a[10];
  # s -> w
  if($b1 eq "a" | $b1 eq "t"){ $nsw += $a[10]; }
  # s->s2
  elsif($a[12] =~ m/cg/){ $nss2 += $a[10];}
  # s->s
  elsif($b1 eq "c" | $b1 eq "g"){ $nss += $a[10];}
  else{ print "sth wrong $a[11]-> $a[12]\n"; exit;}
  
 }else{
  print "sth wrong $a[11]-> $a[12]\n"; exit;}

} # end of filtering




} $line++; }


print "category gc0 gc1 n n(WS)  n(WSp) n(SW) n(SSp) n(SpW) n(SpS)  nAT nCG nCG2\n";
print "$category $GC0 $GC1 $ntotal $nws $nws2 $nsw $nss2 $ns2w $ns2s $nat $ncg $ncg2\n";



