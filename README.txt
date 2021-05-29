
Here are some selected data and scripts for the work summarized in the paper:

"Revisiting the neutral dynamics derived limiting guanine-cytosine content
using the human de novo point mutation data",

Wentian Li, Yannis Almirantis, Astero Provata

May 2021

***************************************************************************
data sheet denovo-db.ALL.SHORT.SORTED.SNP.COM1.COM2.much.expanded.tsv.gz
***************************************************************************

This file is created by these steps:

 * download the two files from https://denovo-db.gs.washington.edu/denovo-db/
denovo-db.ssc-samples.variants.tsv.gz and denovo-db.non-ssc-samples.variants.tsv.gz

 * pick up only a few relevant columns from these two files
(these columns are: 1: SampleID, 2:StudyName, 7:PrimaryPhenotype, 9: Chr, 10: Position, 11: Variant,
21: Gene, 22: FunctionClass, 29:CaddScore )

> zcat denovo-db.non-ssc-samples.variants.tsv.gz |grep -v version| cut -f 1,2,7,9,10,11,21,22,29 > denovo-db.ALL.SHORT.tsv
> zcat denovo-db.ssc-samples.variants.tsv.gz |grep -v version|grep -v PubmedID| cut -f 1,2,7,9,10,11,21,22,29 >> denovo-db.ALL.SHORT.tsv

 one more column is added to indicate if the mutation is from SSC (Simons Simplex Collection) or not.

 * sort by chromosome and chromosome positions (the coordinate is in human reference genome version GRCh37/hg19)

> cat denovo-db.ALL.SHORT.tsv | sort -k4,4V -k5,5n > denovo-db.ALL.SHORT.SORTED.tsv

  * filter out SNP only (i.e., no indels)

>head -1 denovo-db.ALL.SHORT.SORTED.tsv > denovo-db.ALL.SHORT.SORTED.SNP.tsv
>cat denovo-db.ALL.SHORT.SORTED.tsv | awk '{ if(length($6) ==3) print $0 }'  >> denovo-db.ALL.SHORT.SORTED.SNP.tsv

  * one chromosome position may appear in the file multiple times, if that position is within a gene and the gene
has multiple transcripts. We use a perl script (not included) to combine multiple lines corresponding to one position.
the resulting file is denovo-db.ALL.SHORT.SORTED.SNP.COM1.tsv

  * a mutation at the same position, with the same target(after-mutation) allele, may occur in multiple persons
in multiple studies. we use another perl script (not included) to combine multiple occurrence of the same mutation.
note that if the target (after-mutation) allele is different, we have the situation of tri-allelic SNP. then we
keep two lines. the resulting file is denovo-db.ALL.SHORT.SORTED.SNP.COM1.COM2.tsv  

   one more column is added to indicate the number of times a de novo mutation occurs in the dataset.

  * finally, we use the GRCh37/hg19 human reference genome to annotate each chromosome position at which a de novo
mutation occurred in denovo-db data. the base in the denovo-db should be consistent with the base in the
hg19 reference genome, otherwise that mutation is removed.

  before the run, these columns exist

  1 SampleID
  2 StudyName
  3 PrimaryPhenotype
  4 Chr
  5 Position
  6 Variant
  7 Gene
  8 FunctionClass
  9 CaddScore

  10 SCC
  11 Count

  after the run, these extra columns are added
  12 triplet context of the mutation (the position where the mutation occurs is in the middle)
  13 the triplet after the mutation
  
  seven columns each for window size 1kb, 2kb, 5kb, 10kb, and 20kb
  14-20 window size (usually it's 1001), number of S (G or C), number of uppercase (unique, non-repetitive,
bases according to RepeatMasker), number of uppercase G or C, number of lowercase (repetitive sequence, 
transposon, according to RepeatMasker), number of lowercase g or c,  number of CpG (ignoring uppercase/lowercase
difference)
  21-27 same for 2kb window
  28-34 same for 5kb window
  25-44 same for 10kb window
  41-48 same for 20kb window

 
***************************************************************************
perl script: count-mut-within-GC-range.pl
***************************************************************************

usage example:
> mutation-count-within-GC-range.pl intergenic control 0.3 0.4 1k

category gc0 gc1 n n(WS)  n(WSp) n(SW) n(SSp) n(SpW) n(SpS)  nAT nCG nCG2
intergenic 0.3 0.4 19304 6080 1151 7215 184 1646 55 8822 8772 1710

the script will count the number of mutations in each mutation type (weak-> strong, or, A/T -> G/C, 
weak -> G or C within CpG, G/C -> A/T, etc.) for (in the above example) intergenic mutations,
control samples only, and the mutation position's 1kb neighborhood has a GC% between 0.3 and 0.4.


***************************************************************************
result file: table2a-2kb-window-intergenic-control.input
***************************************************************************

this file is produced using the above perl script -- identical to Table 2(i) in the paper (first 8 coolumns).

 
***************************************************************************
result file: table2b-2kb-intergenic.input
***************************************************************************

this file is produced by other scripts (not included) which takes 2kb windows from
human genome's intergenic regions to get GC%s, then from these GC%, mean and median
of GC% are obtained for each GC content quantile.  same calculation is carried
out for CpG/CG. this file is identical to Table 2 (ii) in the paper (first 3 columns).


***************************************************************************
R script:
***************************************************************************







