#!/bin/bash
thisdir=$(dirname $0)
if [ $# -lt 5 ]; then
  echo "$0 <tsv file> <fasta file> <name prefix> <result dir> <chrom name1> [chrom name2 ...]" >&2
fi

tsv=$1
shift
fasta=$1
shift
name=$1
shift
dir=$1
shift
rnames=$*
#[-c|--chrom <chrom name>] [-n|--name <sv chrom name>] <tsv file> <fasta file>
if [ ! -d $dir ]; then
  echo "$dir : No such directory." >&2
  exit
fi


for chrom in $rnames;
do
  outfile=$dir/$name.$chrom.sv.fasta
  #echo "node $thisdir/SVGenerator.js --chrom $chrom --name $chrom $tsv $fasta > $outfile"
  node $thisdir/SVGenerator.js --nonstop --chrom $chrom --name $chrom $tsv $fasta > $outfile &
  echo $outfile
done
