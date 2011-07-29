SVGenerator 1.1.0
==========
generate DNA Structural Variation from FASTA format (Node.js)

Change Log
----------------
* [0.1.0]: release
* [1.0.0]: support multi fasta
* [1.1.0]: support rnames option

Overview
----------------
### Installation ###
    git clone git://github.com/shinout/svgen.git
    cd svgen
    git submodule init
    git submodule update

### Usage ###
#### write JavaScript ####
    var SVGen = require('/path/to/svgen.js');
    var svgen = new SVGen('/path/to/single_fasta.fa');
    // rname, start, len, type, extra
    svgen.register('chr5', 5400001, 320, 'DEL');         // deletion at chr5:5400001  (length: 320bp)
    svgen.register('chr22', 123499, 200, 'INS', sg.getRandomFragment(200));          // insertion at left side of chr22:123499 with a random sequence (length: 200bp)
    svgen.register('chrX', 9876543, 624, 'INV');         // inversion at chrX:9876543 (length: 625bp)
    svgen.register('chr2', 2286120, 30, 'DUP', 11);         // tandem duplictation at chr2:2286120 (length: 625bp)
    svgen.register('chr9', 34091822, 300, 'TRA', 'chr1:91152048');         // translocation at 34091822 from chr1:91152048 (length:300bp)
    svgen.run(); // then result will come up to your STDOUT

#### command line ####
    # get information sv and snp randomly generated.
    node svbedgen.js sample.fasta > sv.bed

    # get fasta with SV and SNP designated in a bed file
    node svgen.js sv.bed sample.fasta  > chr11_sv.fasta

    **** svbedgen ****
    [synopsis]
      node svbedgen.js <fasta file>
    [options]
      ** length: the lengths of SVs.
      --svlen <int> mean length of SVs. default:1500
      --inslen <int>  mean length of insertions. default: the same as "svlen" option. if 0, then no insertions are registered.
      --dellen <int>  mean length of deletions. default: the same as "svlen" option. if 0, then no deletions are registered.
      --invlen <int>  mean length of inversions. default: the same as "svlen" option. if 0, then no inversions are registered.
      --duplen <int>  mean length of tandem duplications. default: the same as "svlen" option. if 0, then no tandem duplications are registered.
      --tralen <int>  mean length of translocations. default: the same as "svlen" option. if 0, then no translocations are registered.

      ** deviation : standard deviation of each SVs.
      --svdev <int> stddev of SVs. default:300
      --insdev <int>  stddev of insertions. default: the same as "svlen" option
      --deldev <int>  stddev of deletions. default: the same as "svlen" option
      --invdev <int>  stddev of inversions. default: the same as "svlen" option
      --dupdev <int>  stddev of tandem duplications. default: the same as "svdev" option.
      --tradev <int>  stddev of translocations. default: the same as "svlen" option

      ** number: total number of SVs.
      --sv <int>  the number of SVs to register. default:10000. If each SV is specifically designated its number, that configuration priors to this.

      ** rate : the rate of each SV type, determined by the following rule.
      *  let "total" be the sum of each num, i.e. insnum, delnum, invnum, dupnum and tranum.
      *  the rate of each SV is  x / total, where x is one of insnum, delnum, invnum, dupnum or tranum.
      --insnum <int>  for insertions. default: 200.
      --delnum <int>  for deletions. default: 200.
      --invnum <int>  for inversions. default: 200.
      --dupnum <int>  for tandem duplications. default: 200.
      --tranum <int>  for translocations. default: 200.

      ** repeat number of tandem duplication.
      --duprep <int>  mean repeat folds of tandem duplications. default: 40.
      --duprdev <int> stddev of repeat folds of tandem duplications. default: 10.

      ** configuration for SNP.
      * for simplicity, only single nucleotide alteration is supported.
      --snprate reciprocal rate to insert SNP default:10000 (1/10000). if 0, then no SNPs are registered.

      ** configuration for a fasta file.
      --rnames <sequence id1>[,sequence_id2,...]   sequence ids to use. default: null (use all rnames in a fasta file.)
      --json <json file>   fasta summary file to shortcut calculation.

    **** svgen ****
    [usage]
      node svgen.js <bed file> <fasta file>
    [options]
      --rnames|-r <sequence id1>[,sequence_id2,...]   sequence ids to use. default: null (use all rnames in a fasta file.)
      --json|-j <json file>  fasta summary file to shortcut calculation.
    [bed file columns]
      rname start-position  end-position  SVtype(DEL|INS|INV|DUP|TRA|SNP) length  extra-info
