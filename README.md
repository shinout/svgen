SVGenerator 0.1.0
==========
generate DNA Structural Variation from FASTA format (Node.js)

Change Log
----------------
* [0.1.0]: release

Overview
----------------
### Installation ###
    git clone git://github.com/shinout/svgen.git

### Usage ###
#### write JavaScript ####
    var SVGen = require('/path/to/SVGenerator.js');
    var sg = new SVGen({
      path: '/path/to/single_fasta.fa',
      chrom: 'chr11' // must be compatible with the fasta file
    });

    sg.registerDel(54001, 320); // deletion at 54001  (length: 320bp)
    sg.registerIns(123499, 200); // insertion at left side of 123499 with a random sequence (length: 200bp)
    sg.registerIns(123499, 7, 'ACCTGTA'); // you can insert your original fragment
    sg.registerInv(9876543, 624); // inversion at 9876543 with a random sequence (624bp)
    sg.registerDup(34091822, 30); // inversion at 34091822

    sg.genotype(); // then result will come up to your STDOUT
#### command line ####
    node SVGenerator.js sv.bed --chrom chr11 --name sv_chr11 sample.bed sample.fasta  > chr11_sv.fasta

[usage]
    svgen <bed file> <fasta file>
[options]
    --test  dryrun. default: null
    --nonstop  even if there's no registered SVs, execute making a new sequence. defualt: null
    --chrom  sequence id in the given file to make SVs and SNPs, default: all (all sequences in the given file)
    --name | -n  new chromosome name with SVs and SNPs, default: (the same as original sequences)
    --json <json file>   fasta summary file to shortcut calculation.
[bed file columns]
    rname start-position  end-position  SVtype(DEL|INS|INV|DUP|SNP) length
