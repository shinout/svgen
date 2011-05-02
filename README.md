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
    sg.registerIns(123499, 7, 'ACCTGTA'); // you can insert your original flagment
    sg.registerInv(9876543, 624); // inversion at 9876543 with a random sequence (624bp)

    sg.genotype(); // then result will come up to your STDOUT
#### command line ####
    node SVGenerator.js sv.tsv --chrom chr11 --name sv_chr11 sample.fasta  > chr11_sv.fasta

    # --chrom or -c : chromosome name in FASTA
    # --name or -n  : new name of chromosome with sv
    # argv 1        : TSV file (SV information in it)
    # argv 2        : FASTA file to make sv from


    # TSV format to specify SV
    #
    # TYPE(DEL|INS|INV)  START_BASE_POSIION  LENGTH
    DEL 10354901  341 
    INS 746913  2912
    INV 90110234  567
    DEL 431931  333 
