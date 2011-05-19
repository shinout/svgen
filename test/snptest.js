const SVGenerator = require('../SVGenerator');
const SVConst     = require('../SVConst');
const test        = require('./shinout.test');
const fs          = require('fs');
const dna         = require('../lib/dna');
const EventEmitter = require('events').EventEmitter;

var snpgen = new SVGenerator({
  path    : __dirname + '/long.fa',
  chrom   : 'test_reference',
  svchrom : 'sv_test'
});

snpgen.registerSVFromTSVFile(__dirname + '/withSNP.tsv')

test('equal', snpgen.snps.length, 4, 'invalid count: SNPs');
test('result', 'SNP register test');

/* SNP genotyping test */

var snpgen = new SVGenerator({
  path    : __dirname + '/statement.fa',
  chrom   : 'statement',
  svchrom : 'sv_test'
});

snpgen.registerSNP(4, 1);
snpgen.registerSNP(8, 3);
snpgen.registerSNP(14, 3);
snpgen.registerSNP(93, 1);
snpgen.registerSNP(103, 4);
snpgen.registerSNP(143, 5);
snpgen.registerSNP(163, 2);
snpgen.registerSNP(169, 3);
snpgen.registerSNP(199, 1);
snpgen.registerSNP(223, 5);
snpgen.registerSNP(239, 4);
snpgen.registerSNP(259, 2);
snpgen.registerSNP(275, 1);
snpgen.registerSNP(290, 5);
snpgen.registerSNP(309, 5);


var stream = {
  d : '',
  write: function(d) {
    this.d += d;
  },
  end: function() {
    console.log(this.d);
  },
  on: function(){
  }
}
snpgen.genotype(stream);

