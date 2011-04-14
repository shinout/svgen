var SVGenerator = require('../SVGenerator');
var SVConst = require('../SVConst');
var test = require('./shinout.test');

/* genotyping */

var svgen = new SVGenerator({
  path : __dirname + '/short.fa',
  chrom: 'chr11',
  bufferSize: 300
});


svgen.registerDel(1, 10);
svgen.registerInv(51, 30);
svgen.registerIns(101, 50);
svgen.registerDel(501, 999);
svgen.genotype();
