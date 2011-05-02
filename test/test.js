var SVGenerator = require('../SVGenerator');
var SVConst     = require('../SVConst');
var test        = require('./shinout.test');
var fs          = require('fs');
var dna         = require('../lib/dna');
var EventEmitter = require('events').EventEmitter;

var errs = [];
SVGenerator.error = function(v) {
  console.log(v);
  errs.push(v);
}

/* lib/dna test */
var svgen = new SVGenerator();
var e = errs.shift();
test('ok', e.match(/such file/), 'No error occurred when no file designated.');
test('result', 'setting up test');


/* random flagment test */
var DNA = dna.getRandomFlagment(1000);
var T   = DNA.length - DNA.split('T').join('').length;
var C   = DNA.length - DNA.split('C').join('').length;
var A   = DNA.length - DNA.split('A').join('').length;
var G   = DNA.length - DNA.split('G').join('').length;

//console.log( Math.max(T,C,A,G) / Math.min(T,C,A,G));
test('equal', DNA.length, 1000, 'invalid length of random flagment'); 
test('equal', T + C + A + G, 1000, 'invalid character of random flagment'); 
test('ok', Math.max(T,C,A,G) / Math.min(T,C,A,G) < 1.5, 'invalid deviation of random flagment'); 
test('result', 'random flagment test');

var result = dna.complStrand('atgcATGC\nNnAATT');
test('equal', result, 'tacgTACG\nNnTTAA', 'dna.complStrand: invalid output .');

var svgen = new SVGenerator({
  path : __dirname + '/long.fa'
});

test('equal', svgen.chrom, 'chr11', 'refname auto detection error');

/* registerSV test */
svgen.registerSV('DEL', 100000, 200);
var e = errs.shift();
test('ok', e && e.match(/Number/), 'No error occurred when string is given to type.');

svgen.registerDel(-3, 200);
var e = errs.shift();
test('ok', e && e.match(/start/), 'No error occurred when minus value is given to start.');

svgen.registerDel(0, 200);
var e = errs.shift();
test('ok', e && e.match(/start/), 'No error occurred when zero value is given to start.');

svgen.registerIns(300, -3);
var e = errs.shift();
test('ok', e && e.match(/length/), 'No error occurred when minus value is given to length.');

svgen.registerDel(300, 0);
var e = errs.shift();
test('ok', e && e.match(/length/), 'No error occurred when zero is given to length.');


svgen.registerDel(1, 200);
test('equal', svgen.svs.length, 1, 'Deletion event didn\'t registered.');
console.log(svgen.svs[0]);
test('equal', svgen.svs[0].start , 7, 'invalid coordinate.');
test('equal', svgen.svs[0].end , 200 + 7 - 1 + 3 + 1, 'invalid coordinate.');

svgen.registerInv(100, 90);
var e = errs.shift();
test('ok', e && e.match(/duplicated/), 'No error occurred when duplicated region is given.');

svgen.registerDel(199, 90);
var e = errs.shift();
test('ok', e && e.match(/duplicated/), 'No error occurred when duplicated region is given.');

svgen.registerIns(200, 90);
var e = errs.shift();
test('ok', e && e.match(/duplicated/), 'No error occurred when duplicated region is given.');

svgen.registerIns(11111200, 290);
var e = errs.shift();
test('ok', e && e.match(/range/), 'No error occurred when out of range is given.');

svgen.registerIns(201, 90);
test('equal', svgen.svs.length, 2, 'Insertion event didn\'t registered.');
test('equal', svgen.svs[1].flagment.length, 90, 'invalid flagment length in registerIns');

svgen.registerIns(1001, 90, 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA');
test('equal', svgen.svs.length, 3, 'Insertion event didn\'t registered.');
test('equal', svgen.svs[2].flagment.length, 90, 'invalid flagment length in registerIns');
test('equal', svgen.svs[2].flagment.slice(0,20), 'AAAAAAAAAAAAAAAAAAAA', 'invalid flagment in registerIns');

svgen.registerIns(10001, 10, 'TTTTTTTTTTTTTTTTTTTTTTTTTTTT');
test('equal', svgen.svs.length, 4, 'Insertion event didn\'t registered.');
test('equal', svgen.svs[3].flagment.length, 10, 'invalid flagment length in registerIns');
test('equal', svgen.svs[3].flagment, 'TTTTTTTTTT', 'invalid flagment in registerIns');

svgen.registerInv(20877, 404);
test('equal', svgen.svs.length, 5, 'Inversion event didn\'t registered.');

svgen.registerDel(20900, 356);
var e = errs.shift();
test('ok', e && e.match(/duplicated/), 'No error occurred when duplicated region is given.');

svgen.registerDel(20522, 355);
var e = errs.shift();
test('ok', e && e.match(/duplicated/), 'No error occurred when duplicated region is given.');

svgen.registerDel(20521, 355);
test('equal', svgen.svs.length, 6, 'Deletion event didn\'t registered.');


// sort check

var svgen = new SVGenerator({
  path    : __dirname + '/long.fa',
  chrom   : 'test_reference',
  svchrom : 'sv_test'
});

svgen.registerDel(1500, 300);
svgen.registerIns(500, 100);
svgen.registerInv(100, 300);
svgen.registerIns(12200, 1300);
svgen.registerIns(2200, 3300);
svgen.registerDel(2000, 300); // duplicated
svgen.registerInv(1, 20);

test('equal', svgen.svs[0].pos, 1, 'invalid order');
test('equal', svgen.svs[1].pos, 100, 'invalid order');
test('equal', svgen.svs[2].pos, 500, 'invalid order');
test('equal', svgen.svs[3].pos, 1500, 'invalid order');
test('equal', svgen.svs[4].pos, 2200, 'invalid order');
test('equal', svgen.svs[5].pos, 12200, 'invalid order');




/* coordinate test */
var svgen = new SVGenerator({
  path    : __dirname + '/long.fa',
  chrom   : 'test_reference',
  svchrom : 'sv_test'
});


test('equal', svgen.prelen, 'test_reference'.length + 2, 'invalid prelen.');
// lines * linelen + idlen + blank + idlen2
test('equal', svgen.startIdx, 3000 * 51 + 7 + 2 + svgen.prelen, 'invalid startIdx.');
test('equal', 'A', fs.readFileSync(__dirname + '/long.fa').toString().substr(svgen.endIdx-2,1), 'invalid endIdx');
svgen.registerDel(1, 20);
test('equal', svgen.svs[0].start, svgen.startIdx, 'invalid index of registered deletion.');
test('equal', svgen.svs[0].end, svgen.startIdx + 20, 'invalid index of registered deletion.');

svgen.registerIns(51, 20);
test('equal', svgen.svs[1].start, svgen.startIdx + 51, 'invalid start of registered insertion.');
test('equal', svgen.svs[1].end, svgen.startIdx + 51 + 20, 'invalid end of registered insertion.');

svgen.registerIns(100, 52);
test('equal', svgen.svs[2].start, svgen.startIdx + 100, 'invalid start of registered insertion.');
test('equal', svgen.svs[2].end, svgen.startIdx + 100 + 52 + 2, 'invalid end of registered insertion.');

svgen.registerInv(222, 50);
test('equal', svgen.svs[3].start, svgen.startIdx + 221 +4, 'invalid start of registered inversion.');
test('equal', svgen.svs[3].end, svgen.startIdx + 221 + 4 + 51, 'invalid end of registered inversion.');

test('result', 'registerSV test');

/* register from TSV file test */
var svgen = new SVGenerator({
  path    : __dirname + '/long.fa',
  chrom   : 'test_reference',
  svchrom : 'sv_test'
});
var bool = svgen.registerSVFromTSVFile(__dirname + '/error.tsv')
test('ok', !bool, 'error.tsv registered incorrectly.');

var bool = svgen.registerSVFromTSVFile(__dirname + '/one_ok.tsv')
test('equal', svgen.svs.length, 1, 'one_ok.tsv registered incorrectly.');
test('ok', bool, 'one_ok.tsv registered incorrectly.');

var bool = svgen.registerSVFromTSVFile(__dirname + '/sample.tsv')
console.log(svgen.svs);
test('equal', svgen.svs.length, 20+1 -1 -1 - 11, 'sample.tsv did not registered correctly.');
// one duplication with one_ok and one dupe with itself and 11 out of range
test('ok', bool, 'sample.tsv did not registered correctly.');
test('equal', svgen.svs[0].type, SVConst.DEL, 'sample.tsv did not registered correctly.');
test('equal', svgen.svs[0].pos, 3385, 'sample.tsv did not registered correctly.');

test('equal', svgen.svs[1].type, SVConst.INS, 'sample.tsv did not registered correctly.');
test('equal', svgen.svs[1].pos, 4912, 'sample.tsv did not registered correctly.');

test('equal', svgen.svs[2].type, SVConst.INV, 'sample.tsv did not registered correctly.');
test('equal', svgen.svs[2].pos, 5912, 'sample.tsv did not registered correctly.');

test('result', 'register sv from tsv test');
process.exit();




/* SVConst test */

var result = SVConst.makeDeletion({start: 3, end: 9}, 'atgcATGCAATTGGCC', 0);
test('equal', result, 'atgATTGGCC', 'SVConst.makeDeletion: invalid output .');

var result = SVConst.makeInsertion({start: 2, end: 6, flagment: 'TTTT'}, 'atgcATGCAATTGGCC', 0);
test('equal', result, 'atTTTTgcATGCAATTGGCC', 'SVConst.makeInsertion: invalid output .');

var result = SVConst.makeInversion({start: 0, end: 4}, 'atgcATGCAATTGGCC', 0);
test('equal', result, 'gcatATGCAATTGGCC', 'SVConst.makeInsertion: invalid output .');


/* more complicated offset and length */
var result = SVConst.makeInversion(
  {start: 2, end: 6}, 'AACcCcTTTTGGG', 0
);
test('equal', result, 'AAgGgGTTTTGGG', 'SVConst.makeInsertion: invalid output .');
test('result', 'SVConst static function test');

var result = SVConst.makeInversion(
  {start: 4002, end: 4006}, 'AACcCcTTTTGGG', 4000
);
test('equal', result, 'AAgGgGTTTTGGG', 'SVConst.makeInsertion: invalid output .');
test('result', 'SVConst static function test');


/* SVStream test */
var svstream = new SVGenerator.SVStream({svs: [{
  type: SVConst.DEL,
  start: 3,
  end: 6
}]});

svstream.on('data', function(data) {
  if (!this.result) this.result = '';
  if (!this.count) this.count = 1;
  this.result += data;
  //console.log((this.count++) + ':'+data);
});

svstream.on('end', function() {
  //console.log(this.result);
  test('equal', this.result, 'AATCAATGGCCAAAAATTTT', 'SVStream deletion failed.');
  test('result', 'SVStream test1');
});

svstream.write('AATGGCCAATGGCCAAAAATTTT');
//                 ^^^
svstream.end();





var svstream2 = new SVGenerator.SVStream({svs: [{
  type: SVConst.INS,
  start: 3,
  end: 3 + 5,
  flagment: 'ctgac'
}]});

svstream2.on('data', function(data) {
  if (!this.result) this.result = '';
  if (!this.count) this.count = 1;
  this.result += data;
});

svstream2.on('end', function() {
  test('equal', this.result, 'AATctgacGGCCAATGGCCAAAAATTTT');
  test('result', 'SVStream test2');
});

svstream2.write('AATGGCCAATGGCCAAAAATTTT');
//                  ^
svstream2.end();

var svstream3 = new SVGenerator.SVStream({svs: [{
  type: SVConst.INV,
  start: 10,
  end: 10 + 6
}]});


svstream3.on('data', function(data) {
  if (!this.result) this.result = '';
  if (!this.count) this.count = 1;
  this.result += data;
  //console.log((this.count++) + ':'+data);
});

svstream3.on('end', function() {
  //console.log(this.result);
  test('equal', this.result, 'AATGGCCAATTTGGCCAAATTTT');
  test('result', 'SVStream test3');
});

svstream3.write('AATGGCCAATGGCCAAAAATTTT');
//                         ^^^^^^
svstream3.end();


/* genotyping */
var svgen = new SVGenerator({
  path    : __dirname + '/long.fa',
  chrom   : 'test_reference',
  svchrom : 'sv_test'
});

svgen.registerDel(1500, 300);
svgen.registerIns(500, 100);
svgen.registerInv(100, 300);
svgen.registerDel(12200, 1300);
svgen.registerIns(2200, 3300);
svgen.registerInv(1, 20);

var vals = svgen.genotype(null, true);
var svs  = vals.svs;
test('equal', svs.length, 6, 'invalid svs length');
test('equal', svs[0].start, 0, 'invalid index of registered sv.');
test('equal', svs[0].end, 20, 'invalid index of registered sv.');
test('equal', svs[1].start, 99 + 1, 'invalid index of registered sv.');
test('equal', svs[1].end, 99 + 1 + 300+6, 'invalid index of registered sv.');

test('equal', svgen.fasta.fetch(1,10), 'CGCGGCTAGC', 'FASTAReader check failed.');

var options = vals.options;
options.end = options.start + 9;
options.bufferSize = 100;
var stream = fs.createReadStream(svgen.path, options);

stream.on('data', function(data) {
  test('equal', data.toString(), 'CGCGGCTAGC', 'FASTAReader check failed.');
});

stream.on('end', function() {
  var wstream = new EventEmitter();
  wstream.result = '';

  wstream.write = function(data) {
    this.result += data;
  }

  wstream.end = function() {

    test('equal', this.result.slice(0,'>sv_test'.length), '>sv_test', 'result: invalid chromname');
    test('equal', this.result.charAt('>sv_test'.length + 1 + 50), '\n', 'result: invalid linelen');

    test('equal', this.result.slice(0,'>sv_test'.length), '>sv_test', 'result: invalid chromname');
    test('equal', this.result.charAt('>sv_test'.length + 1 + 50), '\n', 'result: invalid linelen');


    var genome = this.result.slice('>sv_test'.length).split('\n').join('');


    /** INVERSION No.1 **/
    test('equal',genome.substr(0, 20), dna.complStrand(svgen.fasta.fetch(1, 20)).split('').reverse().join('') , 'inversion check failed');

    // modify genome for inversion 
    genome = svgen.fasta.fetch(1,20) + genome.slice(20);
    test('equal', genome.substr(0, 30), svgen.fasta.fetch(1, 30), 'inversion modification failed');


    /** INVERSION No.2 **/
    test('equal',genome.substr(99, 300), dna.complStrand(svgen.fasta.fetch(100, 300)).split('').reverse().join('') , 'inversion check failed');

    // modify genome for inversion 
    genome = genome.slice(0, 99) + svgen.fasta.fetch(100, 300) + genome.slice(99+300);
    test('equal', genome.substr(0, 400), svgen.fasta.fetch(1, 400), 'inversion modification failed');


    /** INSERTION No.1 **/
    test('equal', genome.substr(399, 100) + genome.substr(599, 100), svgen.fasta.fetch(400, 200), 'insertion check failed');

    // modify genome for insertion
    genome = genome.slice(0, 499) + genome.slice(599);
    test('equal', genome.substr(399, 200), svgen.fasta.fetch(400, 200), 'insertion modification failed');


    /** DELETION No.1 **/
    test('equal', genome.substr(1400-1, 200), svgen.fasta.fetch(1400, 100) + svgen.fasta.fetch(1800, 100), 'deletion check failed');

    // modify genome for deletion 
    genome = genome.slice(0, 1499) + svgen.fasta.fetch(1500, 300) + genome.slice(1499);
    test('equal', genome.substr(1400-1, 500), svgen.fasta.fetch(1400, 500), 'deletion modification failed');


    /** INSERTION No.2 **/
    test('equal', genome.substr(2200-100-1, 100) + genome.substr(2200+3300-1, 100), svgen.fasta.fetch(2200-100, 200), 'insertion check failed');

    // modify genome for insertion
    genome = genome.slice(0, 2200-1) + genome.slice(2200+3300-1);
    test('equal', genome.substr(2100-1, 200), svgen.fasta.fetch(2100, 200), 'insertion modification failed');

    /** DELETION No.2 **/
    test('equal', genome.substr(12200-100-1, 200), svgen.fasta.fetch(12200-100, 100) + svgen.fasta.fetch(12200+1300, 100), 'deletion check failed');

    // modify genome for deletion 
    genome = genome.slice(0, 12200-1) + svgen.fasta.fetch(12200, 1300) + genome.slice(12200-1);
    test('equal', genome.substr(12200-100-1, 1500), svgen.fasta.fetch(12200-100, 1500), 'deletion modification failed');

    /** RECONSTRUCTION CHECK **/
    test('equal', genome, svgen.fasta.fetch(1, genome.length));


    test('result', 'genotype test');
  }

  svgen.genotype(wstream);
});


