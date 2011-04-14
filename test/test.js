var SVGenerator = require('../SVGenerator');
var SVConst = require('../SVConst');
var test = require('./shinout.test');
var fs = require('fs');

var errs = [];
SVGenerator.error = function(v) {
  console.log(v);
  errs.push(v);
}

/* setting up test */
var svgen = new SVGenerator();
var e = errs.shift();
test('ok', e.match(/such file/), 'No error occurred when no file designated.');
test('equal', svgen.chrom, 'chr1', 'default chrom name error');
test('result', 'setting up test');


/* random flagment test */
var dna = SVGenerator.getRandomFlagment(1000);
var T = dna.length - dna.split('T').join('').length;
var C = dna.length - dna.split('C').join('').length;
var A = dna.length - dna.split('A').join('').length;
var G = dna.length - dna.split('G').join('').length;

console.log( Math.max(T,C,A,G) / Math.min(T,C,A,G));
test('equal', dna.length, 1000, 'invalid length of random flagment'); 
test('equal', T + C + A + G, 1000, 'invalid character of random flagment'); 
test('ok', Math.max(T,C,A,G) / Math.min(T,C,A,G) < 1.5, 'invalid deviation of random flagment'); 
test('result', 'random flagment test');

/* pos2index test */
var svgen = new SVGenerator({
  //path : __dirname + '/../../../data/references/chr11.fa',
  path : __dirname + '/long.fa',
  chrom: 'chr11'
});

function getChar(pos, len) {
  var idx = svgen.pos2index(pos);
  var fd = fs.openSync(svgen.path, 'r');
  var ret = fs.readSync(fd, len, idx)[0];
  fs.closeSync(fd);
  return ret;
}
test('equal', getChar(3,8), '34567890', 'unexpected position value');
test('equal', getChar(103,8), '34567890', 'unexpected position value');
test('equal', getChar(105143,1), 'h', 'unexpected position value');
test('equal', getChar(3000*50,1), 'o', 'unexpected position value');
test('result', 'pos2index test');


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

svgen.registerInv(130130877, 404);
test('equal', svgen.svs.length, 5, 'Inversion event didn\'t registered.');

svgen.registerDel(130130900, 356);
var e = errs.shift();
test('ok', e && e.match(/duplicated/), 'No error occurred when duplicated region is given.');

svgen.registerDel(130130522, 355);
var e = errs.shift();
test('ok', e && e.match(/duplicated/), 'No error occurred when duplicated region is given.');

svgen.registerDel(130130521, 355);
test('equal', svgen.svs.length, 6, 'Deletion event didn\'t registered.');

test('result', 'registerSV test');


/* GenomeStream test */

var mock = {
  result: '',
  write: function(v) {
    this.result+=v;
  },
  emit: function() {
  }
}

var gs = new SVGenerator.GenomeStream({prelen: 7, nextream: mock});
gs.write('>chr11\n'+(function(){var ret = "";var i = 0;while(i <1000) {ret +="a"+i+"\n"; i++;} return ret;})());
gs.write('>chr11\n'+(function(){var ret = "";var i = 0;while(i <1000) {ret +="a"+i; i++;} return ret;})());
gs.end();
test('equal', mock.result.slice(0, 7), '>chr11\n', 'GenomeStream.write : invalid output .');
test('ok', !mock.result.substr(7, 50).match(/\n/), 'GenomeStream.write : invalid output .');
test('ok', mock.result.substr(7, 51).match(/\n/), 'GenomeStream.write : invalid output .');
test('result', 'GenomeStream test');
mock.result = '';


/* SVStream test */
/* unit test */
var result = SVGenerator.SVStream.complStrand('atgcATGC\nNnAATT');

test('equal', result, 'tacgTACG\nNnTTAA', 'SVStream.complStrand: invalid output .');

var result = SVGenerator.SVStream.makeDeletion({start: 3, end: 9}, 'atgcATGCAATTGGCC', 0);
test('equal', result, 'atgATTGGCC', 'SVStream.makeDeletion: invalid output .');

var result = SVGenerator.SVStream.makeInsertion({start: 2, end: 6, flagment: 'TTTT'}, 'atgcATGCAATTGGCC', 0);
test('equal', result, 'atTTTTgcATGCAATTGGCC', 'SVStream.makeInsertion: invalid output .');

var result = SVGenerator.SVStream.makeInversion({start: 0, end: 4}, 'atgcATGCAATTGGCC', 0);
test('equal', result, 'gcatATGCAATTGGCC', 'SVStream.makeInsertion: invalid output .');


/* more complicated offset and length */
var result = SVGenerator.SVStream.makeInversion(
  {start: 2, end: 6}, 'AACcCcTTTTGGG', 0
);
test('equal', result, 'AAgGgGTTTTGGG', 'SVStream.makeInsertion: invalid output .');
test('result', 'SVStream static function test');

var result = SVGenerator.SVStream.makeInversion(
  {start: 4002, end: 4006}, 'AACcCcTTTTGGG', 4000
);
test('equal', result, 'AAgGgGTTTTGGG', 'SVStream.makeInsertion: invalid output .');
test('result', 'SVStream static function test');


var svstream = new SVGenerator.SVStream({nextream: mock, svs: [{
  type: SVConst.DEL,
  start: 6 + 3,
  end: 6 + 3 + 5 -1 + 1
}]});
svstream.write('>chr11\nAATGGCCAATGGCCAAAAATTTT');
//                        ^^^^^
svstream.end();
test('equal', mock.result, '>chr11\nAAAATGGCCAAAAATTTT', 'SVStream deletion failed.');
mock.result = '';


var svstream = new SVGenerator.SVStream({nextream: mock, svs: [{
  type: SVConst.INS,
  start: 6 + 4,
  end: 6 + 4 + 5 -1 + 1,
  flagment: 'ctgac'
}]});
svstream.write('>chr11\nAATGGCCAATGGCCAAAAATTTT');
//                         ^
svstream.end();
test('equal', mock.result, '>chr11\nAATctgacGGCCAATGGCCAAAAATTTT');
mock.result = '';


var svstream = new SVGenerator.SVStream({nextream: mock, svs: [{
  type: SVConst.INV,
  start: 6 + 10,
  end: 6 + 10 + 6 -1 + 1
}]});
svstream.write('>chr11\nAATGGCCAATGGCCAAAAATTTT');
//                               ^^^^^^
svstream.end();
test('equal', mock.result, '>chr11\nAATGGCCAATGGCCAAAAATTTT');
//                                           ^^^^^^
mock.result = '';

test('result', 'SVStream test');


/* genotyping */
require('./genotype');
