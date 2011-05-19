const fs          = require('fs');
const ArgParser   = require('argparser');
const nrand       = require('./lib/normal_random');
const XORShift    = require('./lib/xorshift');
const random      = new XORShift(new Date().getTime(), true); // function
const SVConst     = require('./SVConst');
const SVGenerator = require('./SVGenerator');
const FASTAReader = require('./lib/FASTAReader/FASTAReader');

const numberize = function(v) {
  return (isNaN(Number(v))) ? null : Number(v);
}

const randomInt = function(max) {
  return Math.floor(random() * max);
}


function stderr(v) {
  process.stderr.write(v + '\n');
}


function main() {
  function showUsage() {
    const cmd = 'node ' + require('path').basename(process.argv[1]);
    stderr('[synopsis]');
    stderr(cmd + ' [--snprate=10000] [--sv=20000] [--rnames=false] [--svlen=1500] [--svdev=300] <fasta file>');
  }

  const p = new ArgParser().addOptions([]).addValueOptions(['snprate', 'sv', 'rnames', 'svlen', 'svdev']).parse();

  /* get arguments, options */
  const fastafile = p.getArgs(0);
  if (!fastafile) {
    showUsage();
    process.exit();
  }

  if (! require('path').existsSync(fastafile)) {
    gen.error(fastafile + ': No such file.');
    process.exit();
  }


  const snprate = numberize(p.getOptions('snprate')) || 10000;
  const svnum   = numberize(p.getOptions('sv')) || 20000;
  const svlen   = numberize(p.getOptions('svlen')) || 1500;
  const svdev   = numberize(p.getOptions('svdev')) || 300;
  stderr('calculating fasta');
  const fastas = new FASTAReader(fastafile);
  const rnames  = (p.getOptions('rnames')) ? p.getOptions('rnames').split(',') : Object.keys(fastas.result);


  stderr('calculating total bases');
  /* calculate total bases */
  const ends  = {};
  const total = (function() {
    var t = 0;
    rnames.forEach(function(rname) {
      const fasta = fastas.result[rname];
      var end = fasta.getEndPos();
      t += end;
      ends[rname] = end;
    });
    return t;
  })();



  stderr('generating SV registration data');
  /* generate SV registration data */
  var svcount = 0;   // the number of svs that have already added
  rnames.forEach(function(rname, i) {
    const fasta        = fastas.result[rname];
    const svgen        = new SVGenerator({freader: fastas, chrom: rname});
    const endpos       = fasta.getEndPos();
    const localSVCount = (rnames.length == i+1)
      ? svnum - svcount 
      : parseInt(svnum * endpos / total);
    const LIMIT = localSVCount * 30;

    svcount += localSVCount;

    var k = 0, counter = 0;
    while( k < localSVCount && counter < LIMIT) {
      var type  = SVConst.types[randomInt(3)];
      var len   = Math.floor(nrand(svlen, svdev, random) + 0.5);
      var start = randomInt(endpos - len -1);

      if (svgen.registerSV(SVConst[type], start, len)) {
        output(type, start, len, rname);
        k++;
      }
      counter++;
    }
  });

  stderr('generating SNP registration data');
  /* generate SNP registration data */
  var snps = {};
  var snpnum = Math.floor(total / snprate);
  const snpcount = 0;

  rnames.forEach(function(rname, i) {
    const fasta        = fastas.result[rname];
    const svgen        = new SVGenerator({freader: fastas, chrom: rname});
    const endpos       = fasta.getEndPos();
    const localSNPCount = (rnames.length == i+1)
      ? snpnum - snpcount 
      : parseInt(snpnum * endpos / total);

    const LIMIT = localSNPCount * 30;
    snpcount += localSNPCount;

    var j = 0, counter = 0;
    while( j < localSNPCount && counter < LIMIT) {
      var type  = 'SNP';
      var to    = (function() {
        var v = randomInt(1001);
        if      (v < 333)  return 1;
        else if (v < 666)  return 2;
        else if (v < 999)  return 3;
        else if (v == 999) return 4;
        else               return 5;
      })();
      var pos   = randomInt(ends[rname] -1) + 1;

      /* convert position, rname */

      if (!snps[pos]) {
        snps[pos] = true;
        output(type, pos, to, rname);
        j++;
      }
      counter++;
    }
  });
}

const output = function() {
  console.log(Array.prototype.join.call(arguments, '\t'));
  console.log.apply(this, arguments);
}

if (process.argv[1] === __filename) {
  main();
}


