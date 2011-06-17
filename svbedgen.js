const fs          = require('fs');
const ArgParser   = require('argparser');
const nrand       = require('./lib/normal_random');
const XORShift    = require('./lib/xorshift');
const random      = new XORShift(new Date().getTime(), true); // function
const SVConst     = require('./SVConst');
const SVGenerator = require('./SVGenerator');
const FASTAReader = require('./lib/FASTAReader/FASTAReader');

const numberize = function(v, _default) {
  return (v === false || v === null || isNaN(Number(v))) ? _default : Number(v);
}

const randomInt = function(max) {
  return Math.floor(random() * max);
}


function main() {
  const p = new ArgParser().addOptions([]).addValueOptions(['snprate', 'sv', 'rnames', 'svlen', 'svdev', 'json', 'exename']).parse();

  function showUsage() {
    const cmd = p.getOptions('exename') || (process.argv[0] + ' ' + require('path').basename(process.argv[1]));
    console.error('[synopsis]');
    console.error('\t' + cmd + ' <fasta file>');
    console.error('[options]');
    console.error('\t--snprate\trate to insert SNP default:10000');
    console.error('\t--sv\tthe number of SVs default:20000');
    console.error('\t--rnames\tdesignate rnames to use (comma separated). default null(using ALL rnames)');
    console.error('\t--svlen\tmean length of SVs. default:1500');
    console.error('\t--svdev\tstddev of SVs. default:300');
    console.error('\t' + '--json <json file>\t fasta summary file to shortcut calculation.');
  }


  /* get arguments, options */
  const fastafile = p.getArgs(0);
  if (!fastafile) {
    showUsage();
    process.exit();
  }

  if (! require('path').existsSync(fastafile)) {
    console.error(fastafile + ': No such file.');
    process.exit();
  }


  const snprate = numberize(p.getOptions('snprate'),10000);
  const svnum   = numberize(p.getOptions('sv'),20000);
  const svlen   = numberize(p.getOptions('svlen'), 1500);
  const svdev   = numberize(p.getOptions('svdev'), 300);
  console.error('calculating fasta');
  const json = (function() {
    var ret = p.getOptions('json');
    if (ret) {
      ret = JSON.parse(fs.readFileSync(ret).toString());
    }
    return ret;
  })();

  const fastas = new FASTAReader(fastafile, json);
  const rnames  = (p.getOptions('rnames')) ? p.getOptions('rnames').split(',') : Object.keys(fastas.result);


  console.error('calculating total bases');
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



  console.error('generating SV registration data');
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
        output(type, start, rname, len);
        k++;
      }
      counter++;
    }
  });

  console.error('generating SNP registration data');
  /* generate SNP registration data */
  if (snprate == 0) {
    return;
  }
  const snpnum = Math.floor(total / snprate);
  var snpcount = 0;

  rnames.forEach(function(rname, i) {
    var snps = {};
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
        output(type, pos, rname, to);
        j++;
      }
      counter++;
    }
  });
}

const output = function(type, pos, rname, op) {
  var end = (type == 'SNP') 
    ? Number(pos) + 1
    : Number(pos) + Number(op);

  //console.log(Array.prototype.join.call(arguments, '\t'));
  console.log([rname, pos, end, type, op].join('\t'));
}

if (process.argv[1] === __filename) {
  main();
}
