const fs          = require('fs');
const ArgParser   = require('argparser');
const nrand       = require('./lib/normal_random');
const XORShift    = require('./lib/xorshift');
const random      = new XORShift(new Date().getTime(), true); // function
const SVConst     = require('./SVConst');
const SVGenerator = require('./SVGenerator');
const FASTAReader = require('./lib/FASTAReader/FASTAReader');
const dna         = require('./lib/dna');
const WSelection  = require('./lib/WeightedSelection');

const numberize = function(v, _default) {
  return (v === false || v === null || isNaN(Number(v))) ? _default : Number(v);
}

const randomInt = function(max) {
  return Math.floor(random() * max);
}


function main() {
  const p = new ArgParser().addOptions([]).addValueOptions([
    'snprate', 'sv', 'rnames', 'svlen', 'svdev', 'json', 'exename',
    'inslen', 'insdev', 'dellen', 'deldev', 'invlen', 'invdev', 'duplen', 'dupdev'
  ]).parse();

  function showUsage() {
    const cmd = p.getOptions('exename') || (process.argv[0] + ' ' + require('path').basename(process.argv[1]));
    console.error('[synopsis]');
    console.error('\t' + cmd + ' <fasta file>');
    console.error('[options]');
    console.error('\t--snprate\treciprocal rate to insert SNP default:10000 (1/10000). if 0, then no SNPs are registered.');
    console.error('\t--sv <int>\tthe number of SVs to register. default:10000');
    console.error('\t--rnames <string>\tdesignate rnames to use (comma separated). default null(using ALL rnames)');
    console.error('\t--svlen <int>\tmean length of SVs. default:1500');
    console.error('\t--svdev <int>\tstddev of SVs. default:300');
    console.error('\t--inslen <int>\tmean length of insertions. default: the same as "svlen" option. if 0, then no insertions are registered.');
    console.error('\t--insdev <int>\tstddev of insertions. default: the same as "svlen" option');
    console.error('\t--dellen <int>\tmean length of deletions. default: the same as "svlen" option. if 0, then no deletions are registered.');
    console.error('\t--deldev <int>\tstddev of deletions. default: the same as "svlen" option');
    console.error('\t--invlen <int>\tmean length of inversions. default: the same as "svlen" option. if 0, then no inversions are registered.');
    console.error('\t--invdev <int>\tstddev of inversions. default: the same as "svlen" option');
    console.error('\t--duplen <int>\tmean length of tandem duplications. default: 50. if 0, then no tandem duplications are registered.');
    console.error('\t--dupdev <int>\tstddev of tandem duplications. default: 20');
    console.error('\t--tralen <int>\tmean length of translocations. default: the same as "svlen" option. if 0, then no translocations are registered.');
    console.error('\t--tradev <int>\tstddev of translocations. default: the same as "svlen" option');
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
  const svnum   = numberize(p.getOptions('sv'),10000);
  const svlen   = numberize(p.getOptions('svlen'), 1500);
  const svdev   = numberize(p.getOptions('svdev'), 300);

  const lens = {
    INS: numberize(p.getOptions('inslen'), svlen),
    DEL: numberize(p.getOptions('dellen'), svlen),
    INV: numberize(p.getOptions('invlen'), svlen),
    DUP: numberize(p.getOptions('duplen'), 50),
    TRA: numberize(p.getOptions('tralen'), svlen)
  };

  const devs = {
    INS: numberize(p.getOptions('insdev'), svdev),
    DEL: numberize(p.getOptions('deldev'), svdev),
    INV: numberize(p.getOptions('invdev'), svdev),
    DUP: numberize(p.getOptions('dupdev'), 20),
    TRA: numberize(p.getOptions('tradev'), svdev)
  };

  console.error('calculating fasta');
  const json = (function() {
    var ret = p.getOptions('json');
    if (ret) {
      ret = JSON.parse(fs.readFileSync(ret).toString());
    }
    return ret;
  })();

  const fastas  = new FASTAReader(fastafile, json);
  const rnames  = (p.getOptions('rnames')) ? p.getOptions('rnames').split(',') : Object.keys(fastas.result);

  const tselect = new WSelection({
    INS : (lens.INS > 0) ? 1 : 0, 
    DEL : (lens.DEL > 0) ? 1 : 0, 
    INV : (lens.INV > 0) ? 1 : 0, 
    DUP : (lens.DUP > 0) ? 1 : 0, 
    TRA : (lens.TRA > 0) ? 1 : 0
  }, random);

  const rselect = new WSelection((function(){
    var ret = {};
    rnames.forEach(function(rname) {
      ret[rname] = fastas.result[rname].getEndPos();
    });
    return ret;
  })() , random); 


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
      var type  = tselect.random();
      var len   = Math.floor(nrand(lens[type], devs[type], random) + 0.5);
      var start = randomInt(endpos - len -1);
      var extra = (function() {
        if ( type == 'INS') {
          return dna.getRandomFragment(len);
        }
        else if ( type == 'TRA') {
          // get fragment from a random position.
          var rn, fa, st;
          do {
            rn = rselect.random();
            fa = fastas.result[rn];
            st = randomInt(fa.getEndPos() - len - 1);
          } while (fastas.hasN(rn, st, len)); // TODO escape loop
          return rn + ':' + st + ':' + len;
        }
        else if ( type == 'DUP') {
          return SVGenerator.getTandemDuplicationRepeatNumber();
        }
        else {
          return '*';
        }
      })();

      if (svgen.registerSV(SVConst[type], start, len)) {
        output(type, start, rname, len, extra);
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
        output(type, pos, rname, 1, to);
        j++;
      }
      counter++;
    }
  });
}

const output = function(type, pos, rname, len, extra) {
  var end = (type == 'SNP') 
    ? Number(pos) + 1
    : Number(pos) + Number(len);

  //console.log(Array.prototype.join.call(arguments, '\t'));
  console.log([rname, pos, end, type, len, extra].join('\t'));
}

if (process.argv[1] === __filename) { main(); }
