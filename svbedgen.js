const fs          = require('fs');
const ArgParser   = require('argparser');
const nrand       = require('./lib/normal_random');
const XORShift    = require('./lib/xorshift');
const random      = new XORShift(new Date().getTime(), true); // function
const SVConst     = require('./SVConst');
const SVGen       = require('./svgen');
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
    'inslen', 'insdev', 'dellen', 'deldev', 'invlen', 'invdev', 
    'duplen', 'dupdev', 'duprep', 'duprdev',
    'insnum', 'delnum', 'invnum', 'dupnum', 'tranum'
  ]).parse();

  function showUsage() {
    const cmd = p.getOptions('exename') || (process.argv[0] + ' ' + require('path').basename(process.argv[1]));
    console.error('[synopsis]');
    console.error('\t' + cmd + ' <fasta file>');
    console.error('[options]');
    console.error('');

    console.error('\t** length: the lengths of SVs.');
    console.error('\t--svlen <int>\tmean length of SVs. default:1500');
    console.error('\t--inslen <int>\tmean length of insertions. default: the same as "svlen" option. if 0, then no insertions are registered.');
    console.error('\t--dellen <int>\tmean length of deletions. default: the same as "svlen" option. if 0, then no deletions are registered.');
    console.error('\t--invlen <int>\tmean length of inversions. default: the same as "svlen" option. if 0, then no inversions are registered.');
    console.error('\t--duplen <int>\tmean length of tandem duplications. default: the same as "svlen" option. if 0, then no tandem duplications are registered.');
    console.error('\t--tralen <int>\tmean length of translocations. default: the same as "svlen" option. if 0, then no translocations are registered.');
    console.error('');

    console.error('\t** deviation : standard deviation of each SVs.');
    console.error('\t--svdev <int>\tstddev of SVs. default:300');
    console.error('\t--insdev <int>\tstddev of insertions. default: the same as "svlen" option');
    console.error('\t--deldev <int>\tstddev of deletions. default: the same as "svlen" option');
    console.error('\t--invdev <int>\tstddev of inversions. default: the same as "svlen" option');
    console.error('\t--dupdev <int>\tstddev of tandem duplications. default: the same as "svdev" option.');
    console.error('\t--tradev <int>\tstddev of translocations. default: the same as "svlen" option');
    console.error('');

    console.error('\t** number: total number of SVs.');
    console.error('\t--sv <int>\tthe number of SVs to register. default:10000. If each SV is specifically designated its number, that configuration priors to this.');
    console.error('');

    console.error('\t** rate : the rate of each SV type, determined by the following rule.');
    console.error('\t*  let "total" be the sum of each num, i.e. insnum, delnum, invnum, dupnum and tranum.');
    console.error('\t*  the rate of each SV is  x / total, where x is one of insnum, delnum, invnum, dupnum or tranum.');
    console.error('\t--insnum <int>\tfor insertions. default: 200.');
    console.error('\t--delnum <int>\tfor deletions. default: 200.');
    console.error('\t--invnum <int>\tfor inversions. default: 200.');
    console.error('\t--dupnum <int>\tfor tandem duplications. default: 200.');
    console.error('\t--tranum <int>\tfor translocations. default: 200.');
    console.error('');

    console.error('\t** repeat number of tandem duplication.');
    console.error('\t--duprep <int>\tmean repeat folds of tandem duplications. default: 40.');
    console.error('\t--duprdev <int>\tstddev of repeat folds of tandem duplications. default: 10.');
    console.error('');

    console.error('\t** configuration for SNP.');
    console.error('\t--snprate\treciprocal rate to insert SNP default:10000 (1/10000). if 0, then no SNPs are registered.');
    console.error('');

    console.error('\t** configuration for a fasta file.');
    console.error('\t' + '--json <json file>\t fasta summary file to shortcut calculation.');
    console.error('');

  }


  /* get fasta */
  const fastafile = p.getArgs(0);
  if (!fastafile) {
    showUsage();
    process.exit();
  }

  if (! require('path').existsSync(fastafile)) {
    console.error(fastafile + ': No such file.');
    process.exit();
  }


  console.error('collecting FASTA information');
  // config data in JSON format for FASTAReader
  const json = (function() {     var ret = p.getOptions('json');
    if (ret) {
      ret = JSON.parse(fs.readFileSync(ret).toString());
    }
    return ret;
  })();

  const fastas  = new FASTAReader(fastafile, json);
  const rnames  = Object.keys(fastas.result);
  
  /* get options */
  const snprate = numberize(p.getOptions('snprate'),10000);
  const svnum   = numberize(p.getOptions('sv'),10000);
  const svlen   = numberize(p.getOptions('svlen'), 1500);
  const svdev   = numberize(p.getOptions('svdev'), 300);

  const lens = {
    INS: numberize(p.getOptions('inslen'), svlen),
    DEL: numberize(p.getOptions('dellen'), svlen),
    INV: numberize(p.getOptions('invlen'), svlen),
    DUP: numberize(p.getOptions('duplen'), svlen),
    TRA: numberize(p.getOptions('tralen'), svlen)
  };

  const devs = {
    INS: numberize(p.getOptions('insdev'), svdev),
    DEL: numberize(p.getOptions('deldev'), svdev),
    INV: numberize(p.getOptions('invdev'), svdev),
    DUP: numberize(p.getOptions('dupdev'), svdev),
    TRA: numberize(p.getOptions('tradev'), svdev)
  };

  const nums = {
    INS: numberize(p.getOptions('insnum'), 200),
    DEL: numberize(p.getOptions('delnum'), 200),
    INV: numberize(p.getOptions('invnum'), 200),
    DUP: numberize(p.getOptions('dupnum'), 200),
    TRA: numberize(p.getOptions('tranum'), 200)
  };

  const duprep  = numberize(p.getOptions('duprep'), 40);
  const duprdev = numberize(p.getOptions('duprdev'), 10);


  // weighted selection of rnames
  const rselect = new WSelection((function(){
    var ret = {};
    rnames.forEach(function(rname) {
      ret[rname] = fastas.result[rname].getEndPos();
    });
    return ret;
  })() , random); 



  console.error('generating SV registration data');
  // weighted selection of SV types
  const tselect = new WSelection({
    INS : (lens.INS > 0) ? nums.INS : 0, 
    DEL : (lens.DEL > 0) ? nums.DEL : 0, 
    INV : (lens.INV > 0) ? nums.INV : 0, 
    DUP : (lens.DUP > 0) ? nums.DUP : 0, 
    TRA : (lens.TRA > 0) ? nums.TRA : 0
  }, random);

  // frequency of SV events for each rname
  const sv_counts = {};
  var rname_rand;
  for (var i=0; i<svnum; i++) {
    rname_rand = rselect.random();
    if (!sv_counts[rname_rand]) {
      sv_counts[rname_rand] = 0;
    }
    sv_counts[rname_rand]++;
  }

  const svgen = new SVGen(fastas);
  rnames.forEach(function(rname) {
    var cnt = sv_counts[rname] || 0;
    const fasta        = fastas.result[rname];
    const endpos       = fasta.getEndPos();
    const LIMIT = cnt * 30;

    for (var i=0, e=0; i<cnt && e<LIMIT; i++) {
      var type  = tselect.random();
      var len   = Math.floor(nrand(lens[type], devs[type], random) + 0.5);
      if (len < 1) {
        e++;
        i--;
        continue;
      }
      var start = randomInt(endpos - len -1);
      var extra = (function() {
        switch (type) {
        case 'INS':
          return dna.getRandomFragment(len);
        case 'TRA':
          // get fragment from a random position.
          var rn, fa, st, extra_canditate;
          do {
            rn = rselect.random();
            fa = fastas.result[rn];
            st = randomInt(fa.getEndPos() - len - 1);
            extra_canditate = rn + ':' + st + ':' + len;
          } while (!SVGen.valid.TRA(fastas, rname, start, len, extra_canditate)); // TODO escape loop
          return extra_canditate;
        case 'DUP':
          return Math.max(Math.floor(nrand(duprep, duprdev, random)), 2);
        default:
          return '*';
        }
      })();

      try {
        svgen.register(rname, start, len, type, extra);
        output(type, start, rname, len, extra);
      }
      catch (err) {
        //console.error(err.message);
        e++;
        i--;
      }
    }
  });

  console.error('generating SNP registration data');
  /* generate SNP registration data */
  if (snprate == 0) {
    return;
  }

  // weighted selection of SNP types
  const snpselect = new WSelection({
    1 : 1000,
    2 : 1000,
    3 : 1000,
    4 : 3,
    5 : 3 
  });

  // frequency of SNP events for each rname
  const snp_counts  = {};
  const snp_total   = rselect.total() / snprate;
  // var rname_rand; // the same variable exists above.
  for (var i=0; i<snp_total; i++) {
    rname_rand = rselect.random();
    if (!snp_counts[rname_rand]) {
      snp_counts[rname_rand] = 0;
    }
    snp_counts[rname_rand]++;
  }

  rnames.forEach(function(rname, i) {
    var cnt    = snp_counts[rname] || 0;
    var elimit = cnt * 30;
    var snps   = {};
    const fasta        = fastas.result[rname];
    const endpos       = fasta.getEndPos();
    for (var i=0, e=0; i<cnt && e<elimit; i++) {
      var type  = 'SNP';
      var extra = snpselect.random(); 
      var pos   = randomInt(endpos -1) + 1;

      if (!snps[pos]) {
        snps[pos] = 1;
        output(type, pos, rname, 1, extra);
      }
      else {
        e++;
        i--;
      }
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
