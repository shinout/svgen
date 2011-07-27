const SVConst     = require('./SVConst');
const FASTAReader = require('./lib/FASTAReader/FASTAReader');
const dna         = require('./lib/dna');
const nrand       = require('./lib/normal_random');
const XORShift    = require('./lib/xorshift');
const random      = new XORShift(new Date().getTime(), true); // function
const AP          = require('argparser');
const LS          = require('linestream');
const fs          = require('fs');
const spawn       = require('child_process').spawn;

const REPEAT_MEAN   = 40;
const REPEAT_STDDEV = 20;

const numberize = function(v, _default, allow_zero) {
  return ((allow_zero && v == 0) || v === false || v === null || isNaN(Number(v))) ? _default : Number(v);
}

/*
 * node SVGenerator.js
 * 
 * arguments
 *    1: bed file name   (required)
 *    2: fasta file name (required)
 *
 * options
 *  -c, --chrom: chrom name
 *  -s, --svchrom: sv chrom name
 *
 * output
 * fasta with sv to stdout
 *
 */
function main() {
  const p = new AP().addValueOptions(['chrom', 'c', 'name', 'n', 'json', 'exename']).addOptions(['nonstop', 'test']).parse();


  function showUsage() {
    const cmd = p.getOptions('exename') || (process.argv[0] + ' ' + require('path').basename(process.argv[1]));
    console.error('[usage]');
    console.error('\t' + cmd + ' <bed file> <fasta file>');
    console.error('[options]');
    console.error('\t' + '--test\tdryrun. default: null');
    console.error('\t' + '--nonstop\t even if there\'s no registered SVs, execute making a new sequence. defualt: null');
    console.error('\t' + '--chrom\t sequence id in the given file to make SVs and SNPs, default: all (all sequences in the given file)');
    console.error('\t' + '--name | -n\t new chromosome name with SVs and SNPs, default: (the same as original sequences)');
    console.error('\t' + '--json <json file>\t fasta summary file to shortcut calculation.');
    console.error('[bed file columns]');
    console.error('\trname\tstart-position\tend-position\tSVtype(DEL|INS|INV|DUP|TRA|SNP)\tlength');
  }


  var bed = p.getArgs(0);
  if (!bed) {
    showUsage();
    process.exit();
  }

  var fastafile = p.getArgs(1);
  if (!fastafile) {
    showUsage();
    process.exit();
  }
  if (! require('path').existsSync(fastafile)) {
    console.error(fastafile+ ': No such file.');
    process.exit();
  }

  var chrom   = p.getOptions('chrom') || p.getOptions('c');

  var svchrom = p.getOptions('name') || p.getOptions('n') || chrom;

  var svgen   = new gen({path: fastafile, chrom: chrom, svchrom: svchrom, json: p.getOptions('json')});
  if (!svgen.valid) {
    process.exit();
  }

  var bedresult = svgen.registerSVFromBEDFile(bed);
  if (!bedresult && !p.getOptions('nonstop')) process.exit();

  const dryrun = p.getOptions('test');
  var ret = svgen.genotype(null, dryrun);
  
  if (dryrun) console.error(ret);

}



/*
 * gen:
 * @param {Object} op
 *    {String} path       : original fasta file name
 *    {Object} freader    : fasta reader object
 *    {String} chrom      : fasta id (if null, then set first fasta id)
 *    {String} svchrom    : sv fasta id (if null, then set to the same as chrom)
 *    {Number} bufferSize : for each buffer size to get fasta data (default 40960)
 *    {String} json       : fasta summary file to shortcut calculation
 *
 * @return
 *  {gen} object.
 */
function gen(op) {
  op = op || {};
  this.path = (function() {
    if (op.path) return op.path;
    return (op.freader instanceof FASTAReader) ? op.freader.fpath : '';
  })();

  if (!(op.rfeader instanceof FASTAReader) && !require('path').existsSync(this.path)) {
    console.error('"'+this.path + '": No such file.');
    return;
  }

  var json = op.json && require('path').existsSync(op.json) ? JSON.parse(fs.readFileSync(op.json).toString()) : false;
  this.fastas = (op.freader instanceof FASTAReader) ? op.freader : new FASTAReader(this.path, json);

  try {
    this.chrom = op.chrom || Object.keys(this.fastas.result)[0];
  }
  catch (e){
    console.error(e);
    console.error(this.path + 'does not seem to be FASTA format.');
    return;
  }

  if (! this.fastas.result[this.chrom]) {
    console.error('"'+this.chrom+ '" is not found in ' + this.path);
    return;
  }


  this.svchrom    = op.svchrom || this.chrom;
  this.fasta      = this.fastas.result[this.chrom];
  this.prelen     = this.fasta.idlen();
  this.linelen    = this.fasta.linelen;
  this.bufferSize = op.bufferSize || 40960;
  this.startIdx   = this.fasta.getStartIndex();
  this.endIdx     = this.fasta.getEndIndex();
  this.svs        = [];
  this.snps       = [];
  this.valid      = true;
}

gen.SVStream     = require('./SVStream');
gen.dna          = dna;
gen.FASTAReader  = FASTAReader;
gen.SVConst      = SVConst;


/**
 * check duplication between registered SVs.
 * (internal)
 */
gen.prototype.checkDuplication = function(idxStart, idxEnd) {
  const len = this.svs.length;
  if (len == 0) { return false;}

  var i = 0;
  while (i == 0 || this.svs[i-1]) {
    var start = (i >= 1) ? this.svs[i-1].end: -1;
    var end   = (i == len) ? this.endIdx - this.startIdx +1 : this.svs[i].start;
    if (end <= idxStart) {
      i++;
      continue;
    }

    return !( start < idxStart && idxEnd < end);
  }
  return false;
}

/**
 * register SV of the given type
 * (internal)
 * @return {Boolean} : succeed or not
 */
gen.prototype.registerSV= function(type, start, len, op) {
  start = Number(start);
  len   = Number(len);
  if (typeof type != 'number' ||  type > SVConst.types.length || type < 0) {
    console.error('SV type error. type must be adequate Number.');
    return false;
  }

  if (start <= 0) {
    console.error('SV['+SVConst.types[type]+']: start must be positive.');
    return false;
  }

  if (len <= 0) {
    console.error('SV['+SVConst.types[type]+']: length must be positive.');
    return false;
  }

  var idxStart = this.fasta.getIndex(start) - this.startIdx;
  var idxEnd   = this.fasta.getIndex(start + len-1) + 1 - this.startIdx;

  /* check range */
  if (idxStart < 0 || this.endIdx - this.startIdx < idxStart) {
    console.error('SV['+SVConst.types[type]+']: position of '+start+' is out of range.');
    return false;
  }
  if (idxEnd < 0 || this.endIdx - this.startIdx < idxEnd) {
    console.error('SV['+SVConst.types[type]+']: position of '+start+' with length('+len+') is out of range.');
    return false;
  }

  /* check duplication */
  if (this.checkDuplication(idxStart, idxEnd)) {
    console.error('SV['+SVConst.types[type]+']: position of '+start+' is duplicated.');
    return false;
  }

  /* check NNNN */
  if (this.fastas.hasN(this.chrom, start, len)) {
    console.error('SV['+SVConst.types[type]+']: position of '+start+' is region of NNNNNNNNNN.');
    return false;
  }


  var svdata = { type: type, start: idxStart, end: idxEnd, pos: start, len: len};
  if (typeof op == "object") {
    Object.keys(op).forEach(function(k) {
      svdata[k] = op[k];
    });
  }
  
  this.svs.push(svdata);

  this.svs.sort(function(a,b){
    return (a.start > b.start) ? 1 : -1;
  });
  return true;
}

/**
 * register deletion
 *
 * @param start  : base position (1-base start) left side will of which be deleted.
 * @param length : deletion length
 */

gen.prototype.registerDel= function(start, len) { return this.registerSV(SVConst.DEL, start, len); }


/**
 * register insertion
 *
 * @param start    : base position (1-base start) left side will of which be inserted.
 * @param length   : insertion length
 * @param fragment : insertion fragment. if incompatible with length, then trimmed.
 */
gen.prototype.registerIns= function(start, len, fragment) {
  if (!fragment) {
    fragment = dna.getRandomFragment(len);
  }
  else if (fragment.length < len) {
    fragment += dna.getRandomFragment(len - fragment.length);
  }
  else if (fragment.length > len) {
    fragment = fragment.slice(0, len);
  }
  return this.registerSV(SVConst.INS, start, len, {fragment: fragment}); 
}

/**
 * register tandem duplication
 *
 * @param start      : base position (1-base start) left side will of which be tandem-duplicated.
 * @param length     : duplication length
 * @param repeat_num : the number of repeats.
 */
gen.prototype.registerDup = function(start, len, repeat_num) {
  var repeat_num = gen.getTandemDuplicationRepeatNumber(repeat_num);
  return this.registerSV(SVConst.DUP, start, len, {repeat_num: repeat_num});
}

gen.getTandemDuplicationRepeatNumber = function(repeat_num) {
  return numberize(repeat_num, Math.max(Math.floor(nrand(REPEAT_MEAN, REPEAT_STDDEV, random)),2) , false);
}

/**
 * register inversion
 *
 * @param start  : base position (1-base start) left side will of which be inverted.
 * @param length : inversion length
 * @return {Boolean} : succeed or not
 */
gen.prototype.registerInv= function(start, len) { return this.registerSV(SVConst.INV, start, len); }


/**
 * register SV from BED data
 *
 * @param bed        : bed file name
 * @return {Boolean} : succeed or not
 */
gen.prototype.registerSVFromBEDFile = function(bed) {
  if (! require('path').existsSync(bed)) {
    console.error(bed + ': No such file.');
    return false;
  }

  var lines = fs.readFileSync(bed, 'utf-8').split('\n');
  var ret = false;
  lines.forEach(function(line) {
    if (!line || line.charAt(0) == '#') return;
    var svinfo = line.split('\t');

    if (svinfo[0] != this.chrom) {
      //console.error(svinfo[0] + ' : different chromosome type. skip registration.');
      return;
    }
    var svtype = (svinfo[3] == 'SNP') ? 'SNP' : SVConst[svinfo[3]];
    var result;
    switch (svtype) {
      case 'SNP':
        result = this.registerSNP(svinfo[1], svinfo[5]);
        break;
      case SVConst.DEL:
        result = this.registerDel(svinfo[1], svinfo[4]);
        break;
      case SVConst.INS: 
      case SVConst.TRA:  // the same process as INS
        result = this.registerIns(svinfo[1], svinfo[4], svinfo[5]);
        break;
      case SVConst.INV: 
        result = this.registerInv(svinfo[1], svinfo[4]);
        break;
      case SVConst.DUP: 
        result = this.registerDup(svinfo[1], svinfo[4], svinfo[5]);
        break;
      default:
        result = false;
        break;
    }
    ret = (result || ret);
  }, this);

  return ret;
}

/**
 * register SNP
 *
 * @param {Number} pos : base position (1-base start) of SNP.
 * @param {Number} to  : 1: A->C->G->T->A  2: A->G->A C->T->C 3: A->T->G->C->A 4: del 5: INS(random)
 * @return {Boolean}   : succeed or not
 */
gen.prototype.registerSNP = function(pos, to) { 
  to = Number(to);
  if (pos <= 0) {
    console.error('SNP: start must be positive.');
    return false;
  }

  if (to < 1 || 5 < to) {
    console.error('SNP: "to" must be in 1 to 5.');
    return false;
  }

  var idxStart = this.fasta.getIndex(pos);

  /* check range */
  if (idxStart < this.startIdx || this.endIdx < idxStart) {
    console.error('SNP: position of '+ pos +' is out of range.');
    return false;
  }

  var snpdata = { to: to, start: idxStart - this.startIdx, pos: pos};
  this.snps.push(snpdata);

  this.snps.sort(function(a,b){
    return (a.start > b.start) ? 1 : -1;
  });
  return true;
}



/**
 * get sv fasta to stdout.
 */
gen.prototype.genotype = function(wstream, dryrun) {

  var options   = {
    flags      : 'r',
    encoding   : 'utf-8',
    bufferSize : this.bufferSize,
    start      : this.startIdx,
    end        : this.endIdx -1,
  };

  if (dryrun) return {svs: this.svs, snps: this.snps, options: options};

  var rstream = fs.createReadStream(this.path, options);

  var out;

  if (typeof wstream == 'string') {
    out = fs.createWriteStream(wstream);
  }
  else {
    out = (wstream &&
           typeof wstream == 'object' && 
           typeof wstream.write == 'function' &&
           typeof wstream.end == 'function' &&
           typeof wstream.on == 'function'
          )
      ? wstream 
      : process.stdout; 
  }

  
  out.write('>' + this.svchrom + '\n');


  var svstream = new gen.SVStream({svs: this.svs, snps: this.snps}); 

  // pipe
  rstream.on('data', function(d) {
    svstream.write(d);
  });
  rstream.on('end', function() {
    svstream.end();
  });

  var linelen = this.linelen;
  var fold = spawn('fold', ['-w', linelen]);

  svstream.pipe(fold.stdin);

  fold.stdout.on('data', function(data) {
    out.write(data.toString());
  });

  fold.stdout.on('end', function() {
    out.end();
  });

  //fold.stdout.pipe(out);

  svstream.on('error', function(e){
    console.error(e);
  });
}

module.exports = gen;


if (__filename == process.argv[1]) {
  main();
}
