const SVConst     = require('./SVConst');
const FASTAReader = require('./lib/FASTAReader/FASTAReader');
const dna         = require('./lib/dna');
const AP          = require('argparser');
const LS          = require('linestream');
const fs          = require('fs');
const spawn       = require('child_process').spawn;

/*
 * node SVGenerator.js
 * 
 * arguments
 *    1: tsv file name   (required)
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
  function showUsage() {
    gen.error('Usage: ' + require('path').basename(process.argv[0]) + ' ' + process.argv[1] + ' [--test] [--nonstop] [-c|--chrom <chrom name>] [-n|--name <sv chrom name>] <tsv file> <fasta file>');
    gen.error('tsv file columns: SVtype(DEL|INS|INV|SNP)\tstart-position\tlength');
  }

  var p = new AP().addValueOptions(['chrom', 'c', 'name', 'n']).addOptions(['nonstop', 'test']).parse();

  var tsv = p.getArgs(0);
  if (!tsv) {
    showUsage();
    process.exit();
  }

  var fastafile = p.getArgs(1);
  if (!fastafile) {
    showUsage();
    process.exit();
  }
  if (! require('path').existsSync(fastafile)) {
    gen.error(fastafile+ ': No such file.');
    process.exit();
  }

  var chrom   = p.getOptions('chrom') || p.getOptions('c');

  var svchrom = p.getOptions('name') || p.getOptions('n') || chrom;

  var svgen   = new gen({path: fastafile, chrom: chrom, svchrom: svchrom});
  if (!svgen.valid) {
    process.exit();
  }

  var tsvresult = svgen.registerSVFromTSVFile(tsv);
  if (!tsvresult && !p.getOptions('nonstop')) process.exit();

  const dryrun = p.getOptions('test');
  var ret = svgen.genotype(null, dryrun);
  
  if (dryrun) console.log(ret);

}



/*
 * gen:
 * @param {Object} op
 *    {String} path       : original fasta file name
 *    {Object} freader    : fasta reader object
 *    {String} chrom      : fasta id (if null, then set first fasta id)
 *    {String} svchrom    : sv fasta id (if null, then set to the same as chrom)
 *    {Number} bufferSize : for each buffer size to get fasta data (default 40960)
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
    gen.error('"'+this.path + '": No such file.');
    return;
  }

  this.fastas = (op.freader instanceof FASTAReader) ? op.freader : new FASTAReader(this.path);

  try {
    this.chrom = op.chrom || Object.keys(this.fastas.result)[0];
  }
  catch (e){
    console.error(e);
    gen.error(this.path + 'does not seem to be FASTA format.');
    return;
  }

  if (! this.fastas.result[this.chrom]) {
    gen.error('"'+this.chrom+ '" is not found in ' + this.path);
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
    gen.error('SV type error. type must be adequate Number.');
    return false;
  }

  if (start <= 0) {
    gen.error('SV['+SVConst.types[type]+']: start must be positive.');
    return false;
  }

  if (len <= 0) {
    gen.error('SV['+SVConst.types[type]+']: length must be positive.');
    return false;
  }

  var idxStart = this.fasta.getIndex(start) - this.startIdx;
  var idxEnd   = this.fasta.getIndex(start + len-1) + 1 - this.startIdx;

  /* check range */
  if (idxStart < 0 || this.endIdx - this.startIdx < idxStart) {
    gen.error('SV['+SVConst.types[type]+']: position of '+start+' is out of range.');
    return false;
  }
  if (idxEnd < 0 || this.endIdx - this.startIdx < idxEnd) {
    gen.error('SV['+SVConst.types[type]+']: position of '+start+' with length('+len+') is out of range.');
    return false;
  }


  /* check duplication */
  if (this.checkDuplication(idxStart, idxEnd)) {
    gen.error('SV['+SVConst.types[type]+']: position of '+start+' is duplicated.');
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
 * @param flagment : insertion flagment. if incompatible with length, then trimmed.
 */
gen.prototype.registerIns= function(start, len, flagment) {
  if (!flagment) {
    flagment = dna.getRandomFlagment(len);
  }
  else if (flagment.length < len) {
    flagment += dna.getRandomFlagment(len - flagment.length);
  }
  else if (flagment.length > len) {
    flagment = flagment.slice(0, len);
  }
  return this.registerSV(SVConst.INS, start, len, {flagment: flagment}); 
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
 * register SV from TSV data
 *
 * @param tsv        : tsv file name
 * @return {Boolean} : succeed or not
 */
gen.prototype.registerSVFromTSVFile = function(tsv) {
  if (! require('path').existsSync(tsv)) {
    gen.error(tsv + ': No such file.');
    return false;
  }

  var lines = fs.readFileSync(tsv, 'utf-8').split('\n');
  var ret = false;
  lines.forEach(function(line) {
    if (!line || line.charAt(0) == '#') return;
    var svinfo = line.split('\t');

    if (svinfo[3] != this.chrom) {
      //gen.error(svinfo[3] + ' : different chromosome type. skip registration.');
      return;
    }
    var svtype = (svinfo[0] == 'SNP') ? 'SNP' : SVConst[svinfo[0]];
    var result;
    switch (svtype) {
      case 'SNP':
        result = this.registerSNP(svinfo[1], svinfo[2]);
        break;
      case SVConst.DEL:
        result = this.registerDel(svinfo[1], svinfo[2]);
        break;
      case SVConst.INS: 
        result = this.registerIns(svinfo[1], svinfo[2]);
        break;
      case SVConst.INV: 
        result = this.registerInv(svinfo[1], svinfo[2]);
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
    gen.error('SNP: start must be positive.');
    return false;
  }

  if (to < 1 || 5 < to) {
    gen.error('SNP: "to" must be in 1 to 5.');
    return false;
  }

  var idxStart = this.fasta.getIndex(pos);

  /* check range */
  if (idxStart < this.startIdx || this.endIdx < idxStart) {
    gen.error('SNP: position of '+ pos +' is out of range.');
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
    out = (typeof wstream == 'object' && 
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
    gen.error(e);
  });
}

/**
 * output error
 *
 */
gen.error = function(v) {
  process.stderr.write(v+"\n");
}


module.exports = gen;


if (__filename == process.argv[1]) {
  main();
}
