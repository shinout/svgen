var SVConst     = require('./SVConst');
var FASTAReader = require('./lib/FASTAReader/FASTAReader');
var dna         = require('./lib/dna');

/* constructor */
function gen(op) {
  op = op || {};
  this.path = op.path || '';

  if (!require('path').existsSync(this.path)) {
    gen.error('"'+this.path + '": No such file.');
    return;
  }

  this.fastas = new FASTAReader(this.path);

  try {
    this.chrom = op.chrom || Object.keys(this.fastas.result)[0];
  }
  catch (e){
    console.log(e);
    gen.error(this.path + 'does not seem to be FASTA format.');
    return;
  }

  if (! this.fastas.result[this.chrom]) {
    gen.error('"'+this.chrom+ '" is not found in ' + this.path);
    return;
  }

  this.fasta      = this.fastas.result[this.chrom];
  this.prelen     = this.fasta.idlen();
  this.linelen    = this.fasta.linelen;
  this.bufferSize = op.bufferSize || 40960;
  this.startIdx   = this.fasta.getStartIndex();
  this.endIdx     = this.fasta.getEndIndex();
  this.svs        = [];
}

/* class variables */
gen.error = function(v) {
  process.stderr.write(v+"\n");
  //process.exit();
}

gen.GenomeStream = require('./GenomeStream');
gen.SVStream     = require('./SVStream');
gen.dna          = dna;
gen.FASTAReader  = FASTAReader;
gen.SVConst      = SVConst;


gen.prototype.checkDuplication = function(idxStart, idxEnd) {
  if (this.svs.length == 0) { return false;}

  var i = 0;
  do {
    if ((this.svs[i].start <= idxStart && idxStart <= this.svs[i].end)
                      ||
        (this.svs[i].start <= idxEnd && idxEnd <= this.svs[i].end) ) {
      return true;
    }
    i++;
  } while (this.svs[i] && this.svs[i].start <= idxEnd);
  return false;
}


gen.prototype.registerSV= function(type, start, len, op) {
  if (typeof type != 'number' ||  type > SVConst.types.length || type < 0) {
    gen.error('SV type error. type must be adequate Number.');
    return;
  }

  if (start <= 0) {
    gen.error('SV['+SVConst.types[type]+']: start must be positive.');
    return;
  }

  if (len <= 0) {
    gen.error('SV['+SVConst.types[type]+']: length must be positive.');
    return;
  }

  var idxStart = this.fasta.getIndex(start);
  var idxEnd   = this.fasta.getIndex(start + len-1) + 1;

  /* check range */
  if (idxStart < this.startIdx || this.endIdx < idxStart) {
    gen.error('SV['+SVConst.types[type]+']: position of '+start+' is out of range.');
    return;
  }
  if (idxEnd < this.startIdx || this.endIdx < idxEnd) {
    gen.error('SV['+SVConst.types[type]+']: position of '+start+' with length('+len+') is out of range.');
    return;
  }


  /* check duplication */
  if (this.checkDuplication(idxStart, idxEnd)) {
    gen.error('SV['+SVConst.types[type]+']: position of '+start+' is duplicated.');
    return;
  }
  var svdata = { type: type, start: idxStart, end: idxEnd};
  if (typeof op == "object") {
    Object.keys(op).forEach(function(k) {
      svdata[k] = op[k];
    });
  }
  
  this.svs.push(svdata);

  this.svs.sort(function(a,b){
    return (a.start > b.start) ? 1 : -1;
  });
}


gen.prototype.registerDel= function(start, len) { this.registerSV(SVConst.DEL, start, len); }

/**
 * @param start: base position (1-base start) left side will of which be inserted.
 * @param length: insertion length
 * @param flagment: insertion flagment. if incompatible with length, then trimmed.
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
  this.registerSV(SVConst.INS, start, len, {flagment: flagment}); 
}
gen.prototype.registerInv= function(start, len) { this.registerSV(SVConst.INV, start, len); }


/**
 * get sv fasta
 *
 */
gen.prototype.genotype = function() {
  var genStream = new gen.GenomeStream({prelen: this.prelen});
  var svStream  = new gen.SVStream({nextream: genStream, svs: this.svs});
  var options   = {
    flags: 'r',
    encoding: 'utf-8',
    bufferSize: this.bufferSize
  };

  var stream = require('fs').createReadStream(this.path, options);

  stream.on('data', function(data){
    svStream.write(data);
  });

  stream.on('error', function(e){
    gen.error(e);
  });

  stream.on('end', function(){
    svStream.end();
  });

  stream.on('close', function(){});
  stream.on('fd', function(fd){});
}

module.exports = gen;
