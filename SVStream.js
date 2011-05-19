var EventEmitter = require('events').EventEmitter;
var SVConst = require('./SVConst');

/* constructor */
function SVStream(op) {
  this.svs = op.svs || [];
  this.snps = op.snps || [];
  this.pos = 0;
  this.i = 0; // counter for SV
  this.j = 0; // counter for SNP
  this.remnant = '';
}


/* extends */
SVStream.prototype = new EventEmitter();

SVStream.prototype.snpize = function(chunk) {
  var snp = this.snps[this.j];
  var end = this.pos + chunk.length;
  var ins = [];
  var del = [];
  var ret = '';
  while (snp && snp.start+1 <= end) {
    var str = chunk.slice(0, snp.start + 1 - this.pos);
    switch (snp.to) {
    case 1:
    case 2:
    case 3:
      ret += SVConst.makeSNP(snp, str, this.pos);
      break;
    case 4:
      ret += SVConst.makeDeletion({start: snp.start, end: snp.start+1}, str, this.pos);
      del.push(snp.start);
      break;
    case 5:
      ret += SVConst.makeInsertion({
        start    : snp.start,
        end      : snp.start+1,
        flagment : SVConst.BASES[Math.floor(Math.random() * 4)]
      }, str, this.pos);
      ins.push(snp.start);
      break;
    default:
      ret += str;
      break;
    }
    chunk = chunk.slice(snp.start + 1 - this.pos);
    this.pos += str.length;
    this.j++;
    snp = this.snps[this.j];
  }
  return {chunk: ret + chunk, ins: ins, del: del};
};

SVStream.prototype.write = function(data) {
  var chunk = this.remnant + data;
  this.remnant = '';
  var snp = this.snpize(chunk);
  chunk = snp.chunk;
  var sv  = this.svs[this.i];
  var end = this.pos + chunk.length;
  
  /* if all svs were applied */
  if (!sv) {
    this.emitDataWithoutLF(chunk);
    this.pos += chunk.length;
    return;
  }

  /* if the current sv not in this range */
  if (end < sv.start) {
    this.emitDataWithoutLF(chunk);
    this.pos += chunk.length;
    return;
  }

  if (sv.start < this.pos) {
    throw new Error('invalid sv start position.');
  }

  /* execute making sv */
  this.makeSV(sv, chunk, end);
};

SVStream.prototype.emitDataWithoutLF = function(data) {
  this.emit('data', data.split('\n').join(''));
};

SVStream.prototype.pipe = function(stream) {
  this.on('data', function(data) {
    stream.write(data);
  });

  this.on('end', function() {
    stream.end();
  });

  this.on('error', function(e) {
    stream.emit('error', e);
  });
};

SVStream.prototype.end = function() {
  if (this.remnant) {
    var sv = this.svs[this.i];
    var end = this.pos + this.remnant.length;
    this.makeSV(sv, this.remnant, end);
    this.emitDataWithoutLF(this.remnant);
  }
  this.emit('end');
}

SVStream.prototype.makeSV = function(sv, chunk, end) {
  while (sv && sv.end <= end) {
    var str = chunk.slice(0, sv.end - this.pos);
    switch (sv.type) {
    case SVConst.DEL:
      this.emitDataWithoutLF(SVConst.makeDeletion(sv, str, this.pos));
      break;
    case SVConst.INS:
      this.emitDataWithoutLF(SVConst.makeInsertion(sv, str, this.pos));
      break;
    case SVConst.INV:
      this.emitDataWithoutLF(SVConst.makeInversion(sv, str, this.pos));
      break;
    default:
      // err
      break;
    }
    chunk = chunk.slice(sv.end - this.pos);
    this.pos += str.length;
    this.i++;
    sv = this.svs[this.i];
  }
  this.remnant = chunk;
}


/* exports */
module.exports = SVStream;
