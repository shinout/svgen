var EventEmitter = require('events').EventEmitter;
var SVConst = require('./SVConst');

/* constructor */
function SVStream(op) {
  this.svs = op.svs || [];
  this.pos = 0;
  this.i = 0;
  this.remnant = '';
}

/* extends */
SVStream.prototype = new EventEmitter();

SVStream.prototype.write = function(data) {
  var chunk = this.remnant + data;
  this.remnant = '';
  var sv = this.svs[this.i];
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
