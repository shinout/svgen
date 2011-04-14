var SVConst = require('./SVConst');
function gen(op) {
  op = op || {};
  this.path = op.path || '';

  if (!require('path').existsSync(this.path)) {
    gen.error('"'+this.path + '": No such file.');
  }

  this.chrom = op.chrom || 'chr1';

  if (typeof op.pos2index == 'function') {
    this.pos2index = op.pos2index;
  }

  if (typeof op.idx2pos == 'function') {
    this.idx2pos = op.idx2pos;
  }

 
  this.prelen = ('>' + this.chrom + '\n').length;

  this.linelen = op.linelen || 50;
  this.bufferSize = op.bufferSize || 40960;

  this.svs = [];
}

gen.error = function(v) {
  process.stderr.write(v+"\n");
}

gen.GenomeStream = require('./GenomeStream');
gen.SVStream = require('./SVStream');

/* static functions */
/* get random DNA flagment */
gen.getRandomFlagment = function(len) {
  var flagment = '';
  for (var i=0; i<len; i++) {
    var p = Math.random();
    if (p > 0.75) {
      flagment += 'A';
    }
    else if (p > 0.5) {
      flagment += 'G';
    }
    else if (p > 0.25) {
      flagment += 'T';
    }
    else {
      flagment += 'C';
    }
  }
  return flagment;
}

gen.idx2pos = function(idx, prelen, linelen) {
  prelen = prelen || 0;
  linelen = linelen || 50;
  return Math.max(0, idx - prelen + 1 - Math.floor((idx - prelen + 1)/(linelen + 1)));
}


gen.pos2index = function(pos, prelen, linelen) {
  prelen = prelen || 0;
  linelen = linelen || 50;
  return prelen + pos -1 + Math.floor( (pos -1)/linelen );
}

/* pos : Nth base -> character index (leftside index of the base)*/
gen.prototype.pos2index = function(pos) {
  return gen.pos2index(pos, this.prelen, this.linelen);
}

gen.prototype.idx2pos = function(idx) {
  return gen.idx2pos(idx, this.prelen, this.linelen);
}


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
  /* check duplication */
  var idxStart = this.pos2index(start);
  var idxEnd   = this.pos2index(start + len -1) + 1; // the last "+ 1" means "right side index"

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
    flagment = gen.getRandomFlagment(len);
  }
  else if (flagment.length < len) {
    flagment += gen.getRandomFlagment(len - flagment.length);
  }
  else if (flagment.length > len) {
    flagment = flagment.slice(0, len);
  }
  this.registerSV(SVConst.INS, start, len, {flagment: flagment}); 
}
gen.prototype.registerInv= function(start, len) { this.registerSV(SVConst.INV, start, len); }

gen.prototype.genotype = function() {
  var genStream = new gen.GenomeStream({prelen: this.prelen});
  var svStream = new gen.SVStream({nextream: genStream, svs: this.svs});
  var options = {
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
