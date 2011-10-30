/***
 * SVStream
 * argument       : list of SV information
 * expected input : just sequence (String). without meta data or whitespaces. 
 * output         : sequence with SV
 *
 *
 *
 ***/
const EventEmitter = require('events').EventEmitter;
const dna          = require('dna');
const BASES        = ['A', 'C', 'G', 'T'];

/* constructor */
function SVStream(regions) {
  // list of sv info
  this.regions = regions;
  this.total   = regions.length;

  // state
  this.offset  = 0;
  this.i       = 0; // counter for SV
  this.remnant = '';

}

/* extends */
SVStream.prototype = new EventEmitter();

SVStream.prototype.write = function(seq) {
  seq = this.remnant + seq;
  const start = 1;
  var end     = seq.length;

  if (this.total <=  this.i) {
    this.emit('data', seq);
    this.remnant = '';
    return;
  }

  for (this.i; this.i < this.total && start <= end; this.i++) {
    var sv      = this.regions[this.i];
    var svstart = sv[0] - this.offset;
    var svend   = sv[1] - this.offset;
    var type    = sv[2];
    var extra   = sv[3];

    if (start > svstart || svend > end) {
      break;
    }

    // emit non effected region
    this.emit('data', seq.slice(0, svstart -1));

    // emit effected region with SV
    this.emit('data', SVStream.sv(seq.slice(svstart -1, svend), type, extra) );
    
    end -= svend;
    seq = seq.slice(svend);
    this.offset += svend;
  }

  this.remnant = seq;
};

SVStream.prototype.end = function() {
  this.emit('data', this.remnant);
  this.emit('end');
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


SVStream.sv = function(seq, type, extra) {
  switch (type) {
  case 'INS':
  case 'TRA':
    return extra + seq;
  case 'DEL':
    return '';
  case 'INV':
    return dna.complStrand(seq, true);
  case 'DUP':
    var times = Number(extra);
    var ret = '';
    for (var i=0; i<times; i++) {
      ret += seq;
    }
    return ret;
  case 'SNP':
    var idx = BASES.indexOf(seq.toUpperCase());
    return (idx >= 0) ? BASES[(idx + extra) % 4] : seq;
  default:
    return seq;
  }
};


module.exports = SVStream;
