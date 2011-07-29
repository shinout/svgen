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
const dna          = require('./lib/dna');
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


/**
 * test
 **/
function test() {
  const WF     = require('./lib/workflow/wflight');
  const wf     = new WF();
  const seq    = 'abcdefghijklmnopqrstuvwxyz1234567890';
  const answer = 'dINSERTEDefThijasrqponmlkuvwxyz1234123412341234123412341234123412341234567890';

  wf.addCommands([1, 2, 4, 8].map(function(num) {
    return function() {
      var len  = seq.length;
      var l    = Math.floor(len / num);
      var seqs = [];
      var st   = 0;
      for (var i=0; i<num-1; i++) {
        seqs.push(seq.substr(st, l));
        st += l;
      }
      seqs.push(seq.substr(st));

      const svstream = new SVStream([
        [1, 3, 'DEL', null],
        [5, 5, 'INS', 'INSERTED'],
        [7, 7, 'SNP', 1],
        [11, 20, 'INV', null],
        [27, 30, 'DUP', 10],
      ]);

      var result = '';
      svstream.on('data', function(data) {
        result += data;
      });

      svstream.on('end', function() {
        console.log(result);
        console.log(result == answer);
        wf.next();
      });

      seqs.forEach(function(s) {
        svstream.write(s);
      });

      svstream.end();
    };
  }));

  wf.run();
}

module.exports = SVStream;
if (__filename == process.argv[1]) { test(); }
