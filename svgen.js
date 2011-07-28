const SVConst     = require('./SVConst');
const FASTAReader = require('./lib/FASTAReader/FASTAReader');
const dna         = require('./lib/dna');
const nrand       = require('./lib/normal_random');
const XORShift    = require('./lib/xorshift');
const random      = new XORShift(new Date().getTime(), true); // function
const SortedList  = require('./lib/SortedList');
const AP          = require('argparser');
const LS          = require('linestream');
const fs          = require('fs');
const pa          = require('path');
const spawn       = require('child_process').spawn;
const numberize = function(v, _default, allow_zero) {
  return ((allow_zero && v == 0) || v === false || v === null || isNaN(Number(v))) ? _default : Number(v);
}


/*
 * node svgen.js
 * 
 * arguments
 *    1: bed file name   (required)
 *    2: fasta file name (required)
 *
 * options
 *  -j, --json: sv chrom name
 *
 * output
 * fasta with sv to stdout
 *
 */
function main() {
  const p = new AP().addValueOptions(['json', 'exename']).addOptions([]).parse();


  function showUsage() {
    const cmd = p.getOptions('exename') || (process.argv[0] + ' ' + pa.basename(process.argv[1]));
    console.error('[usage]');
    console.error('\t' + cmd + ' <bed file> <fasta file>');
    console.error('[options]');
    console.error('\t' + '--json|-j <json file>\t fasta summary file to shortcut calculation.');
    console.error('[bed file columns]');
    console.error('\trname\tstart-position\tend-position\tSVtype(DEL|INS|INV|DUP|TRA|SNP)\tlength\textra info');
  }

  /* check requirements */
  var bedfile = p.getArgs(0);
  if (!bedfile) {
    showUsage();
    process.exit();
  }

  var fastafile = p.getArgs(1);
  if (!fastafile) {
    showUsage();
    process.exit();
  }
  if (! pa.existsSync(fastafile)) {
    console.error(fastafile+ ': No such file.');
    process.exit();
  }

  var svgen     = new SVGen(fastafile, {json: json});
  var bedstream = new LS(bedfile, {trim: true});

  /* register from BED */
  bedstream.on('data', function(line) {
    if (!line || line.charAt(0) == '#') return;
    var svinfo = line.split('\t');
    if (svinfo.length < 6) return;

    var rname  = svinfo[0];
    var start  = svinfo[1];
    //var end    = svinfo[2];
    var type   = svinfo[3];
    var len    = svinfo[4];
    var extra  = svinfo[5];
    try {
      svgen.register(rname, start, len, type, extra);
    }
    catch (e) {
      console.error(e);
    }
  });

  /* run */
  bedstream.on('end', function() {
    svgen.run();
  });
}

// options for SortedList
const listOptions = {
  SV: {
    filter: function(val, pos) {
      return (this.arr[pos]   == null || (this.arr[pos]   != null && this.arr[pos][1]  <  val[0])) 
        &&   (this.arr[pos+1] == null || (this.arr[pos+1] != null && val[1] < this.arr[pos+1][0]));
    },
    compare: function(a, b) {
      if (a == null) return -1;
      if (b == null) return  1;
      var c = a[0] - b[0];
      return (c > 0) ? 1 : (c == 0)  ? 0 : -1;
    }
  },

  SNP: {
    filter : function(val, pos) {
      return (this.arr[pos] == null) || (this.arr[pos] != null && this.arr[pos] != val);
    }
  }
};

function SVGen(fasta, options) {
  this.fastas = (fasta instanceof FASTAReader) ? fasta : (function() {
    if (pa.existsSync(fasta)) {
      return new FASTAReader(fasta, options.json);
    }
    else {
      throw new Error(fasta + ' : no such file.');
    }
  })();
  if (!this.fastas) { return false; }

  this.regions = {};

  Object.keys(this.fastas.result).forEach(function(rname) {
    //this.regions[rname] = {SV: new SortedList(null ,sv_options), SNP: new SortedList(null, snp_options)};
    this.regions[rname] = {SV: null, SNP: null};
  }, this);
}

SVGen.prototype.register = function(rname, start, len, type, extra) {
  if (! this.fastas[rname] instanceof FASTAReader.FASTA) {
    throw new Error(rname + ' : invalid rname.');
  }

  if (typeof SVGen.valid[type] != 'function') {
    throw new Error(type + ' : invalid type of SV.');
  }

  if (!SVGen.validRange(this.fastas, rname, start, len)) {
    throw new Error('invalid range.');
  }

  if (!SVGen.noNRegion(this.fastas, rname, start, len)) {
    throw new Error('region of NNN...');
  }

  var snp_or_sv = (type == 'SNP') ? 'SNP' : 'SV';
  if (this.regions[rname][snp_or_sv] == null) {
    this.regions[rname][snp_or_sv] = new SortedList(null, listOptions[snp_or_sv]);
  }
  var regions = this.regions[rname][snp_or_sv];

  if (!SVGen.valid[type](this.fastas, rname, start, len, extra)) {
    throw new Error('invalid format for ' + type + '.');
  }

  var bool = regions.insert([
    this.fastas.result[rname].getIndex(start),
    this.fastas.result[rname].getIndex(start + len),
    type,
    extra
  ]);

  if (!bool) {
    throw new Error('overlapped region.');
  }
};


SVGen.prototype.run = function() {
};

/*** static functions ***/
/**
 *
 * return boolean
 */
SVGen.validRange = function(fastas, rname, start, len) {
  var fasta = fastas.result[rname];
  return (1 <= start) && ((start + len) <= fasta.getEndPos()); 
};

/**
 *
 * return boolean
 */
SVGen.noNRegion = function(fastas, rname, start, len) {
  return !fastas.hasN(rname, start, len);
};

// specific validation of each SV
SVGen.valid = {
  DEL: function(fastas, rname, start, len, extra) {
    return true;
  },

  INS: function(fastas, rname, start, len, extra) {
    return extra.length == len;
  },

  INV: function(fastas, rname, start, len, extra) {
    return true;
  },

  DUP: function(fastas, rname, start, len, extra) {
    return true;
  },

  TRA: function(fastas, rname, start, len, extra) {
    var trinfo = extra.split(':');
    if (trinfo.length < 3) return false;
    var trname = trinfo[0];
    var tstart = Number(trinfo[1]);
    var tlen   = Number(trinfo[2]);
    if (! fastas.result[trname] instanceof FASTAReader.FASTA) return false;
    if (!SVGen.validRange(fastas, trname, tstart, tlen)) return false;
    if (!SVGen.noNRegion(fastas, trname, tstart, tlen)) return false;
    return true; 
  },

  SNP: function(fastas, start, len, extra) {
    return ([1,2,3,4,5].indexOf(extra) >= 0);
  }
};

module.exports = SVGen;
if (__filename == process.argv[1]) { main(); }
