const SVStream    = require('./SVStream');
const FASTAReader = require('./lib/FASTAReader/FASTAReader');
const dna         = require('./lib/dna');
const nrand       = require('./lib/normal_random');
const XORShift    = require('./lib/xorshift');
const random      = new XORShift(new Date().getTime(), true); // function
const SortedList  = require('./lib/SortedList');
const AP          = require('argparser');
const LS          = require('linestream');
const Flow        = require('./lib/workflow/wflight');
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
    console.error('\trname\tstart-position\tend-position\tSVtype(DEL|INS|INV|DUP|TRA|SNP)\tlength\textra-info');
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

  const jsonfile = p.getOptions('json');

  const json = (function() {
    if (pa.existsSync(jsonfile)) {
      try {
        return JSON.parse(fs.readFileSync(jsonfile).toString());
      }
      catch (e) {
        return null;
      }
    }
    else {
      return null;
    }
  })();

  var svgen     = new SVGen(fastafile, {json: json});
  var bedstream = new LS(bedfile, {trim: true});

  /* register from BED */
  bedstream.on('data', function(line) {
    if (!line || line.charAt(0) == '#') return;
    var svinfo = line.split('\t');
    if (svinfo.length < 6) return;

    var rname  = svinfo[0];
    var start  = Number(svinfo[1]);
    //var end    = svinfo[2];
    var type   = svinfo[3];
    var len    = Number(svinfo[4]);
    var extra  = svinfo[5];
    try {
      svgen.register(rname, start, len, type, extra);
    }
    catch (e) {
      console.error(e.message);
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
    },

    compare: function(a, b) {
      if (a == null) return -1;
      if (b == null) return  1;
      var c = a[0] - b[0];
      return (c > 0) ? 1 : (c == 0)  ? 0 : -1;
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
  this.fastafile = this.fastas.fpath;
  if (!this.fastas) { return false; }

  this.regions = {};

  Object.keys(this.fastas.result).forEach(function(rname) {
    this.regions[rname] = {SV: new SortedList(null ,listOptions.SV), SNP: new SortedList(null, listOptions.SNP)};
    //this.regions[rname] = {SV: null, SNP: null};
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

  var data = SVGen.valid[type](this.fastas, rname, start, len, extra);
  if (!data) {
    throw new Error('invalid format for ' + type + '.');
  }

  var snp_or_sv = (type == 'SNP') ? 'SNP' : 'SV';
  /*
  if (this.regions[rname][snp_or_sv] == null) {
    this.regions[rname][snp_or_sv] = new SortedList(null, listOptions[snp_or_sv]);
  }
  */
  var regions = this.regions[rname][snp_or_sv];
  var bool = regions.insert(data);

  if (!bool) {
    throw new Error('overlapped region.');
  }
};


SVGen.prototype.run = function(wstream) {
  const rnames = Object.keys(this.fastas.result);
  const wf     = new Flow();
  wstream      = (function() {
    if (typeof wstream == 'string') { return fs.createWriteStream(wstream); }
    return  (wstream &&
           typeof wstream == 'object' &&
           typeof wstream.write == 'function' &&
           typeof wstream.end == 'function' &&
           typeof wstream.on == 'function'
          )
      ? wstream 
      : process.stdout; 
  })();


  const that = this;
  wf.addCommands(rnames.map(function(rname) {
    return function() {
      const fasta       = that.fastas.result[rname];
      const rstream     = fs.createReadStream(that.fastafile, {
        flags      : 'r',
        encoding   : 'utf-8',
        bufferSize : 40960,
        start      : fasta.getStartIndex(),
        end        : fasta.getEndIndex() -1
      });
      const snpstream   = new SVStream(that.regions[rname].SNP.toArray());
      const svstream    = new SVStream(that.regions[rname].SV.toArray());
      const fold        = spawn('fold', ['-w', fasta.linelen]);

      snpstream.pipe(svstream);
      svstream.pipe(fold.stdin);

      wstream.write('>' + rname + '\n');
      rstream.on('data', function(data) {
        snpstream.write(data.toString().split('\n').join(''));
      });

      rstream.on('end', function() {
        snpstream.end();
      });

      fold.stdout.on('data', function(data) {
        wstream.write(data.toString());
      });

      fold.stdout.on('end', function() {
        wstream.write('\n');
        wf.next();
      });
    };
  }));
  wf.run();
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

/**
 * specific validation of each SV
 * returns valid data for SVStream.
 * 
 **/
SVGen.valid = {
  DEL: function(fastas, rname, start, len, extra) {
    return [start, start + len -1, 'DEL', null];
  },

  INS: function(fastas, rname, start, len, extra) {
    return (extra.length == len) ? [start, start, 'INS', extra] : false;
  },

  INV: function(fastas, rname, start, len, extra) {
    return [start, start + len -1, 'DEL', null];
  },

  DUP: function(fastas, rname, start, len, extra) {
    extra = Number(extra);
    return (!isNaN(extra) && extra >= 2) ? [start, start + len -1, 'DUP', extra] : false;
  },

  TRA: function(fastas, rname, start, len, extra) {
    var trinfo = extra.split(':');
    if (trinfo.length < 2) return false;
    var trname = trinfo[0];
    var tstart = Number(trinfo[1]);
    if (! fastas.result[trname] instanceof FASTAReader.FASTA) return false;
    if (!SVGen.validRange(fastas, trname, tstart, len)) return false;
    if (!SVGen.noNRegion(fastas, trname, tstart, len)) return false;
    return [start, start, 'INS', fastas.fetch(trname, tstart, len)];  // converted to insertion info.
  },

  SNP: function(fastas, rname, start, len, extra) {
    extra = Number(extra);
    return ([1,2,3].indexOf(extra) >= 0) ? [start, start, 'SNP', extra] : false;
  }
};


SVGen.getRandomFragment = dna.getRandomFragment;
SVGen.version = "1.0.0";
module.exports = SVGen;
if (__filename == process.argv[1]) { main(); }
