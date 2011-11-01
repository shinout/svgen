#!/usr/bin/env node
const SVStream    = require('./SVStream');
const FASTAReader = require('fastareader');
const dna         = require('dna');
const SortedList  = require('sortedlist');
const AP          = require('argparser');
const LS          = require('linestream');
const Junjo       = require('junjo');
const fs          = require('fs');
const pa          = require('path');
const spawn       = require('child_process').spawn;
const numberize = function(v, _default, allow_zero) {
  return ((allow_zero && v == 0) || v === false || v === null || isNaN(Number(v))) ? _default : Number(v);
};


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
  const p = new AP().addValueOptions(['j', 'json', 'exename', 'r', 'rnames']).addOptions([]).parse();


  function showUsage() {
    const cmd = p.getOptions('exename') || (process.argv[0] + ' ' + pa.basename(process.argv[1]));
    console.error('[usage]');
    console.error('\t' + cmd + ' <bed file> <fasta file>');
    console.error('[options]');
    console.error('\t' + '--rnames|-r <sequence id1>[,sequence_id2,...]\t sequence ids to use. default: null (use all rnames in a fasta file.)');
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

  const rnames = (p.getOptions('rnames', 'r')) ? p.getOptions('rnames', 'r').split(',') : null;
  const svgen     = new SVGen(fastafile, {json: json, rnames: rnames});
  const bedstream = new LS(bedfile, {trim: true});

  /* register from BED */
  bedstream.on('data', function(line) {
    if (!line || line.charAt(0) == '#') return;
    var svinfo = line.trim().split('\t');
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
    //console.log(svgen.getRegions("chr1"));
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


/**
 * constructor
 * 
 * @param (String) fasta: file name of FASTA to use.
 * @param (Object) options: 
 *   (Object) json  : json data of fasta file.
 *   (Array) rnames : sequence ids to use.
 *   (Function) callback : called when all sequences are written (emitted).
 *
 **/
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

  this.rnames = [];

  if (options.rnames instanceof Array) {
    options.rnames.forEach(function(rname) {
      if (this.fastas.result[rname] instanceof FASTAReader.FASTA) {
        this.rnames.push(rname);
      }
      else {
        console.error('rname "' + rname + '" is not in FASTA file, so skipped this name.');
      }
    }, this);
  }
  else {
    this.rnames = Object.keys(this.fastas.result);
  }

  this.callback = options.callback;

  this.regions = {};

  Object.keys(this.fastas.result).forEach(function(rname) {
    this.regions[rname] = {SV: new SortedList(null ,listOptions.SV), SNP: new SortedList(null, listOptions.SNP)};
    //this.regions[rname] = {SV: null, SNP: null};
  }, this);
}

Object.defineProperty(SVGen.prototype, 'callback', {
  get: function() {
    return this._callback || function() {};
  },
  set: function(v) {
    if (typeof v == 'function') this._callback = v;
  }
});

/**
 * register SV/SNP
 *
 * @param (String) rname : rname to put SV/SNP
 * @param (Integer) start : start position of SV/SNP
 * @param (Integer) len: length of SV/SNP
 * @param (String) type: type of SV/SNP (must be one of DEL, INS, INV, DUP, TRA, SNP)
 * @param (mixed) extra: extra information for each SV/SNP
 *  switch (type)
 *   case DEL, INV : extra is not used.
 *   case INS      : extra is a sequence insert.
 *   case DUP      : extra is the number of repeat.
 *   case TRA      : extra is the position information from which to insert in the following format.
 *     "rname:start_position"
 *   case SNP      : extra is one of [1,2,3], which represents the type of alteration.
 *
 * @param (boolean) suspend: suspend registration until the next registration is succeeded.
 *@throws Error
 *@returns
 **/
SVGen.prototype.register = function(rname, start, len, type, extra, suspend) {
  try {
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

    var data = SVGen.valid[type].call(this, this.fastas, rname, start, len, extra);
    if (!data) {
      throw new Error('invalid format for ' + type + '.\n' + [rname, start, len, extra].join('\n') );
    }

    var snp_or_sv = (type == 'SNP') ? 'SNP' : 'SV';
    /*
    if (this.regions[rname][snp_or_sv] == null) {
      this.regions[rname][snp_or_sv] = new SortedList(null, listOptions[snp_or_sv]);
    }
    */
    var regions = this.regions[rname][snp_or_sv];
    var pos = regions.insert(data);

    if (pos === false) {
      throw new Error('overlapped region.');
    }
    if (suspend) this.pending = { regions: regions, pos: pos };
  }
  catch (e) {
    if (this.pending) this.pending.regions.remove(this.pending.pos);
    throw e;
  }
  finally {
    if (!suspend) this.pending = null;
    return pos;
  }
};

SVGen.prototype.getRegions = function(rname, type) {
  if (!type) type = "SV";
  return this.regions[rname][type];
};

/**
 * get a sequence with SV in FASTA format.
 * @param (mixed) wstream : writable stream. When string given, a file the name of which is wstream, is created.
 *                          If null, then wstream = process.stdout
 * @returns
 **/
SVGen.prototype.run = function(wstream) {
  const rnames = this.rnames;
  wstream = (function() {
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

  var $j = new Junjo({noTimeout: true});

  rnames.forEach(function(rname, num) {
    $j("sv." + rname, function() {
      var fasta   = that.fastas.result[rname];
      var rstream = fs.createReadStream(that.fastafile, {
        flags      : 'r',
        encoding   : 'utf-8',
        bufferSize : 40960,
        start      : fasta.getStartIndex(),
        end        : fasta.getEndIndex() -1
      });
      var snpstream = new SVStream(that.regions[rname].SNP.toArray());
      var svstream  = new SVStream(that.regions[rname].SV.toArray());
      var fold      = spawn('fold', ['-w', fasta.linelen]);

      // rstream -> snpstream -> svstream -> fold.stdin -> wstream

      snpstream.pipe(svstream);
      svstream.pipe(fold.stdin);


      rstream.setEncoding("utf8");
      rstream.on('data', function(data) {
        snpstream.write(data.split('\n').join(''));
      });

      rstream.on('end', function() {
        snpstream.end();
      });

      fold.stdout.setEncoding('utf8');

      this.absorb(fold.stdout, "data", function(data) {
        wstream.write(data);
      });

      // after previous fasta, put LF
      if (num) {
        wstream.write('\n');
      }

      wstream.write('>' + rname + '\n');

    }).after();
  });

  $j(function() {
    wstream.end();
    wstream.on('close', this.cb);
  })
  .after();

  $j.on("end", that.callback);

  $j.run();
};

/*** static functions ***/
/**
 * @param (FASTAReader) fastas
 * @param (String) rname
 * @param (Integer) start 
 * @param (Integer) len 
 * @returns boolean
 */
SVGen.validRange = function(fastas, rname, start, len) {
  var fasta = fastas.result[rname];
  return (1 <= start) && ((start + len) <= fasta.getEndPos()); 
};

/**
 * @param (FASTAReader) fastas
 * @param (String) rname
 * @param (Integer) start 
 * @param (Integer) len 
 * @returns boolean
 */
SVGen.noNRegion = function(fastas, rname, start, len) {
  return !fastas.hasN(rname, start, len);
};

/**
 * specific validation of each SV
 * of each function, 
 * @param (FASTAReader) fastas
 * @param (String) rname
 * @param (Integer) start 
 * @param (Integer) len 
 * @returns valid data for SVStream.
 **/
SVGen.valid = {
  DEL: function(fastas, rname, start, len, extra) {
    var ret = [start, start + len -1, 'DEL', null];
    if (Array.isArray(extra)) ret.push(extra);
    return ret;
  },

  INS: function(fastas, rname, start, len, extra) {
    return (extra.length == len) ? [start, start, 'INS', extra] : false;
  },

  INV: function(fastas, rname, start, len, extra) {
    return [start, start + len -1, 'INV', null];
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

    try {
      // pending the registration of DEL until INS is successfully registered.
      this.register(trname, tstart, len, 'DEL', [rname, start], true); 
    }
    catch (e) {
      return false;
    }
    return [start, start, 'INS', fastas.fetch(trname, tstart, len), trinfo];  // converted to insertion info.
  },

  SNP: function(fastas, rname, start, len, extra) {
    extra = Number(extra);
    return ([1,2,3].indexOf(extra) >= 0) ? [start, start, 'SNP', extra] : false;
  }
};

module.exports = SVGen;
if (process.argv[1].match('/([^/]+?)(\.js)?$')[1] == __filename.match('/([^/]+?)(\.js)?$')[1]) main();
