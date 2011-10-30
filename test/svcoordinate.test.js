const T           = require('./load.test');
const spawn       = require('child_process').spawn;
const exec        = require('child_process').exec;
const FR          = require('fastareader');
const randomInt   = require('random-tools').randomInt;
const ArrayStream = require('arrayemitter');
const Junjo       = require('junjo');
const $svcoord    = require('../svcoordinate');
const SVGen       = require('../svgen');
const fs          = require('fs');
const dna         = require('dna');

function tst() {
  var $j = new Junjo({run : true});


  $j.start(function() {
    this.$ = {
      SV_NUM      : 100,
      SV_FILENAME : 'sv.fasta',
      COORD_NUM   : 50,
      DEBUG_MODE  : false 
    };
  });

  // $j.assign('fasta', __dirname + '/chrM.fa', 'chrM'); // in future
  $j('fasta', function() {
    return Junjo.multi(__dirname + '/chrM.fa', 'chrM');
  });


  $j('svbedgen', function(fa) {
    var f_svbedgen = __dirname + '/../svbedgen.js';
    var options = {
      svlen   : 200,
      svdev   : 20,
      snprate : 0,
      sv      : this.$.SV_NUM
    };

    var op_str = (function(op) {
      return Object.keys(op).map(function(k) { return '--' + k + ' ' + op[k] }).join(' ');
    })(options);

    exec(['node', f_svbedgen, fa, op_str].join(' '), this.cb);
  })
  .after('fasta')
  .firstError('shift')
  .post(function(out, err) {
    if (err.split('\n').length > 8) throw new Error(err);

    return new ArrayStream(out.split('\n').filter(function(v) {
      return v.charAt(0) != '#';
    }), {pause: true});
  });


  $j('randomCoords', function(fa, rname) {
    var fastas = new FR(fa);
    var fasta  = fastas.result[rname];
    var endpos = fasta.getEndPos();
    var ret = [];
    for (var i=0, l = this.$.COORD_NUM; i<l; i++) {
      var a = randomInt(endpos, 1, "xorshift"), b = a + randomInt(400, 1, "xorshift");
      if (b > endpos) { i--; continue;}
      var strand = (randomInt(1,0,"xorshift")) ? "+" : "-";
      ret.push([rname, a, b, strand, "anno" + i]);
    }

    this.$.coords = ret.sort(function(a, b) {
      return (a[1] > b[1]) ? 1 : -1;
    });


    return new ArrayStream(
      this.$.coords.map(function(v) { return v.join('\t') })
    , {pause: true}
    );
  })
  .after('fasta');


  $j('svcoord', function(coords, svbed, fasta, rname) {
    var $co = $svcoord(this.$.DEBUG_MODE);

    var self = this;
    $co.on('end', function(err, svgen) {
      self.$.svgen = svgen;
    });

    this.absorb($co, 'data', function(data, result) {
      if (!result) return [data];
      result.push(data);
    });

    $co.run(fasta, svbed, coords, null);

  })
  .after('randomCoords', 'svbedgen', 'fasta')
  .firstError('shift');

  $j('svfasta', function() {
    var svgen = this.$.svgen;
    svgen.callback = this.cb;
    svgen.run(fs.createWriteStream(this.$.SV_FILENAME));
  }).after('svcoord');


  $j('check', function(bed, fasta, rname) {
    var fastas   = new FR(fasta);
    var svfastas = new FR(this.$.SV_FILENAME);

    var coords = this.$.coords;
    bed.forEach(function(v) {
      var l = v.length - 1;
      var rn = v[0], start = v[1], end = v[2], strand = v[3];
      var part = [l-5], type = v[l-4], prname = v[l-3], pstart= v[l-2], pend = v[l-1], pstrand = v[l];
      if (type == 'INS') return;

      var newSeq = svfastas.fetch(rname, start, end - start + 1, strand == '-');
      var oriSeq = fastas.fetch(rname, pstart, pend - pstart + 1, pstrand == '-');

      //if (type == 'INV') newSeq = dna.complStrand(newSeq, true);

      if (type == "TRA") console.cyan(rname, start, pstart, end, pend, newSeq, oriSeq, part, type)
      var bool = T.equal(newSeq, oriSeq, part);
      if (!bool) {
        console.yellow("==========================================================");
        console.yellow(v);
        console.yellow("new", newSeq);
        console.yellow("ori", oriSeq);
        console.yellow("==========================================================");
      }
    });
  })
  .after('svcoord', 'fasta', 'svfasta');


  $j.catcher = function(e, args) {
    console.yellow('in', this.label);
    console.red(e.stack);
    this.terminate();
  };

}

tst();
