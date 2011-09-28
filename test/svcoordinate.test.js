if (typeof global != 'undefined') require('../lib/Junjo/test/load.test').load(global);

const spawn  = require('child_process').spawn;
const exec   = require('child_process').exec;
const FR     = require('../lib/FASTAReader/FASTAReader');
const random = new require('../lib/xorshift')(new Date().getTime(), true);

function randomInt(max) {
  return Math.floor(random() * max);
}

function tst() {
  var $j = new Junjo({run : true});


  $j.catcher = function(e, args) {
    console.eyellow('in', this.label);
    console.ered(e.stack);
    $j.terminate();
  };

  $j('const', function() {
    return {
      SV_NUM    : 100,
      COORD_NUM : 100
    };
  });

  $j('fasta', function() {
    return Junjo.multi(__dirname + '/chrM.fa', 'chrM');
  });



  $j('svbedgen', function(fa) {
    var f_svbedgen = __dirname + '/../svbedgen.js';
    var options = {
      svlen : 200,
      svdev : 20,
      sv    : $j.results('const').SV_NUM
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
    return out.split('\n').filter(function(v) {
      return v.charAt(0) != '#';
    });
  });

  $j('randomCoords', function(fa, rname) {
    var fasta  = new FR(fa).result[rname];
    var endpos = fasta.getEndPos();
    var ret = [];
    for (var i=0, l = $j.results('const').COORD_NUM; i<l; i++) {
      var a = randomInt(endpos), b = a + randomInt(400);
      if (b > endpos) { i--; continue;}
      ret.push([rname, a, b, "anno" + i]);
    }
    return ret.sort(function(a, b) {
      return (a[1] > b[1]) ? 1 : -1;
    });
  })
  .after('fasta');

  $j('svcoord', function(coords, svbed, fa, rname) {
  })
  .after('randomCoords', 'svbedgen', 'fasta');

}

if (node) tst();
