var SVStream = require("../SVStream");
/**
 * test
 **/
function test() {
  const Junjo = require('junjo');
  const $j     = new Junjo();
  const seq    = 'abcdefghijklmnopqrstuvwxyz1234567890';
  const answer = 'dINSERTEDefThijasrqponmlkuvwxyz1234123412341234123412341234123412341234567890';

  [1, 2, 4, 8].forEach(function(num) {
    $j("svstream" + num, function() {
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

      this.absorbData(svstream);

      seqs.forEach(function(s) {
        svstream.write(s);
      });
      svstream.end();
    })
    .post(function(err, out) {
      console.log(out);
      console.assert(out == answer);
    });
  });

  $j.run();
}


if (__filename == process.argv[1]) { test(); }
