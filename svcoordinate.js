const fs    = require('fs');
const AP    = require('argparser');
const LS    = require('linestream');
const Junjo = require('./lib/Junjo/Junjo');
const SVGen = require('./svgen');
const cl    = require('./lib/Junjo/lib/termcolor').define()
const pa    = require('path');
const spawn = require('child_process').spawn;
//cl.prod();


function svcoordinate() {
  var $j = new Junjo();

  $j.inputs({
    fasta       : 0,
    bedStream   : 1,
    coordStream : 2,
    json        : 3
  });

  $j('svgen', function(fasta, bedStream, json) {
    var svgen = new SVGen(fasta, {json: json, rnames: null});

    this.absorb(bedStream, 'data', function(line, result) {
      if (!line || line.charAt(0) == '#') return;
      var svinfo = line.split('\t');
      if (svinfo.length < 6) return;

      var rname  = svinfo[0];
      var start  = Number(svinfo[1]);
      //var end    = svinfo[2];
      var type   = svinfo[3];
      var len    = Number(svinfo[4]);
      var extra  = svinfo[5];
      svgen.register(rname, start, len, type, extra);
      return svgen;
    });
    //this.absorbError(bedStream);
  })
  .firstError('shift')
  .after('fasta', 'bedStream', 'json');

  $j('convert', function(coordStream, svgen) {
    this.$.n = 0;

    console.egreen("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
    console.egreen("------------------------     START       ----------------------------------");
    console.egreen("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
    this.absorb(coordStream, 'data', function(line, result, $f) {
      if (!line || line.charAt(0) == '#') return;
      var svinfo = line.split('\t');
      if (svinfo.length < 3) return;

      var rname  = svinfo[0];
      //var regions = svgen.regions[rname];
      var regions = svgen.regions;
      var svpos = regions.SV;
      var snvs  = regions.SV;
      var start  = Number(svinfo[1]);
      var end    = Number(svinfo[2]);
      if (!result) {
        result = { offset: 0, diff: 0, outputs: [] };
      }
      if (result.outputs.length) {
        // console.log(result.outputs)
        result.outputs.forEach(function(r) {
          console.log([r.rname || rname, r.start, r.end, r.part, r.type || '*', r.pstart, r.pend, $f.$.n].join('\t'));
        }, $f);
      }
      $f.$.n++;
      return getNewCoordinate(rname, start, end, regions, result.offset, result.diff);
    });
  })
  .firstError('shift')
  .after('coordStream', 'svgen')
  .next(function(out) {
    //console.log(out);
  });

  $j.catchesAbove(function(e, args) {
    console.ered('[Error] : ' +  e.message + ' in label: ' + cl.yellow(this.label));
    console.ered(e.stack || e);
    this.terminate();
  });

  return $j.prepare();
}

function main() {
  var p  = new AP().addOptions([]).addValueOptions([]).parse();

  var $j = new Junjo();
  var $s = svcoordinate();

  function showUsage() {
    const cmd = p.getOptions('exename') || (process.argv[0] + ' ' + require('path').basename(process.argv[1]));
    console.error('[synopsis]');
    console.error('\t' + cmd + ' <fasta file> <sv bed file> <coordinate bed file>');
    console.error('[options]');
  }

  $j('inputCheck', function() {
    var fasta    = p.getArgs(0);
    var svbed    = p.getArgs(1);
    var coordbed = p.getArgs(2);
    if (!fasta || !svbed || !coordbed) throw new Error('requires three arguments.');
    if (! pa.existsSync(fasta)) throw new Error(fasta + ' : No such file.');
    if (! pa.existsSync(svbed)) throw new Error(svbed + ' : No such file.');
    if (! pa.existsSync(coordbed)) throw new Error(svbed + ' : No such file.');

    $s.shortcut('fasta', fasta);
    $s.shortcut('bedStream', new LS(svbed, {trim: true}));

    var sort = spawn('sortBed'); 
    fs.createReadStream(coordbed).pipe(sort.stdin);
    $s.shortcut('coordStream', new LS(sort.stdout, {trim: true}));
  })
  .catches(function(e) {
    console.ered('[Error] : ' +  e.message);
    showUsage();
    this.terminate();
  });

  $j('json', function() {
    var jsonfile = p.getOptions('json');
    return JSON.parse(fs.readFileSync(jsonfile).toString());
  })
  .failSafe(null)
  .post(function(v) {
    $s.shortcut('json', v);
  });

  $j(function() {
    $s.run();
  }).afterAbove();

  $j.run();
}

function getNewPos(rname, pos, regions, offset, diff) {
  var results = [];
  var svlist = regions[rname].SV;

  offset || (offset = 0);
  diff   || (diff= 0);
  var idx = svlist.bsearch([pos]);
  for (var i=offset; i<=idx; i++) {
    var svinfo  = svlist.get(i);
    var svstart = Number(svinfo[0]);
    var svend   = Number(svinfo[1]);
    var type    = svinfo[2];
    var extra   = (type == 'DUP') ? Number(svinfo[3]) : svinfo[3];

    if (svend >= pos) {
      break;
    }
    console.eyellow(type, svstart, svend, (type == 'INS') ? extra.length : extra);

    switch (type) {
      case 'DEL':
        diff -= svend - svstart + 1;
        break;
      case 'DUP':
        diff += (svend - svstart + 1) * (extra - 1);
        break;
      case 'INS':
        diff += extra.length;
        break;
      case 'INV':
        // no diff change
        break;
    }
  }
  return {pos: pos + diff, offset: i, diff: diff};
}

function getNewCoordinate(rname, start, end, regions, offset, diff) {
  var results = [];
  var svlist = regions[rname].SV;

  console.ecyan('--------------------------------------------------------------');
  console.ecyan('Coordinates:',start, end);
  var posinfo = getNewPos(rname, start, regions, offset, diff);
  offset = posinfo.offset;
  diff   = posinfo.diff;
  console.epurple('diff:', diff, 'offset', offset);

  var n = offset, sv = svlist.get(n), part = 0, nextype = null;
  var A = start, B = end, D = diff;
  while (sv && sv[0] <= end) {

    /* a: sv start, b: sv end, A: coord start, B: coord end **/
    var a = sv[0], b = sv[1];
    var abAB = [a, b, A, B], AB = ['A', 'B'], ab = ['a', 'b'];

    var ar = abAB.sort(function(i, j) {
      return (Number(i) > Number(j)) ? 1 : -1;
    })
    .map(function(v) { return (v == A || v == B) ? cl.cyan(v) : v });

    /** order of four points. aAbB means sv start => coord start => sv end => coord end **/
    var str = (a == B) ? 'AaBb'
             :(A == b) ? 'aAbB'
      : abAB.map(function(v) {
        return (v == A || v == B) 
          ? (AB.length) ? AB.shift() : ab.shift()
          : (ab.length) ? ab.shift() : AB.shift();
      }).join('');

    var type  = sv[2];
    var extra = sv[3];
    var trans = sv[4];

    switch (str) {
      case 'aAbB':
      /***
       *   A ---- B
       * a --- b
       **/
        switch (type) {
          case 'DEL':
            if (trans) {
              console.eblue('-----------TRA-------------');
              var trname = trans[0], tstart = trans[1];
              var D2 = getNewPos(trname, tstart, regions).pos - a;
              results.push({ rname: trname, start: A + D2, end: b + D2,  type: 'TRA', part: ++part, pstart: A, pend: b});
              console.eblue('---------------------------');
            }
            D -= b - a + 1;
            A = b + 1;
            nextype = 'LDEL';
            break;
          case 'INS':
            D += extra.length;
            break;
          case 'INV':
            results.push({ start: a + b - b + D, end: a + b - A + D, type: 'INV', part: ++part, pstart: A, pend: b});
            nextype = 'LINV';
            A = b + 1;
          case 'DUP':
            var len = b - a + 1;
            part++;
            for (var i=0; i<extra-1; i++) {
              results.push({ start: A + D, end: b + D, type: 'DUP', part: part, pstart: A, pend: b});
              D += len;
            }
            nextype = 'LDUP';
            // A = b + 1;
            break;
        }
        break;


      case 'aABb':
      /**
       *   A -- B
       * a ------- b
       **/
        switch (type) {
          case 'DEL':
            // totally deleted
            break;
          case 'INS':
            // impossible
            break;

          case 'INV':
            results.push({ start: a + b - B + D, end: a + b - A + D, type: 'INV', part: ++part, pstart: A, pend: B});
            break;
          case 'DUP':
            results.push({ start: A + D, end: B + D, part: ++part, pstart: A, pend: B });
            var len = b - a + 1;
            for (var i=0; i<extra-1; i++) {
              D += len;
              results.push({ start: A + D, end: B + D, type: 'DUP', part: part, pstart: A, pend: B});
            }
            break;
        }
        A = B = null;
        break;


      case 'AaBb':
      /**
       * A ----- B
       *    a ----- b
       **/
        switch (type) {
          case 'DEL':
            results.push({ start: A + D, end: a - 1 + D, type: 'RDEL', part: ++part, pstart: A, pend: a - 1});
            if (trans) {
              console.eblue('-----------TRA-------------');
              var trname = trans[0], tstart = trans[1];
              var D2 = getNewPos(trname, tstart, regions).pos - a;
              results.push({ rname: trname, start: a + D2, end: B + D2,  type: 'TRA', part: ++part, pstart: a, pend: B})
              console.eblue('---------------------------');
            }
            break;
          case 'INS':
            results.push({ start: A + D, end: B - 1 + D, type: 'RINS', part: ++part, pstart: A, pend: B - 1});
            results.push({ start: B + D, end: B + D + extra.length, type: 'INS', part: ++part, pstart: B, pend: B});
            break;

          case 'INV':
            results.push({ start: A + D, end: a + D - 1, type: 'RINV', part: ++part, pstart: A, pend: a - 1});
            results.push({ start: a + b - B + D, end: a + b - a + D, type: 'INS', part: ++part, pstart: a, pend: B});
            break;
          case 'DUP':
            results.push({ start: A + D, end: B + D, part: ++part, pstart: A, pend: B});
            var len = b - a + 1;
            part++;
            for (var i=0; i<extra-1; i++) {
              D += len;
              results.push({ start: a + D, end: B + D, type: 'DUP', part: part, pstart: a, pend: B});
            }
            break;
        }
        A = B = null;
        break;



      case 'AabB':
      /**
       * A --------- B
       *    a -- b
       **/
        switch (type) {
          case 'DEL':
            results.push({ start: A + D, end: a - 1 + D, type: 'RDEL', part: ++part, pstart: A, pend: a - 1});
            if (trans) {
              console.eblue('-----------TRA-------------');
              var trname = trans[0], tstart = trans[1];
              var D2 = getNewPos(trname, tstart, regions).pos - a;
              results.push({ rname: trname, start: a + D2, end: b + D2,  type: 'TRA', part: ++part, pstart: a, pend: b})
              console.eblue('---------------------------');
            }
            D -= b - a + 1;
            A = b + 1;
            nextype = 'LDEL';
            break;
          case 'INS':
            results.push({ start: A + D, end: a - 1 + D, type: 'RINS', part: ++part, pstart: A, pend: a - 1 });
            results.push({ start: a + D, end: a + D + extra.length, type: 'INS', part: ++part, pstart: a, pend: a});
            D += extra.length;
            A = b + 1;
            nextype = 'LINS';
            break;

          case 'INV':
            results.push({ start: A + D, end: a - 1 + D, type: 'RINV', part: ++part, pstart: A, pend: a - 1});
            results.push({ start: a + b - b + D, end: a + b - a + D, type: 'INV', part: ++part, pstart: a, pend: b});
            A = b + 1;
            nextype = 'LINV';
            break;
          case 'DUP':
            var len = b - a + 1;
            results.push({ start: A + D, end: b + D, type: 'RDUP', part: ++part, pstart: A, pend: b});
            part++
            for (var i=0; i<extra-1; i++) {
              D += len;
              results.push({ start: a + D, end: b + D, type: 'DUP', part: part, pstart: a, pend: b});
            }
            A = b + 1;
            nextype = 'LDUP';
            break;
        }
        break;

      default:
        console.ered("INVALID PATTERN", type + str);
        break;
    }

    ar.unshift(cl.red(type), cl.red(type == 'INS' ? extra.length : extra || ''), cl.red(n + ' OVERLAPPING'));
    console.ewhite.apply(console, ar);
    sv = svlist.get(++n);
  }

  if ( A != null & B != null) {
    results.push({ start: A + D, end: B + D, type: nextype, part: ++part, pstart: A, pend: B });
  }
  return {offset: offset, diff: diff, outputs: results};
}


if (process.argv[1] === __filename) { main() }
