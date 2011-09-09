const fs    = require('fs');
const AP    = require('argparser');
const LS    = require('linestream');
const Junjo = require('./lib/Junjo/Junjo');
const SVGen = require('./svgen');
const cl    = require('./lib/termcolor').define();
const pa    = require('path');
const spawn = require('child_process').spawn;

function main() {
  var p  = new AP().addOptions([]).addValueOptions([]).parse();
  var $j = new Junjo();

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
    return Junjo.multi(fasta, svbed, coordbed);
  })
  .catches(function(e) {
    console.ered('[Error] : ' +  e.message);
    showUsage();
    $j.terminate();
  });

  $j('json', function() {
    var jsonfile = p.getOptions('json');
    return JSON.parse(fs.readFileSync(jsonfile).toString());
  })
  .failSafe(null);

  $j('registerBeds', function(fasta, svbed, coordbed, json) {
    var svgen = new SVGen(fasta, {json: json, rnames: null});
    var bedstream = new LS(svbed, {trim: true});

    this.absorb(bedstream, 'data', function(line, result) {
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
    this.absorbError(bedstream);
  })
  .firstError('shift')
  .after('inputCheck', 'json');

  $j('sort', function(fasta, svbed, coordbed, svgen) {
    var sort = spawn('sortBed'); 
    fs.createReadStream(coordbed).pipe(sort.stdin);

    var coord = new LS(sort.stdout, {trim: true});

    sort.stderr.on('data', function(data) {
      console.log(data.toString());
    });
    this.absorbError(sort.stderr, 'data');

    this.absorb(coord, 'data', function(line, result) {
      if (!line || line.charAt(0) == '#') return;
      var svinfo = line.split('\t');
      if (svinfo.length < 3) return;

      var rname  = svinfo[0];
      var regions = svgen.regions[rname];
      var svpos = regions.SV;
      var snvs  = regions.SV;
      var start  = Number(svinfo[1]);
      var end    = svinfo[2];
      console.log(rname, start, end);
    });
  })
  .firstError('shift')
  .after('inputCheck', 'registerBeds')
  .next(function(out) {
    console.log(out);
  });

  $j.catchesAbove(function(e, args) {
    console.ered('[Error] : ' +  e.message + ' in label: ' + cl.yellow(this.label()));
    console.error(e.stack || e);
    $j.terminate();
  });

  $j.run();
}


if (process.argv[1] === __filename) { main() }
