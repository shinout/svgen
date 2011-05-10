const LS = require('linestream');
const AP = require('argparser');

function main() {
  var p = new AP().parse();

  var tsvfile = p.getArgs(0);
  
  if (! require('path').existsSync(tsvfile)) {
    process.stderr.write(tsvfile + ': No such file.');
    process.exit();
  }

  var ls = new LS(tsvfile, {trim : true});
  var svtypes = ['DEL', 'INS', 'INV'];

  ls.on('data', function(line) {
    if (!line || line.charAt(0) == '#') return;
    tsvdata = line.split('\t');
    if (tsvdata.length < 4) return;

    var type   = tsvdata[0];
    var start  = Number(tsvdata[1]);
    var length = Number(tsvdata[2]);
    var rname  = tsvdata[3];

    if (start <= 0 || length <= 0) return;

    var bpL, bpR;

    switch (type) {
    default: return;
    case 'DEL':
      bpR = start; 
      bpL = start + length; 
      break;
    case 'INS':
      bpL = start; 
      bpR = start; 
      break;
    case 'INV':
      process.stdout.write('R' + '\t' + start + '\t' + type + '\t' + rname + '\n');
      process.stdout.write('L' + '\t' + start + '\t' + type + '\t' + rname + '\n');

      process.stdout.write('R' + '\t' + (start + length) + '\t' + type + '\t' + rname + '\n');
      process.stdout.write('L' + '\t' + (start + length) + '\t' + type + '\t' + rname + '\n');
      return;
    }
    process.stdout.write('R' + '\t' + bpR + '\t' + type + '\t' + rname + '\n');
    process.stdout.write('L' + '\t' + bpL + '\t' + type + '\t' + rname + '\n');
  });

}

if (__filename == process.argv[1]) {
  main();
}

