if (typeof global != 'undefined') require('../lib/Junjo/test/load.test').load(global);
const FR          = require('../lib/FASTAReader/FASTAReader');


var f1 = new FR('chrM.fa');
var e1 = f1.result['chrM'].getEndPos();


var f2 = new FR('sv.fasta');
var e2 = f2.result['chrM'].getEndPos();

console.log(e1, e2);
