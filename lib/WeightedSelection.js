function WeightedSelection(table, random) {
  // TODO validate first.
	
  this.names = Object.keys(table).filter(function(name) {
		return (table[name] > 0);
	});

  this.hists = this.names.map((function() {
    var total = 0;
    return function(v) {
      total += Number(table[v]);
      return total;
    };
  })());
  this.hists.unshift(0);

  if (typeof random == 'function') {
    this.fn = random;
  }
}

WeightedSelection.prototype.fn = function() {
  return Math.random();
};

WeightedSelection.prototype.total = function(val) {
  return this.hists[this.hists.length -1];
};

WeightedSelection.prototype.count = function(val) {
  return this.names.length;
};

WeightedSelection.prototype.random = function(val) {
  var stpos  = 0;
  var endpos = this.hists.length -1;

  // TODO filter NaN, over max
  var val = (val != null) ? Number(val) : Math.floor(this.fn() * this.hists[endpos]);

  var cenpos, cenval, sub;
  while ( endpos - stpos > 1) {
    cenpos = Math.floor((stpos + endpos) / 2);
    cenval = this.hists[cenpos];
    sub = val - cenval;
    if (sub > 0) {
      stpos = cenpos;
    }
    else if (sub < 0) {
      endpos = cenpos;
    }
    // exact match
    else {
      return this.names[cenpos];
    }
  }
  return this.names[stpos];
};


function test() {
  const chrs = {
    notselected1 : 0,
    chr1         : 200000,
    chr2         : 180000,
    chr3         : 170000,
    chr4         : 150000,
    chr5         : 140000,
    notselected2 : 0,
    chr6         : 120000
  };


  var sel = new WeightedSelection(chrs);
	console.log(sel);
  console.log(sel.random(0) == 'chr1');
  console.log(sel.random(1) == 'chr1');
  console.log(sel.random(199999) == 'chr1');
  console.log(sel.random(200000) == 'chr2');
  console.log(sel.random(200001) == 'chr2');
  console.log(sel.random(379999) == 'chr2');
  console.log(sel.random(380000) == 'chr3');
  console.log(sel.random(920000) == 'chr6');

  console.log(sel.random(1920000) == 'chr6'); // OVER
  console.log(sel.random(-1) == 'chr1'); // OVER

  var count = 0;
  const N = 100000;
  for (var i=0; i<N; i++) {
    if (sel.random() == 'chr1') {
      count++;
    }
  }
  console.log(count / N);
  console.log(chrs['chr1'] / sel.hists[sel.hists.length -1]);
}


if (process.argv[1] === __filename) { test(); }
module.exports = WeightedSelection;
