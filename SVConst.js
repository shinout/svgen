var a = {
  DEL: 0,
  INS: 1,
  INV: 2,
  types: []
}

a.types[a.DEL] = 'DEL';
a.types[a.INS] = 'INS';
a.types[a.INV] = 'INV';

/* static functions */
/* get random DNA flagment */
a.getRandomFlagment = function(len) {
  var flagment = '';
  for (var i=0; i<len; i++) {
    var p = Math.random();
    if (p > 0.75) {
      flagment += 'A';
    }
    else if (p > 0.5) {
      flagment += 'G';
    }
    else if (p > 0.25) {
      flagment += 'T';
    }
    else {
      flagment += 'C';
    }
  }
  return flagment;
}

a.idx2pos = function(idx, prelen, linelen) {
  prelen = prelen || 0;
  linelen = linelen || 50;
  return Math.max(0, idx +1 - prelen - Math.floor((idx +1 - prelen)/(linelen + 1)));
}


a.pos2index = function(pos, prelen, linelen) {
  prelen = prelen || 0;
  linelen = linelen || 50;
  return prelen + pos -1 + Math.floor( (pos -1)/linelen );
}

a.makeDeletion = function(sv, str, pos) {
  return str.slice(0, sv.start - pos) + str.slice(sv.end - pos);
}

a.makeInsertion = function(sv, str, pos) {
  return str.slice(0, sv.start - pos) + sv.flagment + str.slice(sv.start - pos);
}

a.makeInversion = function(sv, str, pos) {
  return str.slice(0, sv.start - pos) + 
  this.complStrand(str.slice(sv.start-pos, sv.end-pos)).split('').reverse().join('') + 
  str.slice(sv.end - pos);
}

a.complStrand = function(str) {
  var ret = [];
  var i = 0;
  str.split('').forEach(function(c) {
    switch (c) {
      case 'a':
        ret[i] = 't';
        break;
      case 'A':
        ret[i] = 'T';
        break;
      case 't':
        ret[i] = 'a';
        break;
      case 'T':
        ret[i] = 'A';
        break;
      case 'c':
        ret[i] = 'g';
        break;
      case 'C':
        ret[i] = 'G';
        break;
      case 'g':
        ret[i] = 'c';
        break;
      case 'G':
        ret[i] = 'C';
        break;
      default:
        ret[i] = c;
        break;
    }
    i++;
  });
  return ret.join('');
}





module.exports = a;
