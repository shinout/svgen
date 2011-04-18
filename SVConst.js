var a = {
  DEL: 0,
  INS: 1,
  INV: 2,
  types: []
}

a.types[a.DEL] = 'DEL';
a.types[a.INS] = 'INS';
a.types[a.INV] = 'INV';

/**
 * SVConst.getRandomFlagment
 * get random DNA flagment
 * @param number  len  : flagment length
 * @return string : DNA flagment
 */
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

/**
 * SVConst.idx2pos
 * convert charcter index to DNA base position
 * @param number  idx     : character index
 * @param number  prelen  : header data length
 * @param number  linelen : one line length
 * @return number : DNA base position
 */
a.idx2pos = function(idx, prelen, linelen) {
  prelen = prelen || 0;
  linelen = linelen || 50;
  return Math.max(0, idx +1 - prelen - Math.floor((idx +1 - prelen)/(linelen + 1)));
}


/**
 * SVConst.pos2index
 * convert DNA base position to character index
 * @param number  pos     : DNA base position
 * @param number  prelen  : header data length
 * @param number  linelen : one line length
 * @return number : character index
 */
a.pos2index = function(pos, prelen, linelen) {
  prelen = prelen || 0;
  linelen = linelen || 50;
  return prelen + pos -1 + Math.floor( (pos -1)/linelen );
}

/**
 * SVConst.makeDeletion
 * make deletion from str with sv object
 * @param object  sv   : sv object (must have 'start' and 'end' keys)
 * @param string  str  : original DNA sequence part of which to be deleted
 * @param number  pos  : delete position (char index)
 * @return string : deleted DNA sequence
 */
a.makeDeletion = function(sv, str, pos) {
  return str.slice(0, sv.start - pos) + str.slice(sv.end - pos);
}

/**
 * SVConst.makeInsertion
 * make insertion to str with sv object
 * @param object  sv   : sv object (must have 'start', 'end' and 'flagment' keys)
 * @param string  str  : original DNA sequence part of which to be inserted
 * @param number  pos  : insertion position (char index)
 * @return string : inserted DNA sequence
 */
a.makeInsertion = function(sv, str, pos) {
  return str.slice(0, sv.start - pos) + sv.flagment + str.slice(sv.start - pos);
}

/**
 * SVConst.makeInversion
 * make inversion to str with sv object
 * @param object  sv   : sv object (must have 'start' and 'end' keys)
 * @param string  str  : original DNA sequence part of which to be inverted
 * @param number  pos  : inverted position (char index)
 * @return string : inverted DNA sequence
 */
a.makeInversion = function(sv, str, pos) {
  return str.slice(0, sv.start - pos) + 
  this.complStrand(str.slice(sv.start-pos, sv.end-pos)).split('').reverse().join('') + 
  str.slice(sv.end - pos);
}

/**
 * SVConst.complStrand
 * make a complementary strand of str
 * @param string  str  : DNA sequence
 * @return string : complementary DNA sequence
 */
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
