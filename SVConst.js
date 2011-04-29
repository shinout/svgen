var dna = require('./lib/dna');
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
  dna.complStrand(str.slice(sv.start-pos, sv.end-pos)).split('').reverse().join('') + 
  str.slice(sv.end - pos);
}

module.exports = a;
