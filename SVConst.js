const BASES = ['A', 'C', 'G', 'T'];
var dna = require('./lib/dna');
var a = {
  DEL   : 0,
  INS   : 1,
  INV   : 2,
  DUP   : 3,
  TRA   : 4,
  types : [],
  BASES : BASES
}

a.types[a.DEL] = 'DEL';
a.types[a.INS] = 'INS';
a.types[a.INV] = 'INV';
a.types[a.DUP] = 'DUP';
a.types[a.TRA] = 'TRA';

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
 * @param object  sv   : sv object (must have 'start', 'end' and 'fragment' keys)
 * @param string  str  : original DNA sequence part of which to be inserted
 * @param number  pos  : insertion position (char index)
 * @return string : inserted DNA sequence
 */
a.makeInsertion = function(sv, str, pos) {
  return str.slice(0, sv.start - pos) + sv.fragment + str.slice(sv.start - pos);
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
/**
 * SVConst.makeTandemDuplication
 * make tandem duplication to str with sv object
 * @param object  sv   : sv object (must have 'start' and 'end' keys)
 * @param string  str  : original DNA sequence part of which to be tandem duplicated
 * @param number  pos  : position (char index)
 * @return string : DNA sequence
 */
a.makeTandemDuplication = function(sv, str, pos) {
  return str.slice(0, sv.start - pos) + (function(){
    ret = "";
    var repeat_seq = str.slice(sv.start-pos, sv.end-pos);
    for (var i=0; i < sv.repeat_num; i++) {
      ret += repeat_seq;
    }
    return ret;
  })() + str.slice(sv.end - pos);
}

/**
 * SVConst.makeSNP
 * make inversion to str with sv object
 */
a.makeSNP = function(snp, str, pos) {
 return str.slice(0, snp.start - pos) + (function() {
   var c = str.charAt(snp.start - pos);
   var idx = BASES.indexOf(c.toUpperCase());
   return (idx >= 0) ? BASES[(idx + snp.to) % 4] : c;
 })() + str.slice(snp.start + 1 - pos);
}

module.exports = a;
