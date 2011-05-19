/**
 * random number generator generator
 * @param (Number) seed 
 * @param (boolean) normalize : if true, range becomes 0 to 1 
 * @return function : random number generator
 *
 */
function XORShift(seed, normalize) {
  /* prepare seed */
  seed = (typeof seed == Number) ? seed : Math.floor(Math.random() * 10000000);
  var sq = Math.sqrt(Math.random() + seed);
  var num = sq.toString().replace('.', '').replace(/^0*/, '');
  var len = num.length;
  num = Number(num);

  var w = (num + new Date().getTime() >>> ((seed + len) % 8)) >>> 0,
      x = 123456789,
      y = 362436069,
      z = 521288629;

  /**
   * random number generator 
   * @return if normalize is false, integer ( 0 to 0x100000000 )
   *         if normalize is true , float ( 0 to 1 )
   */
  return function() {
    var t = x ^ (x << 11);
    x = y; y = z; z = w;
    w = ((w ^ (w >>> 19)) ^ (t ^ (t>>>8))) >>> 0;
    return (normalize) ? w/0x100000000 : w;
  }
}

module.exports = XORShift;
