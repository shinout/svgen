function normal_random(mean, dev, fn) {
  fn = (typeof fn == 'function') ? fn : function(v) {return Math.random(v)};
  var a = fn();
  var b = fn();
  with (Math) {
    return dev * sqrt(-2 * log(a)) * sin(2 * PI * b) + mean;
  }
}

module.exports = normal_random;
