var EventEmitter = require('events').EventEmitter;
function GenomeStream(op) {
  op = op || {};
  this.nextream = op.nextream || process.stdout;
  this.interval = op.interval || 50;
  this.prelen = op.prelen || 0;
  this.delimiter = op.delimiter || '\n'; 

  this.pos = 0;
  this.remnant = '';

  this.on('data', function(data) {
    if (this.pos < this.prelen) {
      var len = Math.min(this.prelen-this.pos, data.length);
      this.nextream.write(data.slice(0, len));
      this.pos += len;
      if (data.length > len) {
        this.write(data.slice(len));
      }
    }

    else {
      var chunk = this.remnant + data.replace(/\n/g, '');
      this.remnant = '';

      while (chunk.length >= this.interval) {
        this.nextream.write(chunk.slice(0, this.interval)+ this.delimiter);
        this.pos += this.interval;
        chunk = chunk.slice(this.interval);
      }
      this.remnant = chunk;
    }
  });

  this.on('end', function() {
    if (this.remnant) {
      this.nextream.write(this.remnant + this.delimiter);
    }
    this.nextream.emit("end");
  });
}

GenomeStream.prototype = new EventEmitter();

GenomeStream.prototype.write = function(data) {
  this.emit('data', data);
}

GenomeStream.prototype.end = function() {
  this.emit('end');
}

module.exports = GenomeStream;
