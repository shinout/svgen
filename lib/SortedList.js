function SortedList(arr, options) {
  if (arr instanceof Array) {
    return SortedList.build(arr, options);
  }
  this.arr = [];
  options = options || {};
  Object.keys(options).forEach(function(k) {
    this[k] = options[k];
  }, this);
}

SortedList.build = function(arr, options) {
  var list = new SortedList(null, options);
  arr.forEach(function(v) {
    list.insert(v);
  });
  return list;
};

SortedList.prototype.bsearch = function(val) {
  var mpos,
      spos = 0,
      epos = this.arr.length;
  while (epos - spos > 1) {
    mpos = Math.floor((spos + epos)/2);
    mval = this.arr[mpos];
    switch (this.compare(val, mval)) {
    case 1  :
    default :
      spos = mpos;
      break;
    case -1 :
      epos = mpos;
      break;
    case 0  :
      return mpos;
    }
  }
  return (this.arr[0] == null || spos == 0 && this.arr[0] != null && this.compare(this.arr[0], val) == 1) ? -1 : spos;
};

SortedList.prototype.get = function(pos) {
  return this.arr[pos];
};

SortedList.prototype.toArray = function(pos) {
  return this.arr.map(function(v) {return v;});
};

SortedList.prototype.size = function() {
  return this.arr.length;
};

SortedList.prototype.head = function() {
  return this.arr[0];
};

SortedList.prototype.tail = function() {
  return (this.arr.length == 0) ? null : this.arr[this.arr.length -1];
};

SortedList.prototype.insert = function(val) {
  var pos = this.bsearch(val);
  if (this.filter(val, pos)) {
    this.arr.splice(pos+1, 0, val);
    return true;
  }
  else {
    return false;
  }
};

SortedList.prototype.filter = function(val, pos) {
  return true;
};

SortedList.prototype.add = SortedList.prototype.insert;

SortedList.prototype.delete = function(pos) {
  this.arr.splice(pos, 1);
};
SortedList.prototype.remove = SortedList.prototype.delete;

SortedList.prototype.compare = function(a, b) {
  var c = a - b;
  return (c > 0) ? 1 : (c == 0)  ? 0 : -1;
};


function test() {
  // sample : ranges with no overlap
  var list = new SortedList([
    [152, 222], [33, 53], [48, 96], [928, 1743], [66, 67], [11, 12], [30, 32], [20,30]
  ],
  {
    filter: function(val, pos) {
      return (this.arr[pos]   == null || (this.arr[pos]   != null && this.arr[pos][1]  <  val[0])) 
        &&   (this.arr[pos+1] == null || (this.arr[pos+1] != null && val[1] < this.arr[pos+1][0]));
    },
    compare: function(a, b) {
      if (a == null) return -1;
      if (b == null) return  1;
      var c = a[0] - b[0];
      return (c > 0) ? 1 : (c == 0)  ? 0 : -1;
    }
  });


  console.log(list.arr);
  var cloned = list.toArray();
  cloned.push("hoge");
  console.log(cloned, list.arr);
}

if (process.argv[1] === __filename) { test(); }
module.exports = SortedList;
