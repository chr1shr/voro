import math, cmath, copy, re, sys, new
import defaults, svg, trans, pathdata, glyphs, _curve

############################### generic curve with marks (tick marks, arrows, etc)

class Curve(svg.SVG):
  attrib = {"stroke": "black", "fill": "none"}
  smooth = False
  marks = []
  random_sampling = True
  random_seed = 12345
  recursion_limit = 15
  linearity_limit = 0.05
  discontinuity_limit = 5.
  text_offsetx = 0.
  text_offsety = -2.5
  text_attrib = {}

  _varlist = ["attrib", "smooth", "marks", "random_sampling", "random_seed", "recursion_limit", "linearity_limit", "discontinuity_limit", "text_offsetx", "text_offsety", "text_attrib"]

  def __init__(self, expr, low, high, **kwds):
    self.__dict__["tag"] = None
    self.__dict__["children"] = []
    self.__dict__["_svg"] = []
    self.__dict__["f"] = svg.cannonical_parametric(expr)
    self.__dict__["low"] = low
    self.__dict__["high"] = high
    
    for var in self._varlist:
      if not callable(eval("self.%s" % var)) and var[:1] != "__" and var[-1:] != "__":
        if var in kwds:
          self.__dict__[var] = kwds[var]
          del kwds[var]
        else:
          self.__dict__[var] = copy.deepcopy(eval("self.%s" % var))

    # needed to set arrow color
    if "stroke" in kwds:
      self.attrib["stroke"] = kwds["stroke"]
      del kwds["stroke"]
    
    if "farrow" in kwds:
      if kwds["farrow"] is True:
        self.add(self.high, glyphs.farrowhead)
      else:
        self.add(self.high, kwds["farrow"])
      self.marks[-1][1]["fill"] = self["stroke"]
      del kwds["farrow"]
    
    if "barrow" in kwds:
      if kwds["barrow"] is True:
        self.add(self.low, glyphs.barrowhead)
      else:
        self.add(self.low, kwds["barrow"])
      self.marks[-1][1]["fill"] = self["stroke"]
      del kwds["barrow"]

    self.clean_arrows()

    # now the rest of the attributes (other than stroke)
    self.__dict__["attrib"].update(kwds)
    self.__dict__["trans"] = []

  def __call__(self, t, transformed=True):
    x, y = self.f(t)
    if transformed:
      for trans in self.trans:
        x, y = trans(x, y)
    return x, y

  def angle(self, t, transformed=True):
    x, y = self(t, transformed)

    tprime = t + trans.epsilon * abs(self.high - self.low)
    xprime, yprime = self(tprime, transformed)

    delx, dely = xprime - x, yprime - y
    return math.atan2(dely, delx)

  ### pickleability and access issues
  def __getattr__(self, name): return self.__dict__[name]
  def __setattr__(self, name, value): self.__dict__[name] = value

  def __getstate__(self):
    mostdict = copy.copy(self.__dict__)
    del mostdict["f"]
    del mostdict["trans"]
    transcode = map(lambda f: (f.func_code, f.func_name), self.trans)
    fcode = self.f.func_code, self.f.func_name
    return (sys.version_info, defaults.version_info, mostdict, transcode, fcode)

  def __setstate__(self, state):
    self.__dict__ = state[2]
    self.__dict__["trans"] = []

    for code, name in state[3]:
      context = globals()
      if "z" in code.co_names:
        context.update(cmath.__dict__)
      else:
        context.update(math.__dict__)
      f = new.function(code, context)
      f.func_name = name
      self.__dict__["trans"].append(f)

    context = globals()
    if "z" in state[4][0].co_names:
      context.update(cmath.__dict__)
    else:
      context.update(math.__dict__)
    self.__dict__["f"] = new.function(state[4][0], context)
    self.__dict__["f"].func_name = state[4][1]

  def __deepcopy__(self, memo={}):
    mostdict = copy.copy(self.__dict__)
    del mostdict["trans"]
    del mostdict["f"]
    if "repr" in mostdict: del mostdict["repr"]
    output = new.instance(self.__class__)
    output.__dict__ = copy.deepcopy(mostdict, memo)
    output.__dict__["trans"] = copy.copy(self.trans)
    output.__dict__["f"] = self.f

    memo[id(self)] = output
    return output

  ### presentation
  def __repr__(self):
    marks = ""
    if len(self.marks) == 1: marks = " (1 mark)"
    elif len(self.marks) > 1: marks = " (%d marks)" % len(self.marks)

    trans = ""
    if len(self.trans) > 0: trans = " (%d trans)" % len(self.trans)

    attrib = ""
    for var in "stroke", "fill":
      if var in self.attrib:
        attrib += " %s=%s" % (var, repr(self.attrib[var]))

    return "<%s %s from %g to %g%s%s%s>" % (self.__class__.__name__, self.f, self.low, self.high, marks, trans, attrib)

  ### transformation is like Delay
  def transform(self, t): self.trans.append(svg.cannonical_transformation(t))

  def bbox(self): return pathdata.bbox(self.d())

  ### construct the SVG path
  def svg(self):
    obj = new.instance(svg.SVG)
    obj.__dict__["tag"] = "path"
    obj.__dict__["attrib"] = self.attrib
    obj.__dict__["children"] = self.children
    obj.__dict__["_svg"] = obj
    obj.attrib["d"] = self.d()

    if self.marks == []:
      output = obj
    else:
      output = svg.SVG("g", obj)

      lowX, lowY = self(self.low)
      highX, highY = self(self.high)

      for item in self.marks:
        if isinstance(item, (int, long, float)):
          t, mark = item, glyphs.tick
        else:
          t, mark = item # marks should be (pos, mark) pairs or just pos

        X, Y = self(t)
        if self.low <= t <= self.high or \
           math.sqrt((X - lowX)**2 + (Y - lowY)**2) < trans.epsilon or \
           math.sqrt((X - highX)**2 + (Y - highY)**2) < trans.epsilon:

          angle = self.angle(t)

          if isinstance(mark, basestring):
            mark = self._render_text(X, Y, angle, mark)

          else:
            mark = trans.transform(lambda x, y: (X + math.cos(angle)*x - math.sin(angle)*y,
                                                 Y + math.sin(angle)*x + math.cos(angle)*y), mark)
          output.append(mark)

    self._svg = output

  def d(self):
    data = _curve.curve(self.f, self.trans, self.low, self.high,
                        self.random_sampling, self.random_seed, self.recursion_limit, self.linearity_limit, self.discontinuity_limit)

    segments = []
    last_d = None
    for d in data:
      if d is not None:
        if last_d is None: segments.append([])
        segments[-1].append(d)
      last_d = d

    output = []
    if self.smooth:
      for seg in segments:
        output.extend(pathdata.smooth(*seg))
    else:
      for seg in segments:
        output.extend(pathdata.poly(*seg))

    return output

  def _render_text(self, X, Y, angle, text):
    text_attrib = {"transform": "translate(%g, %g) rotate(%g)" %
                   (X + self.text_offsetx*math.cos(angle) - self.text_offsety*math.sin(angle),
                    Y + self.text_offsetx*math.sin(angle) + self.text_offsety*math.cos(angle), 180./math.pi*angle),
                   "text-anchor": "middle"}
    text_attrib.update(self.text_attrib)
    return svg.SVG("text", 0., 0., **text_attrib)(text)

  ### lots of functions for adding/removing marks
  def _matches(self, matching, mark):
    if matching is None: return True

    try:
      if isinstance(mark, matching): return True
    except TypeError: pass
    
    if isinstance(mark, svg.SVG) and isinstance(matching, basestring):
      if re.search(matching, mark.tag): return True
      if "repr" in mark.__dict__ and re.search(matching, mark.repr): return True

    if isinstance(matching, basestring) and isinstance(mark, basestring):
      if re.search(matching, mark): return True

    return matching == mark

  def wipe(self, low=None, high=None, matching=None):
    if low is None: low = self.low
    if high is None: high = self.high

    newmarks = []
    for item in self.marks:
      if isinstance(item, (int, long, float)):
        if self._matches(matching, item):
          if not low <= item <= high: newmarks.append(item)
        else:
          newmarks.append(item)

      else:
        pos, mark = item # marks should be (pos, mark) pairs or just pos
        if self._matches(matching, mark):
          if not low <= pos <= high: newmarks.append(item)
        else:
          newmarks.append(item)

    self.marks = newmarks

  def keep(self, low=None, high=None, matching=None):
    if low is None: low = self.low
    if high is None: high = self.high

    newmarks = []
    for item in self.marks:
      if isinstance(item, (int, long, float)):
        if self._matches(matching, item):
          if low <= item <= high: newmarks.append(item)

      else:
        pos, mark = item # marks should be (pos, mark) pairs or just pos
        if self._matches(matching, mark):
          if low <= pos <= high: newmarks.append(item)

    self.marks = newmarks

  def drop(self, t, tolerance=None, matching=None):
    if tolerance is None: tolerance = trans.epsilon * abs(self.high - self.low)
    self.wipe(t - tolerance, t + tolerance, matching=matching)

  def add(self, t, mark, angle=0., dx=0., dy=0.):
    if not isinstance(mark, basestring):
      mark = trans.transform(lambda x, y: (dx + math.cos(angle)*x - math.sin(angle)*y,
                                           dy + math.sin(angle)*x + math.cos(angle)*y), mark)
    self.marks.append((t, mark))

  def tick(self, t, mark=None):
    if mark is None:
      self.add(t, glyphs.tick)
      self.marks[-1][1]["stroke"] = self["stroke"]
    elif isinstance(mark, basestring):
      self.add(t, glyphs.tick)
      self.marks[-1][1]["stroke"] = self["stroke"]
      self.add(t, mark)
    else:
      self.add(t, mark)
      self.marks[-1][1]["stroke"] = self["stroke"]

  def minitick(self, t): self.tick(t, glyphs.minitick)

  def _markorder(self, a, b):
    if isinstance(a, (int, long, float)):
      posa, marka = a, None
    else:
      posa, marka = a # marks should be (pos, mark) pairs or just pos
    if isinstance(b, (int, long, float)):
      posb, markb = b, None
    else:
      posb, markb = b # marks should be (pos, mark) pairs or just pos

    if marka is None and markb is not None: return 1
    if marka is not None and markb is None: return -1
    return cmp(posa, posb)

  def sort(self, order=None):
    if order is None: order = lambda a, b: self._markorder(a, b)
    self.marks.sort(order)

  def closest(self, t, tolerance=None, matching=None):
    if tolerance is None: tolerance = trans.epsilon * abs(self.high - self.low)

    candidates = []
    for item in self.marks:
      if isinstance(item, (int, long, float)):
        if self._matches(matching, item) and abs(t - item) < tolerance:
          candidates.append(item)
      else:
        pos, mark = item
        if self._matches(matching, mark) and abs(t - pos) < tolerance:
          candidates.append(item)

    def closecmp(a, b):
      if isinstance(a, (int, long, float)):
        posa, marka = a, None
      else:
        posa, marka = a # marks should be (pos, mark) pairs or just pos
      if isinstance(b, (int, long, float)):
        posb, markb = b, None
      else:
        posb, markb = b # marks should be (pos, mark) pairs or just pos
      return cmp(abs(posa - t), abs(posb - t))

    candidates.sort(closecmp)
    return candidates

  def clean_arrows(self):
    end1, end2 = None, None
    for item in self.marks:
      if not isinstance(item, (int, long, float)):
        pos, mark = item # marks should be (pos, mark) pairs or just pos
        if mark not in (glyphs.tick, glyphs.minitick, glyphs.frtick, glyphs.frminitick) and \
           not isinstance(mark, basestring) and not (isinstance(mark, svg.SVG) and self.tag == "text"):
          if pos == self.low: end1 = item
          if pos == self.high: end2 = item

    newmarks = []
    for item in self.marks:
      if isinstance(item, (int, long, float)):
        if (item != self.low and item != self.high) or \
           (item == self.low and end1 is None) or \
           (item == self.high and end2 is None):
          newmarks.append(item)

      else:
        pos, mark = item
        if (mark not in (glyphs.tick, glyphs.minitick, glyphs.frtick, glyphs.frminitick)) or \
           (pos != self.low and pos != self.high) or \
           (pos == self.low and end1 is None) or \
           (pos == self.high and end2 is None):
          newmarks.append(item)

    self.marks = newmarks
        
  ### act like a list
  def append(self, other): self.marks.append(other)
  def prepend(self, other): self.marks[0:0] = [other]
  def insert(self, i, other): self.marks.insert(i, other)
  def remove(self, other): self.marks.remove(other)
  def __len__(self): return len(self.marks)

  def extend(self, other):
    if isinstance(other, SVG):
      self.marks.extend(other.children)
    elif isinstance(other, basestring):
      self.marks.append(other)
    else:
      self.marks.extend(other)

  def __add__(self, other):
    output = copy.deepcopy(self)
    output += other
    return output

  def __iadd__(self, other):
    self.marks.append(other)
    return self

  def __mul__(self, other):
    output = copy.deepcopy(self)
    output *= other
    return output

  def __rmul__(self, other):
    return self * other

  def __imul__(self, other):
    self.marks *= other
    return self

  def count(self, *args, **kwds): return self.marks.count(*args, **kwds)
  def index(self, *args, **kwds): return self.marks.index(*args, **kwds)
  def pop(self, *args, **kwds): return self.marks.pop(*args, **kwds)
  def reverse(self, *args, **kwds): return self.marks.reverse(*args, **kwds)

############################### plot axes

class XAxis(Curve):
  text_offsety = 2.5 + 3.
  xlogbase = None

  _varlist = Curve._varlist + ["xlogbase"]

  def _reassign_marks(self):
    if self.xlogbase is not None:
      return logticks(self.low, self.high)
    else:
      return ticks(self.low, self.high)

  def _reassign_f(self):
    output = eval("lambda t: (t, %s)" % repr(self.y))
    output.func_name = "x-value"
    return output
  
  def __init__(self, low, high, y, **kwds):
    self.__dict__["low"] = low
    self.__dict__["high"] = high
    self.y = y

    if "marks" not in kwds:
      kwds["marks"] = self._reassign_marks()

    Curve.__init__(self, self._reassign_f(), low, high, **kwds)

  def _render_text(self, X, Y, angle, text):
    text_attrib = {"transform": "translate(%g, %g) rotate(%g)" %
                   (X + self.text_offsetx*math.cos(angle) - self.text_offsety*math.sin(angle),
                    Y + self.text_offsetx*math.sin(angle) + self.text_offsety*math.cos(angle), 180./math.pi*angle),
                   "text-anchor": "middle"}
    text_attrib.update(self.text_attrib)
    return svg.SVG("text", 0., 0., **text_attrib)(text)

  def __setattr__(self, name, value):
    self.__dict__[name] = value
    if name == "xlogbase":
      self.marks = self._reassign_marks()
    if name in ("xlogbase", "low", "high"):
      self.f = self._reassign_f()

class YAxis(Curve):
  text_offsetx = -2.5
  text_offsety = 1.5 # when dominant-baseline is implemented everywhere, this hack will no longer be necessary
  ylogbase = None

  _varlist = Curve._varlist + ["ylogbase"]

  def _reassign_marks(self):
    if self.ylogbase is not None:
      return logticks(self.low, self.high)
    else:
      return ticks(self.low, self.high)

  def _reassign_f(self):
    output = eval("lambda t: (t, %s)" % repr(self.y))
    output.func_name = "y-value"
    return output
  
  def __init__(self, low, high, y, **kwds):
    self.__dict__["low"] = low
    self.__dict__["high"] = high
    self.y = y

    if "marks" not in kwds:
      kwds["marks"] = self._reassign_marks()

    Curve.__init__(self, self._reassign_f(), low, high, **kwds)

  def _render_text(self, X, Y, angle, text):
    angle += math.pi/2.
    text_attrib = {"transform": "translate(%g, %g) rotate(%g)" %
                   (X + self.text_offsetx*math.cos(angle) - self.text_offsety*math.sin(angle),
                    Y + self.text_offsetx*math.sin(angle) + self.text_offsety*math.cos(angle), 180./math.pi*angle),
                   "text-anchor": "end"}
    text_attrib.update(self.text_attrib)
    return svg.SVG("text", 0., 0., **text_attrib)(text)

  def __setattr__(self, name, value):
    self.__dict__[name] = value
    if name == "ylogbase":
      self.marks = self._reassign_marks()
    if name in ("ylogbase", "low", "high"):
      self.f = self._reassign_f()

############################### functions for making ticks

def format_number(x, format="%g", scale=1.):
  eps = trans.epsilon * abs(scale)
  if abs(x) < eps: return "0"
  return format % x

def unicode_number(x, scale=1.):
  """Converts numbers to a Unicode string, taking advantage of special
Unicode characters to make nice minus signs and scientific notation."""
  output = format_number(x, u"%g", scale)

  if output[0] == u"-":
    output = u"\u2013" + output[1:]

  index = output.find(u"e")
  if index != -1:
    uniout = unicode(output[:index]) + u"\u00d710"
    saw_nonzero = False
    for n in output[index+1:]:
      if n == u"+": pass # uniout += u"\u207a"
      elif n == u"-": uniout += u"\u207b"
      elif n == u"0":
        if saw_nonzero: uniout += u"\u2070"
      elif n == u"1":
        saw_nonzero = True
        uniout += u"\u00b9"
      elif n == u"2":
        saw_nonzero = True
        uniout += u"\u00b2"
      elif n == u"3":
        saw_nonzero = True
        uniout += u"\u00b3"
      elif u"4" <= n <= u"9":
        saw_nonzero = True
        if saw_nonzero: uniout += eval("u\"\\u%x\"" % (0x2070 + ord(n) - ord(u"0")))
      else: uniout += n

    if uniout[:2] == u"1\u00d7": uniout = uniout[2:]
    return uniout

  return output

def ticks(low, high, maximum=None, exactly=None, format=unicode_number):
  if exactly is not None:
    output = []
    t = low
    for i in xrange(exactly):
      output.append((t, glyphs.tick))
      output.append((t, format(t, scale=abs(high - low))))
      t += (high - low)/(exactly - 1.)
    return output

  if maximum is None: maximum = 10

  counter = 0
  granularity = 10**math.ceil(math.log10(max(abs(low), abs(high))))
  lowN = math.ceil(1.*low / granularity)
  highN = math.floor(1.*high / granularity)

  def subdivide(counter, granularity, low, high, lowN, highN):
    countermod3 = counter % 3
    if countermod3 == 0: granularity *= 0.5
    elif countermod3 == 1: granularity *= 0.4
    elif countermod3 == 2: granularity *= 0.5
    counter += 1
    lowN = math.ceil(1.*low / granularity)
    highN = math.floor(1.*high / granularity)
    return counter, granularity, low, high, lowN, highN
    
  while lowN > highN:
    counter, granularity, low, high, lowN, highN = \
             subdivide(counter, granularity, low, high, lowN, highN)

  last_granularity = granularity
  last_trial = None

  while True:
    trial = []
    for n in range(int(lowN), int(highN)+1):
      t = n * granularity
      trial.append(t)

    if len(trial) > maximum:
      if last_trial is None:
        v1, v2 = low, high
        return [(v1, format(v1, scale=abs(high - low))), (v2, format(v2, scale=abs(high - low)))]

      else:
        if counter % 3 == 2:
          counter, granularity, low, high, lowN, highN = \
                   subdivide(counter, granularity, low, high, lowN, highN)
        trial = []
        for n in range(int(lowN), int(highN)+1):
          t = n * granularity
          trial.append(t)

        output = []
        for t in last_trial:
          output.append((t, glyphs.tick))
          output.append((t, format(t, scale=abs(high - low))))
        for t in trial:
          if t not in last_trial:
            output.append((t, glyphs.minitick))
        return output

    last_granularity = granularity
    last_trial = trial

    counter, granularity, low, high, lowN, highN = \
             subdivide(counter, granularity, low, high, lowN, highN)

def logticks(low, high, base=10., maximum=None, format=unicode_number):
  if maximum is None: maximum = 10

  lowN = math.floor(math.log(low, base))
  highN = math.ceil(math.log(high, base))

  trial = []
  for n in range(int(lowN), int(highN)+1):
    trial.append(base**n)

  output = []

  # don't need every decade if the ticks cover too many
  for i in range(1, len(trial)):
    subtrial = trial[::i]
    if len(subtrial) <= maximum:
      for t in trial:
        output.append((t, glyphs.tick))
        if t in subtrial:
          output.append((t, format(t)))
      break

  if len(trial) <= 2:
    output2 = ticks(low, high, maximum=maximum, format=format)

    lowest = min(output2)
    for t, mark in output:
      if t < lowest: output2.append((t, mark))

    return output2

  for n in range(int(lowN), int(highN)+1):
    t = base**n
    for m in range(2, int(math.ceil(base))):
      output.append((m * t, glyphs.minitick))
        
  return output

