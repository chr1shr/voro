import math, cmath, copy, new, sys
import defaults, svg, glyphs, trans, curve

############################### class Fig

class Fig(trans.Delay):
  xmin = None
  xmax = None
  ymin = None
  ymax = None
  xlogbase = None
  ylogbase = None
  x = 10.
  y = 10.
  width = 80.
  height = 80.
  flipx = False
  flipy = False
  clip = False

  _varlist = ["xmin", "xmax", "ymin", "ymax", "xlogbase", "ylogbase", "x", "y", "width", "height", "flipx", "flipy", "clip"]

  def __init__(self, *args, **kwds):
    trans.Delay.__init__(self, *args, **kwds)

    for var in self._varlist:
      if var in kwds:
        self.__dict__[var] = kwds[var]
        del kwds[var]

    if len(kwds) > 0:
      raise TypeError, "Unrecognized keywords " + ", ".join(map(lambda word: "\"%s\"" % word, kwds.keys()))

    self.trans = None

  def __repr__(self):
    ren = "ren"
    if len(self.children) == 1: ren = ""
    clip = ""
    if self.clip: clip = " clip"
    return "<Fig (%d child%s) xmin=%s xmax=%s ymin=%s ymax=%s%s>" % (len(self.children), ren, self.xmin, self.xmax, self.ymin, self.ymax, clip)

  def transform(self, t):
    t = svg.cannonical_transformation(t)
    x1, y1 = t(self.x, self.y)
    x2, y2 = t(self.x + self.width, self.y + self.height)
    self.x, self.y = x1, y1
    self.width, self.height = x2 - x1, y2 - y1

  def bbox(self): return defaults.BBox(self.x, self.x + self.width, self.y, self.y + self.height)

  def svg(self):
    if self.xmin is not None and self.xmax is not None and \
       self.ymin is not None and self.ymax is not None:
      self.trans = trans.window(self.xmin, self.xmax, self.ymin, self.ymax,
                                x=self.x, y=self.y, width=self.width, height=self.height,
                                xlogbase=self.xlogbase, ylogbase=self.ylogbase,
                                minusInfinityX=(self.x - 10.*self.width), minusInfinityY=(self.y - 10.*self.height),
                                flipx=self.flipx, flipy=self.flipy)
    else:
      self.fit()

    self._svg = new.instance(svg.SVG)
    self._svg.__dict__["tag"] = "g"
    self._svg.__dict__["attrib"] = self.attrib
    self._svg.__dict__["_svg"] = self._svg

    self._svg.__dict__["children"] = []
    for child in self.children:
      self._svg.__dict__["children"].append(trans.transform(self.trans, child))

    if self.clip:
      clipPath = svg.SVG("clipPath", id=svg.randomid("clip-"))(svg.SVG("rect", self.x, self.y, self.width, self.height))
      self._svg["clip-path"] = "url(#%s)" % clipPath["id"]
      self._svg = svg.SVG("g", clipPath, self._svg)

  def __getstate__(self):
    mostdict = copy.copy(self.__dict__)
    if self.trans is not None:
      del mostdict["trans"]
      transcode = self.trans.func_code, self.trans.func_name
    else:
      transcode = None
    return (sys.version_info, defaults.version_info, mostdict, transcode)

  def __setstate__(self, state):
    self.__dict__ = state[2]
    self.__dict__["trans"] = []
    if state[3] is not None:
      code, name = state[3]
      context = globals()
      if "z" in code.co_names:
        context.update(cmath.__dict__)
      else:
        context.update(math.__dict__)
      self.__dict__["trans"] = new.function(code, context)
      self.__dict__["trans"].func_name = name
    else:
      self.__dict__["trans"] = None

  def __deepcopy__(self, memo={}):
    mostdict = copy.copy(self.__dict__)
    del mostdict["trans"]
    if "repr" in mostdict: del mostdict["repr"]
    output = new.instance(self.__class__)
    output.__dict__ = copy.deepcopy(mostdict, memo)
    output.__dict__["trans"] = self.trans

    memo[id(self)] = output
    return output

  def fit(self):
    bbox = defaults.BBox(None, None, None, None)
    for child in self.children:
      bbox += child.bbox()

    if self.xmin is not None: bbox.xmin = self.xmin
    if self.xmax is not None: bbox.xmax = self.xmax
    if self.ymin is not None: bbox.ymin = self.ymin
    if self.ymax is not None: bbox.ymax = self.ymax

    self.trans = trans.window(bbox.xmin, bbox.xmax, bbox.ymin, bbox.ymax,
                              x=self.x, y=self.y, width=self.width, height=self.height,
                              xlogbase=self.xlogbase, ylogbase=self.ylogbase,
                              minusInfinityX=(self.x - 10.*self.width), minusInfinityY=(self.y - 10.*self.height),
                              flipx=self.flipx, flipy=self.flipy)
    
############################### class Canvas

class Canvas(Fig):
  x = 0.
  y = 0.

  def __init__(self, width, height, *args, **kwds):
    Fig.__init__(self, *args, **kwds)
    self.width = width
    self.height = height
    self._tag = "svg"

  def __repr__(self):
    ren = "ren"
    if len(self.children) == 1: ren = ""
    return "<Canvas (%d child%s) width=%g height=%g xmin=%s xmax=%s ymin=%s ymax=%s>" % (len(self.children), ren, self.width, self.height, self.xmin, self.xmax, self.ymin, self.ymax)

  def svg(self):
    Fig.svg(self)
    output = self._svg
    self._svg = svg.SVG("svg", self.width, self.height, [0, 0, self.width, self.height])
    self._svg.__dict__["children"] = output.children
    self._svg.__dict__["attrib"].update(output.attrib)
