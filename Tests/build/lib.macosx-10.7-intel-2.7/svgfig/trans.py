import math, cmath, sys, copy, new, warnings
import svg, defaults

epsilon = 1e-5

############################### copy-and-convert versions of main operations

clone = copy.deepcopy

def tonumber(obj):
  obj = copy.deepcopy(obj)
  obj.tonumber()
  return obj

def transform(trans, obj):
  if isinstance(trans, basestring):
    trans = svg.cannonical_transformation(trans)

  obj = copy.deepcopy(obj)
  if callable(trans):
    obj.transform(trans)
  else:
    for t in trans: obj.transform(t)
  return obj

def evaluate(obj):
  obj = copy.copy(obj) # start with a shallow copy
  if isinstance(obj, svg.SVG) and obj.tag is None:
    obj.svg()
    obj = obj._svg

  obj.__dict__["attrib"] = copy.deepcopy(obj.__dict__["attrib"])

  replacements = {}
  for i in xrange(len(obj.children)):
    obj.children[i] = evaluate(obj.children[i])

  return obj

############################### groups with special transformation properties

class Freeze(svg.SVG):
  def __init__(self, *args, **kwds):
    self.__dict__["tag"] = None
    self.__dict__["attrib"] = kwds
    self.__dict__["children"] = list(args)
    self.__dict__["_svg"] = None

  def __repr__(self):
    if len(self.children) == 1:
      return "<Freeze (1 child)>"
    else:
      return "<Freeze (%d children)>" % len(self.children)

  def transform(self, trans): pass

  def svg(self):
    self._svg = new.instance(svg.SVG)
    self._svg.__dict__["tag"] = "g"
    self._svg.__dict__["attrib"] = self.attrib
    self._svg.__dict__["children"] = self.children
    self._svg.__dict__["_svg"] = self._svg

class Delay(svg.SVG):
  def __init__(self, *args, **kwds):
    self.__dict__["tag"] = None
    self.__dict__["attrib"] = kwds
    self.__dict__["children"] = list(args)
    self.__dict__["_svg"] = None
    self.__dict__["trans"] = []

  def __repr__(self):
    if len(self.children) == 1:
      return "<Delay (1 child) (%d trans)>" % len(self.trans)
    else:
      return "<Delay (%d children) (%d trans)>" % (len(self.children), len(self.trans))

  def transform(self, trans):
    self.trans.append(svg.cannonical_transformation(trans))

  def bbox(self):
    self.svg()
    return self._svg.bbox()

  def svg(self):
    self._svg = new.instance(svg.SVG)
    self._svg.__dict__["tag"] = "g"
    self._svg.__dict__["attrib"] = self.attrib
    self._svg.__dict__["_svg"] = self._svg

    self._svg.__dict__["children"] = []
    for child in self.children:
      self._svg.__dict__["children"].append(transform(self.trans, child))

  def __getstate__(self):
    mostdict = copy.copy(self.__dict__)
    del mostdict["trans"]
    transcode = map(lambda f: (f.func_code, f.func_name), self.trans)
    return (sys.version_info, defaults.version_info, mostdict, transcode)

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

  def __deepcopy__(self, memo={}):
    mostdict = copy.copy(self.__dict__)
    del mostdict["trans"]
    if "repr" in mostdict: del mostdict["repr"]
    output = new.instance(self.__class__)
    output.__dict__ = copy.deepcopy(mostdict, memo)
    output.__dict__["trans"] = copy.copy(self.trans)
    memo[id(self)] = output
    return output

class Pin(svg.SVG):
  def __init__(self, x, y, *args, **kwds):
    self.__dict__["x"] = x
    self.__dict__["y"] = y
    if "rotate" in kwds:
      self.__dict__["rotate"] = kwds["rotate"]
      del kwds["rotate"]
    else:
      self.__dict__["rotate"] = False

    self.__dict__["tag"] = None
    self.__dict__["attrib"] = kwds
    self.__dict__["children"] = list(args)
    self.__dict__["_svg"] = None

  def __repr__(self):
    rotate = ""
    if self.rotate: rotate = "and rotate "
    ren = "ren"
    if len(self.children) == 1: ren = ""

    return "<Pin %sat %g %g (%d child%s)>" % (rotate, self.x, self.y, len(self.children), ren)

  def transform(self, trans):
    trans = svg.cannonical_transformation(trans)

    oldx, oldy = self.x, self.y
    self.x, self.y = trans(self.x, self.y)

    if self.rotate:
      shiftx, shifty = trans(oldx + epsilon, oldy)
      angle = math.atan2(shifty, shiftx)
      trans = eval("lambda x, y: (%(newx)s + cos(%(angle)s)*(x - %(oldx)s) - sin(%(angle)s)*(y - %(oldy)s), %(newy)s + sin(%(angle)s)*(x - %(oldx)s) + cos(%(angle)s)*(y - %(oldy)s))" % {"newx": repr(self.x), "newy": repr(self.y), "oldx": repr(oldx), "oldy": repr(oldy), "angle": repr(angle)}, math.__dict__)

    else:
      trans = eval("lambda x, y: (x + %(newx)s - %(oldx)s, y + %(newy)s - %(oldy)s)" %
                   {"newx": repr(self.x), "newy": repr(self.y), "oldx": repr(oldx), "oldy": repr(oldy)})

    for child in self.children:
      if isinstance(child, svg.SVG): child.transform(trans)

  def svg(self):
    self._svg = new.instance(svg.SVG)
    self._svg.__dict__["tag"] = "g"
    self._svg.__dict__["attrib"] = self.attrib
    self._svg.__dict__["children"] = self.children
    self._svg.__dict__["_svg"] = self._svg

############################### operations on transformations

def transformation_angle(expr, x, y, scale=1.):
  func = svg.cannonical_transformation(expr)
  eps = epsilon
  if scale != 0.: eps *= scale

  xprime, yprime = func(x + eps, y)
  x, y = func(x, y)

  delx, dely = xprime - x, yprime - y
  return math.atan2(dely, delx)

def transformation_jacobian(expr, x, y, scale=1.):
  func = svg.cannonical_transformation(expr)
  eps = epsilon
  if scale != 0.: eps *= scale

  X0, Y0 = func(x, y)
  xhatx, xhaty = func(x + eps, y)
  yhatx, yhaty = func(x, y + eps)

  return (1.*(xhatx - X0)/eps, 1.*(xhaty - Y0)/eps), (1.*(yhatx - X0)/eps, 1.*(yhaty - Y0)/eps)

############################### standard transformations

def window(xmin, xmax, ymin, ymax, x=0, y=0, width=100, height=100, xlogbase=None, ylogbase=None, minusInfinityX=-1000, minusInfinityY=-1000, flipx=False, flipy=False):
  if flipx:
    ox1 = x + width
    ox2 = x
  else:
    ox1 = x
    ox2 = x + width
  if flipy:
    oy1 = y + height
    oy2 = y
  else:
    oy1 = y
    oy2 = y + height
  ix1 = xmin
  iy1 = ymin
  ix2 = xmax
  iy2 = ymax
  
  if xlogbase != None and (ix1 <= 0. or ix2 <= 0.): raise ValueError, "x range incompatible with log scaling: (%g, %g)" % (ix1, ix2)

  if ylogbase != None and (iy1 <= 0. or iy2 <= 0.): raise ValueError, "y range incompatible with log scaling: (%g, %g)" % (iy1, iy2)

  xlogstr, ylogstr = "", ""

  if xlogbase == None:
    xfunc = "%(ox1)s + 1.*(x - %(ix1)s)/(%(ix2)s - %(ix1)s) * (%(ox2)s - %(ox1)s)" % \
            {"ox1": repr(ox1), "ox2": repr(ox2), "ix1": repr(ix1), "ix2": repr(ix2)}
  else:
    xfunc = "x <= 0 and %(minusInfinityX)s or %(ox1)s + 1.*(log(x, %(logbase)s) - log(%(ix1)s, %(logbase)s))/(log(%(ix2)s, %(logbase)s) - log(%(ix1)s, %(logbase)s)) * (%(ox2)s - %(ox1)s)" % \
            {"ox1": repr(ox1), "ox2": repr(ox2), "ix1": repr(ix1), "ix2": repr(ix2), "minusInfinityX": repr(minusInfinityX), "logbase": xlogbase}
    xlogstr = " xlog=%g" % xlogbase

  if ylogbase == None:
    yfunc = "%(oy1)s + 1.*(y - %(iy1)s)/(%(iy2)s - %(iy1)s) * (%(oy2)s - %(oy1)s)" % \
            {"oy1": repr(oy1), "oy2": repr(oy2), "iy1": repr(iy1), "iy2": repr(iy2)}
  else:
    yfunc = "y <= 0 and %(minusInfinityY)s or %(oy1)s + 1.*(log(y, %(logbase)s) - log(%(iy1)s, %(logbase)s))/(log(%(iy2)s, %(logbase)s) - log(%(iy1)s, %(logbase)s)) * (%(oy2)s - %(oy1)s)" % \
            {"oy1": repr(oy1), "oy2": repr(oy2), "iy1": repr(iy1), "iy2": repr(iy2), "minusInfinityY": repr(minusInfinityY), "logbase": ylogbase}
    ylogstr = " ylog=%g" % ylogbase

  output = eval("lambda x,y: (%s, %s)" % (xfunc, yfunc), math.__dict__)
  output.func_name = "(%g, %g), (%g, %g) -> (%g, %g), (%g, %g)%s%s" % (ix1, ix2, iy1, iy2, ox1, ox2, oy1, oy2, xlogstr, ylogstr)
  return output

def rotation(angle, cx=0, cy=0):
  output = eval("lambda x,y: (%(cx)s + cos(%(angle)s)*(x - %(cx)s) - sin(%(angle)s)*(y - %(cy)s), %(cy)s + sin(%(angle)s)*(x - %(cx)s) + cos(%(angle)s)*(y - %(cy)s))" % {"cx": repr(cx), "cy": repr(cy), "angle": repr(angle)}, math.__dict__)
  output.func_name = "rotation %g around %g %g" % (angle, cx, cy)
  return output






