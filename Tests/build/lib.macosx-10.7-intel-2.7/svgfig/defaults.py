import math, re, os, platform, warnings
import pathdata

version = "SVGFig 2.0.0alpha2"
version_info = (2, 0, 0, "alpha2")

############################### default filenames 

class VersionWarning(UserWarning): pass
warnings.filterwarnings("default", category=VersionWarning)

if re.search("windows", platform.system(), re.I):
  try:
    import _winreg
    directory = _winreg.QueryValueEx(_winreg.OpenKey(_winreg.HKEY_CURRENT_USER,
                r"Software\Microsoft\Windows\Current Version\Explorer\Shell Folders"), "Desktop")[0]
  except:
    directory = os.path.expanduser("~") + os.sep + "Desktop"

def expand_fileName(fileName):
  if re.search("windows", platform.system(), re.I) and not os.path.isabs(fileName):
    fileName = defaults.directory + os.sep + fileName
  return fileName

############################### Defaults for each SVG element type

xml_header = """\
<?xml version="1.0" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
"""

##### a
##### altGlyph
##### altGlyphDef
##### altGlyphItem
##### animate
##### animateColor
##### animateMotion
##### animateTransform
##### circle
signature_circle = ["cx", "cy", "r", "stroke", "fill"]
require_circle = ["cx", "cy", "r"]
defaults_circle = {"stroke": "black", "fill": "none"}
def tonumber_circle(svg):
  svg.cx, svg.cy, svg.r = tonumber(svg.cx), tonumber(svg.cy), tonumber(svg.r)

def transform_circle(trans, svg):
  if isnumber(svg.cx) and isnumber(svg.cy):
    x1, y1 = trans(svg.cx, svg.cy)
    if isnumber(svg.r):
      x2, y2 = trans(svg.cx + svg.r, svg.cy)
      svg.r = math.sqrt((x1 - x2)**2 + (y1 - y2)**2)
    svg.cx, svg.cy = x1, y1

def bbox_circle(svg):
  if isnumber(svg.cx) and isnumber(svg.cy):
    if isnumber(svg.r):
      return BBox(svg.cx - svg.r, svg.cx + svg.r, svg.cy - svg.r, svg.cy + svg.r)
    else:
      return BBox(svg.cx, svg.cx, svg.cy, svg.cy)
  else:
    return BBox(None, None, None, None)

##### clipPath
##### color-profile
##### cursor
##### definition-src
##### defs
signature_defs = None

##### desc
##### ellipse
##### feBlend
##### feColorMatrix
##### feComponentTransfer
##### feComposite
##### feConvolveMatrix
##### feDiffuseLighting
##### feDisplacementMap
##### feDistantLight
##### feFlood
##### feFuncA
##### feFuncB
##### feFuncG
##### feFuncR
##### feGaussianBlur
##### feImage
##### feMerge
##### feMergeNode
##### feMorphology
##### feOffset
##### fePointLight
##### feSpecularLighting
##### feSpotLight
##### feTile
##### feTurbulence
##### filter
##### font
##### font-face
##### font-face-format
##### font-face-name
##### font-face-src
##### font-face-uri
##### foreignObject
##### g
signature_g = None

##### glyph
##### glyphRef
##### hkern
##### image
##### line
signature_line = ["x1", "y1", "x2", "y2", "stroke"]
require_line = ["x1", "y1", "x2", "y2"]
defaults_line = {"stroke": "black"}

def tonumber_line(svg):
  svg.x1, svg.y1, svg.x2, svg.y2 = tonumber(svg.x1), tonumber(svg.y1), tonumber(svg.x2), tonumber(svg.y2)

def transform_line(trans, svg):
  if isnumber(svg.x1) and isnumber(svg.y1):
    svg.x1, svg.y1 = trans(svg.x1, svg.y1)

  if isnumber(svg.x2) and isnumber(svg.y2):
    svg.x2, svg.y2 = trans(svg.x2, svg.y2)

def bbox_line(svg):
  isnumber1 = (isnumber(svg.x1) and isnumber(svg.y1))
  isnumber2 = (isnumber(svg.x2) and isnumber(svg.y2))

  if isnumber1 and isnumber2: return BBox(svg.x1, svg.x2, svg.y1, svg.y2)
  elif isnumber1 and not isnumber2: return BBox(svg.x1, svg.x1, svg.y1, svg.y1)
  elif not isnumber1 and isnumber2: return BBox(svg.x2, svg.x2, svg.y2, svg.y2)
  else: return BBox(None, None, None, None)

##### linearGradient
##### marker
signature_marker = None

##### mask
##### metadata
##### missing-glyph
##### mpath
##### path
signature_path = ["d", "stroke", "fill"]
require_path = []
defaults_path = {"d": [], "stroke": "black", "fill": "none"}

def tonumber_path(svg):
  svg.d = pathdata.parse(svg.d)

def transform_path(trans, svg):
  svg.d = pathdata.transform(trans, svg.d)

def bbox_path(svg):
  return pathdata.bbox(pathdata.parse(svg.d))
  
##### pattern
##### polygon
##### polyline
##### radialGradient
##### rect
signature_rect = ["x", "y", "width", "height", "stroke", "fill"]
require_rect = ["x", "y", "width", "height"]
defaults_rect = {"stroke": "black", "fill": "none"}

def transform_rect(trans, svg):
  if isnumber(svg.x) and isnumber(svg.y):
    if isnumber(svg.width) and isnumber(svg.height):
      x1, y1 = trans(svg.x, svg.y)
      x2, y2 = trans(svg.x + svg.width, svg.y + svg.height)
      svg.x, svg.y = x1, y1
      svg.width, svg.height = x2 - x1, y2 - y1
    else:
      svg.x, svg.y = trans(svg.x, svg.y)

def bbox_rect(svg):
  if isnumber(svg.x) and isnumber(svg.y):
    if isnumber(svg.width) and isnumber(svg.height):
      return BBox(svg.x, svg.x + svg.width, svg.y, svg.y + svg.height)
    else:
      return BBox(svg.x, svg.x, svg.y, svg.y)
  else:
    return BBox(None, None, None, None)

##### script
##### set
##### stop
##### style
##### svg
signature_svg = ["width", "height", "viewBox"]
require_svg = []
defaults_svg = {"width": 400, "height": 400, "viewBox": (0, 0, 100, 100),
                "style": {"stroke-width": "0.5pt", "font-size": "4px", "text-anchor": "middle"},
                "font-family": ["Helvetica", "Arial", "FreeSans", "Sans", "sans", "sans-serif"],
                "xmlns": "http://www.w3.org/2000/svg", "xmlns:xlink": "http://www.w3.org/1999/xlink", "version":"1.1",
                }

def tonumber_svg(svg):
  svg.width = tonumber(svg.width)
  svg.height = tonumber(svg.height)
  svg.viewBox = tonumberlist(svg.viewBox)
  svg["style"] = tostringmap(svg["style"])
  svg["font-family"] = tostringlist(svg["font-family"])

##### switch
##### symbol
signature_symbol = None

##### text
signature_text = ["x", "y", "stroke", "fill"]
require_text = ["x", "y"]
defaults_text = {"stroke": "none", "fill": "black"}

def tonumber_text(svg):
  svg.x = tonumber(svg.x)
  svg.y = tonumber(svg.y)

def transform_text(trans, svg):
  if isnumber(svg.x) and isnumber(svg.y):
    svg.x, svg.y = trans(svg.x, svg.y)

def bbox_text(svg):
  if isnumber(svg.x) and isnumber(svg.y):
    return BBox(svg.x, svg.x, svg.y, svg.y) # how to calculate text size???
  else:
    return BBox(None, None, None, None)

##### textPath
##### title
##### tref
##### tspan
signature_tspan = None

##### use
signature_use = ["x", "y", "xlink:href"]
require_use = ["x", "y", "xlink:href"]

def tonumber_use(svg):
  svg.x = tonumber(svg.x)
  svg.y = tonumber(svg.y)

def transform_use(trans, svg):
  if isnumber(svg.x) and isnumber(svg.y):
    svg.x, svg.y = trans(svg.x, svg.y)

def bbox_use(svg):
  if isnumber(svg.x) and isnumber(svg.y):
    return BBox(svg.x, svg.x, svg.y, svg.y)
  else:
    return BBox(None, None, None, None)

##### view
##### vkern

############################### utility functions for default actions

def tonumber(x):
  try:
    return float(x)
  except ValueError:
    return x

def tonumberlist(x):
  if isinstance(x, basestring):
    try:
      return tuple(map(float, re.split("[, \t]+", x)))
    except ValueError:
      return x
  return x

def tostringlist(x):
  if isinstance(x, basestring):
    return re.split("[, \t]+", x)
  return x

def tostringmap(x):
  if isinstance(x, basestring):
    try:
      return dict(map(lambda word: word.split(":"), re.split("[; \t]+", x)))
    except ValueError:
      return x
  else:
    return x

def isnumber(x): return isinstance(x, (int, long, float))

############################### BBox class

class BBox:
  def __init__(self, xmin, xmax, ymin, ymax):
    self.xmin, self.xmax, self.ymin, self.ymax = xmin, xmax, ymin, ymax

  def __repr__(self):
    return "<BBox xmin=%g xmax=%g ymin=%g ymax=%g>" % (self.xmin, self.xmax, self.ymin, self.ymax)

  def insert(self, x, y):
    if self.xmin == None or x < self.xmin: self.xmin = x
    if self.ymin == None or y < self.ymin: self.ymin = y
    if self.xmax == None or x > self.xmax: self.xmax = x
    if self.ymax == None or y > self.ymax: self.ymax = y

  def __add__(self, other):
    output = BBox(self.xmin, self.xmax, self.ymin, self.ymax)
    output += other
    return output

  def __iadd__(self, other):
    if self.xmin is None: self.xmin = other.xmin
    elif other.xmin is None: pass
    else: self.xmin = min(self.xmin, other.xmin)

    if self.xmax is None: self.xmax = other.xmax
    elif other.xmax is None: pass
    else: self.xmax = max(self.xmax, other.xmax)

    if self.ymin is None: self.ymin = other.ymin
    elif other.ymin is None: pass
    else: self.ymin = min(self.ymin, other.ymin)

    if self.ymax is None: self.ymax = other.ymax
    elif other.ymax is None: pass
    else: self.ymax = max(self.ymax, other.ymax)

    return self

  def __eq__(self, other):
    return self.xmin == other.xmin and self.xmax == other.xmax and self.ymin == other.ymin and self.ymax == other.ymax

  def __ne__(self, other): return not (self == other)
