import math, cmath, random, re, os, sys, copy, itertools, codecs, tempfile, new, types, copy_reg, warnings
import defaults

saved = [] # keep track of all fileNames saved for the user's convenience

############################### convenient functions for dealing with SVG

def randomid(prefix="", characters=10):
  return prefix + "".join(random.sample("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789", characters))

# This rgb function could be a lot better... something to think about...
def rgb(r, g, b, maximum=1.):
  return "#%02x%02x%02x" % (max(0, min(r*255./maximum, 255)), max(0, min(g*255./maximum, 255)), max(0, min(b*255./maximum, 255)))

############################### class SVG

def shortcut(tag): return eval("lambda *args, **kwds: SVG(\"%s\", *args, **kwds)" % tag)

class SVG:
  def _preprocess_attribname(self, name):
    name_colon = re.sub("__", ":", name)
    if name_colon != name: name = name_colon

    name_dash = re.sub("_", "-", name)
    if name_dash != name: name = name_dash

    return name

  def __init__(self, tag, *signature_attrib, **more_attrib):
    self.__dict__["tag"] = tag
    self.__dict__["attrib"] = dict(getattr(defaults, "defaults_%s" % tag, {}))
    self.__dict__["children"] = []
    self.__dict__["_svg"] = self

    signature = getattr(defaults, "signature_%s" % tag, None)

    # if there is no signature, inline arguments are interpreted as children
    if signature is None:
      self.children.extend(signature_attrib)

    else:
      if len(signature_attrib) > len(signature):
        raise TypeError, "Tag '%s' expects no more than %d signature attributes (saw %d)" % (tag, len(signature), len(signature_attrib))

      for name, value in zip(signature, signature_attrib):
        self.attrib[name] = value

    # preprocess more_attrib names
    for name, value in more_attrib.items():
      processed_name = self._preprocess_attribname(name)
      if processed_name != name:
        del more_attrib[name]
        more_attrib[processed_name] = value

    self.attrib.update(more_attrib)

    require = getattr(defaults, "require_%s" % tag, None)
    if require is not None:
      for name in require:
        if name not in self.attrib:
          raise TypeError, "Tag '%s' requires a '%s' attribute" % (tag, name)

  ### construct trees inline
  def __call__(self, *children):
    self.children.extend(children)
    return self

  ### recursively tonumber, transform, bbox, and svg
  def tonumber(self):
    if self.tag is not None:
      tonumber_tag = getattr(defaults, "tonumber_%s" % self.tag, None)
      if tonumber_tag is not None: tonumber_tag(self)
    
    for child in self.children:
      if isinstance(child, SVG): child.tonumber()

  def transform(self, t):
    t = cannonical_transformation(t)

    if self.tag is not None:
      tonumber_tag = getattr(defaults, "tonumber_%s" % self.tag, None)
      if tonumber_tag is not None: tonumber_tag(self)

      transform_tag = getattr(defaults, "transform_%s" % self.tag, None)
      if transform_tag is not None: transform_tag(t, self)

    for child in self.children:
      if isinstance(child, SVG): child.transform(t)

  def bbox(self):
    if self.tag is not None:
      tonumber_tag = getattr(defaults, "tonumber_%s" % self.tag, None)
      if tonumber_tag is not None: tonumber_tag(self)

      bbox_tag = getattr(defaults, "bbox_%s" % self.tag, None)
      if bbox_tag is not None:
        output = bbox_tag(self)
      else:
        output = defaults.BBox(None, None, None, None)

    for child in self.children:
      if isinstance(child, SVG): output += child.bbox()
    return output

  def svg(self): self._svg = self

  ### signature attributes are accessible as member data
  def __getattr__(self, name):
    if self.__dict__["tag"] is None: return self.__dict__[name]

    signature = getattr(defaults, "signature_%s" % self.__dict__["tag"], None)
    if signature is not None and name in signature:
      return self.attrib[name]
    else:
      raise AttributeError, "Tag '%s' has no signature attrib '%s'" % (self.tag, name)

  def __setattr__(self, name, value):
    if self.__dict__["tag"] is None or name == "repr" or name in self.__dict__:
      self.__dict__[name] = value

    else:
      signature = getattr(defaults, "signature_%s" % self.__dict__["tag"], None)
      if signature is not None and name in signature:
        self.attrib[name] = value
      else:
        raise AttributeError, "Tag '%s' has no signature attrib '%s'" % (self.tag, name)

  def __nonzero__(self): return True

  ### support access to deep children with tree indexes
  def _treeindex_descend(self, obj, treeindex):
    if isinstance(treeindex, (list, tuple)):
      for i in treeindex[:-1]:
        obj = obj[i]
      treeindex = treeindex[-1]
    return treeindex, obj

  def __getitem__(self, treeindex):
    treeindex, obj = self._treeindex_descend(self, treeindex)

    if isinstance(treeindex, (int, long, slice)): return obj.children[treeindex]
    elif isinstance(treeindex, basestring): return obj.attrib[treeindex]
    else:
      raise IndexError, "treeindex must be [#, #, ... #] or [#, #, ... \"str\"]"

  def __setitem__(self, treeindex, value):
    treeindex, obj = self._treeindex_descend(self, treeindex)

    if isinstance(treeindex, (int, long, slice)): obj.children[treeindex] = value
    elif isinstance(treeindex, basestring): obj.attrib[treeindex] = value
    else:
      raise IndexError, "treeindex must be [#, #, ... #] or [#, #, ... \"str\"]"

  def __delitem__(self, treeindex):
    treeindex, obj = self._treeindex_descend(self, treeindex)

    if isinstance(treeindex, (int, long, slice)): del obj.children[treeindex]
    elif isinstance(treeindex, basestring): del obj.attrib[treeindex]
    else:
      raise IndexError, "treeindex must be [#, #, ... #] or [#, #, ... \"str\"]"

  ################ nested class for walking the tree
  class _SVGDepthIterator:
    def __init__(self, svg, treeindex, depth_limit, attrib, attrib_first):
      self.current = svg
      self.treeindex = treeindex
      self.shown = False
      self.depth_limit = depth_limit
      self.attrib = attrib
      self.attrib_first = attrib_first

    def __iter__(self): return self

    def make_children_iterators(self):
      if getattr(self.current, "children", None) is not None:
        for i, s in enumerate(self.current.children):
          self.iterators.append(self.__class__(s, self.treeindex + (i,), self.depth_limit, self.attrib, self.attrib_first))

    def make_attrib_iterators(self):
      if getattr(self.current, "attrib", None) is not None:
        items = self.current.attrib.items()
        items.sort()
        for k, s in items:
          self.iterators.append(self.__class__(s, self.treeindex + (k,), self.depth_limit, self.attrib, self.attrib_first))

    def next(self):
      if not self.shown:
        self.shown = True
        if self.treeindex != ():
          return self.treeindex, self.current

      if self.depth_limit is not None and len(self.treeindex) >= self.depth_limit: raise StopIteration

      if "iterators" not in self.__dict__:
        self.iterators = []

        if self.attrib and self.attrib_first: self.make_attrib_iterators()
        self.make_children_iterators()
        if self.attrib and not self.attrib_first: self.make_attrib_iterators()

        self.iterators = itertools.chain(*self.iterators)

      return self.iterators.next()
  ################ end nested class

  ### walk the tree or show it (uses the nested class)
  def walk(self, depth_limit=None, attrib=False, attrib_first=False):
    return self._SVGDepthIterator(self, (), depth_limit, attrib, attrib_first)

  def tree(self, depth_limit=None, attrib=False, attrib_first=False, index_width=20, showtop=True, asstring=False):
    if showtop:
      output = [("%s %s" % (("%%-%ds" % index_width) % repr(None), repr(self)))]
    else:
      output = []

    for treeindex, element in self.walk(depth_limit, attrib, attrib_first):
      if isinstance(element, basestring):
        if len(element) > 13:
          repr_element = "'%s...'" % element[0:10]
        else:
          repr_element = "'%s'" % element
      else:
        repr_element = repr(element)

      output.append(("%s %s%s" % (("%%-%ds" % index_width) % repr(list(treeindex)), ". . " * len(treeindex), repr_element)))

    if asstring: return "\n".join(output)
    else: print "\n".join(output)

  ### how to present SVG objects on the commandline (used in tree)
  def __repr__(self):
    if "repr" in self.__dict__: return self.repr

    output = ["%s" % self.tag]

    remaining = copy.copy(self.attrib)  # shallow copy

    value = remaining.pop("id", None)
    if value is not None:
      output.append("id='%s'" % value)

    # special handling of a text child: print it out and truncate if it's too long
    if self.tag in ("text", "tspan") and len(self.children) == 1 and isinstance(self.children[0], basestring):
      value = re.sub("\n", "\\\\n", self.children[0])
      if len(value) > 13:
        repr_value = "'%s...'" % value[0:10]
      else:
        repr_value = "'%s'" % value
      output.append(repr_value)

    signature = getattr(defaults, "signature_%s" % self.tag, None)
    if signature is not None:
      for name in signature:
        try:
          value = remaining.pop(name)

          # special handling of path data: truncate it if it's too long
          if name == "d":
            if isinstance(value, basestring):
              value = re.sub("\n", "\\\\n", value)
              if len(value) > 13:
                repr_value = "'%s...'" % value[0:10]
              else:
                repr_value = "'%s'" % value

            elif isinstance(value, (list, tuple)):
              if len(value) > 3:
                repr_value = repr(value[0:3])
                repr_value = "%s, ...%s" % (repr_value[0:-1], repr_value[-1])
              else:
                repr_value = repr(value)

            else:
              repr_value = repr(value)

          # special handling of floats: use __str__ instead of __repr__
          elif isinstance(value, float):
            repr_value = "%s" % str(value)

          # normal handling
          else:
            repr_value = repr(value)

          output.append("%s=%s" % (name, repr_value))

        # if the signature item is not present, don't print it
        except KeyError:
          pass

    lenchildren = len(self.children)
    if lenchildren == 1:
      # special handling of a text child: already printed
      if self.tag in ("text", "tspan") and isinstance(self.children[0], basestring):
        pass
      else:
        output.append("(1 child)")

    elif lenchildren > 1:
      output.append("(%d children)" % lenchildren)

    lenremaining = len(remaining)
    if lenremaining == 1:
      output.append("(1 other attribute)")
    elif lenremaining > 1:
      output.append("(%d other attributes)" % lenremaining)

    return "<%s>" % " ".join(output)

  ### convert to XML, view, and save
  def xml(self, indent=u"    ", newl=u"\n"):
    # need a parent node
    if self.tag == "svg" or "_tag" in self.__dict__ and self._tag == "svg":
      svg = self
    else:
      svg = SVG("svg")(self)

    output = [defaults.xml_header] + svg_to_xml(svg, indent) + [u""]
    return unicode(newl.join(output))

  def view(self): # no writing-to-disk needed!
    import _viewer
    _viewer.str(self.xml())

  def save(self, fileName, encoding="utf-8", compresslevel=None):
    fileName = defaults.expand_fileName(fileName)

    if compresslevel is not None or re.search("\.svgz$", fileName, re.I) or re.search("\.gz$", fileName, re.I):
      import gzip
      if compresslevel is None:
        f = gzip.GzipFile(fileName, "w")
      else:
        f = gzip.GzipFile(fileName, "w", compresslevel)

      f = codecs.EncodedFile(f, "utf-16", encoding)
      f.write(self.xml())
      f.close()

    else:
      f = codecs.open(fileName, "w", encoding=encoding)
      f.write(self.xml())
      f.close()

    saved.append(fileName)

  def _write_tempfile(self, fileName=None, encoding="utf-8"):
    if fileName is None:
      fd, fileName = tempfile.mkstemp(".svg", "svgfig-")
      os.close(fd)
    else:
      fileName = defaults.expand_fileName(fileName)

    self.save(fileName, encoding)
    return fileName
    
  def inkview(self, fileName=None, encoding="utf-8"):
    fileName = self._write_tempfile(fileName, encoding)
    os.spawnvp(os.P_NOWAIT, "inkview", ("inkview", fileName))

  def inkscape(self, fileName=None, encoding="utf-8"):
    fileName = self._write_tempfile(fileName, encoding)
    os.spawnvp(os.P_NOWAIT, "inkscape", ("inkscape", fileName))

  def firefox(self, fileName=None, encoding="utf-8"):
    fileName = self._write_tempfile(fileName, encoding)
    os.spawnvp(os.P_NOWAIT, "firefox", ("firefox", fileName))

  ### pickleability and value-based equality
  def __getstate__(self):
    return (sys.version_info, defaults.version_info, self.__dict__)

  def __setstate__(self, state):
    python_version = state[0]
    svgfig_version = state[1]
    if svgfig_version != defaults.version_info:
      warnings.warn("Object created in SVGFig %s, but this is SVGFig %s" % (".".join(map(str, svgfig_version)), ".".join(map(str, defaults.version_info))), defaults.VersionWarning, 5)
    self.__dict__ = state[2]

  def __eq__(self, other):
    if id(self) == id(other): return True
    if self.__class__ != other.__class__: return False
    selfdict = copy.copy(self.__dict__)
    otherdict = copy.copy(other.__dict__)
    del selfdict["_svg"]
    del otherdict["_svg"]
    return selfdict == otherdict

  def __ne__(self, other): return not (self == other)

  def __deepcopy__(self, memo={}):
    output = new.instance(self.__class__)
    output.__dict__ = copy.deepcopy(self.__dict__, memo)
    if "repr" in output.__dict__: del output.__dict__["repr"]
    memo[id(self)] = output
    return output

  ### act like a list
  def append(self, other): self.children.append(other)
  def prepend(self, other): self.children[0:0] = [other]
  def insert(self, i, other): self.children.insert(i, other)
  def remove(self, other): self.children.remove(other)
  def __len__(self): return len(self.children)

  def extend(self, other):
    if isinstance(other, SVG):
      self.children.extend(other.children)
    elif isinstance(other, basestring):
      self.children.append(other)
    else:
      self.children.extend(other)

  def __add__(self, other):
    output = copy.deepcopy(self)
    output += other
    return output

  def __iadd__(self, other):
    self.children.append(other)
    return self

  def __mul__(self, other):
    output = copy.deepcopy(self)
    output *= other
    return output

  def __rmul__(self, other):
    return self * other

  def __imul__(self, other):
    self.children *= other
    return self

  def count(self, *args, **kwds): return self.children.count(*args, **kwds)
  def index(self, *args, **kwds): return self.children.index(*args, **kwds)
  def pop(self, *args, **kwds): return self.children.pop(*args, **kwds)
  def reverse(self, *args, **kwds): return self.children.reverse(*args, **kwds)

  ### act like a dict
  def clear(self, *args, **kwds):
    self.children = []
    self.attrib.clear(*args, **kwds)

  def update(self, other):
    if isinstance(other, SVG):
      self.attrib.update(other.attrib)
    else:
      self.attrib.update(other)

  def __contains__(self, other):
    return other in self.attrib or other in self.children

  def fromkeys(self, *args, **kwds): return self.attrib.fromkeys(*args, **kwds)
  def has_key(self, *args, **kwds): return self.attrib.has_key(*args, **kwds)
  def items(self, *args, **kwds): return self.attrib.items(*args, **kwds)
  def keys(self, *args, **kwds): return self.attrib.keys(*args, **kwds)
  def values(self, *args, **kwds): return self.attrib.values(*args, **kwds)
  def get(self, *args, **kwds): return self.attrib.get(*args, **kwds)
  def setdefault(self, *args, **kwds): return self.attrib.setdefault(*args, **kwds)
  def iteritems(self, *args, **kwds): return self.attrib.iteritems(*args, **kwds)
  def iterkeys(self, *args, **kwds): return self.attrib.iterkeys(*args, **kwds)
  def itervalues(self, *args, **kwds): return self.attrib.itervalues(*args, **kwds)
  def pop(self, *args, **kwds): return self.attrib.pop(*args, **kwds)
  def popitem(self, *args, **kwds): return self.attrib.popitem(*args, **kwds)
  def copy(self): return copy.copy(self)
  def deepcopy(self): return copy.deepcopy(self)

############################### rules for converting into XML

# how to convert SVG objects into XML (as a list of lines to be joined later)
def svg_to_xml(svg, indent, depth=0):
  # if the tag is None, it's a dynamic object that needs to be turned into _svg
  if isinstance(svg, SVG) and svg.tag is None:
    svg.svg()
    svg = svg._svg  # follow that? good.

  if isinstance(svg, basestring):
    return [svg]

  elif isinstance(svg, SVG):
    line = [indent * depth, u"<", svg.tag, u" "]
    remaining = copy.copy(svg.attrib)  # shallow copy that we can pop

    try:
      line.append(u"id=\"%s\" " % remaining.pop("id"))
    except KeyError: pass

    # signature attributes first, for readability
    signature = getattr(defaults, "signature_%s" % svg.tag, None)
    if signature is not None:
      for name in signature:
        try:
          line.append(u"%s=\"%s\" " % (name, attrib_to_xml(svg.tag, name, remaining.pop(name))))
        except KeyError: pass

    remainingkeys = remaining.keys()
    remainingkeys.sort() # for reproducible XML (maybe also helps readability)
    for name in remainingkeys:
      line.append(u"%s=\"%s\" " % (name, attrib_to_xml(svg.tag, name, remaining[name])))

    if len(svg.children) == 0:
      line.append(u"/>")
      return [u"".join(line)]

    else:
      line.append(u">")

      # no indenting for text
      if svg.tag in ("text", "tspan"):
        for i in svg.children:
          line.extend(svg_to_xml(i, indent, 0))
        line.append(u"</%s>" % (svg.tag))
        return [u"".join(line)]

      else:
        lines = [u"".join(line)]
        for i in svg.children:
          lines.extend(svg_to_xml(i, indent, depth+1))
        lines.append(u"%s</%s>" % (indent * depth, svg.tag))
        return lines

  else:
    if type(svg) == types.InstanceType:
      raise TypeError, "SVG contains an unrecognized object: instance of class %s" % svg.__class__.__name__
    else:
      raise TypeError, "SVG contains an unrecognized object: %s" % type(svg)

# how to convert different attribute types into XML
def attrib_to_xml(tag, name, value):
  if isinstance(value, basestring):
    return value

  elif isinstance(value, (int, long, float)):
    return repr(value)  # more precise

  elif isinstance(value, (list, tuple)) and tag == "path" and name == "d":
    def numbertostr(x):
      if isinstance(x, (int, long, float)): return repr(x)  # more precise
      else: return x

    line = []
    lastcommand = None
    for datum in value:
      if not isinstance(datum, (list, tuple)):
        raise TypeError, "Pathdata elements must be lists/tuples"

      command = datum[0]
      args = map(numbertostr, datum[1:])

      if lastcommand == command:
        line.append(u" ")
        line.append(u" ".join(args))
        lastcommand = command
      else:
        line.append(command)
        line.append(u" ".join(args))

      lastcommand = command

    return u"".join(line)

  elif isinstance(value, (list, tuple)):
    line = []
    for v in value:
      line.append(attrib_to_xml(tag, name, v))
    return u", ".join(line)

  elif isinstance(value, dict):
    line = []
    for n, v in value.items():
      line.append(u"%s:%s" % (n, attrib_to_xml(tag, name, v)))
    return u"; ".join(line)

  else:
    return unicode(value)

############################### XML preprocessor instructions and comments

class Instruction(SVG):
  def __init__(self, tag, text):
    self.__dict__["tag"] = tag
    self.__dict__["attrib"] = {}
    self.__dict__["children"] = []
    self.__dict__["_svg"] = self
    self.__dict__["text"] = text

  def xml(self): return "<?%s %s?>" % (self.tag, self.text)

  def __repr__(self):
    value = re.sub("\n", "\\\\n", self.text)
    if len(value) > 23:
      return "<?%s %s...?>" % (self.tag, value[0:20])
    else:
      return "<?%s %s?>" % (self.tag, value)

class Comment(SVG):
  def __init__(self, text):
    if text.find("--") != -1:
      raise ValueError, "SVG comments must not include '--'"
    self.__dict__["tag"] = "comment"
    self.__dict__["attrib"] = {}
    self.__dict__["children"] = []
    self.__dict__["_svg"] = self
    self.__dict__["text"] = text
  
  def __eq__(self, other): return SVG.__eq__(self, other) and self.text == other.text

  def xml(self): return "<!-- %s -->" % self.text

  def __repr__(self):
    value = re.sub("\n", "\\\\n", self.text)
    if len(value) > 23:
      return "<!-- %s... -->" % value[0:20]
    else:
      return "<!-- %s -->" % value

class CDATA(SVG):
  def __init__(self, text):
    if text.find("]]>") != -1:
      raise ValueError, "CDATA must not include ']]>'"
    self.__dict__["tag"] = "CDATA"
    self.__dict__["attrib"] = {}
    self.__dict__["children"] = []
    self.__dict__["_svg"] = self
    self.__dict__["text"] = text
  
  def __eq__(self, other): return SVG.__eq__(self, other) and self.text == other.text

  def xml(self): return "<![CDATA[%s]]>" % self.text

  def __repr__(self):
    value = re.sub("\n", "\\\\n", self.text)
    if len(value) > 23:
      return "<![CDATA[%s...]]>" % value[0:20]
    else:
      return "<![CDATA[%s]]>" % value

############################### reading SVG from a file

def load(fileName):
  if re.search("\.svgz$", fileName, re.I) or re.search("\.gz$", fileName, re.I):
    import gzip
    f = gzip.GzipFile(fileName)
  else:
    f = file(fileName)
  return load_stream(f)

def template(fileName, svg, replaceme="REPLACEME"):
  output = load(fileName)
  for treeindex, i in output.walk():
    if isinstance(i, SVG) and i.tag == replaceme:
      output[treeindex] = svg
  return output

def load_stream(stream):
  from xml.sax import handler, make_parser
  from xml.sax.handler import feature_namespaces, feature_external_ges, feature_external_pes

  # With a little thought, this could be streamlined.  In its current
  # state, it works (with processing instructions, comments, and CDATA).
  class ContentHandler(handler.ContentHandler):
    def __init__(self):
      self.stack = []
      self.output = None
      self.all_whitespace = re.compile("^\s*$")
      self.CDATA = False

    def startElement(self, tag, attrib):
      s = SVG(tag)
      s.attrib = dict(attrib.items())
      if len(self.stack) > 0:
        last = self.stack[-1]
        last.children.append(s)
      self.stack.append(s)

    def characters(self, ch):
      if self.CDATA:
        last = self.stack[-1]
        last.text += ch

      else:
        if not isinstance(ch, basestring) or self.all_whitespace.match(ch) is None:
          if len(self.stack) > 0:
            last = self.stack[-1]
            if len(last.children) > 0 and isinstance(last.children[-1], basestring):
              last.children[-1] = last.children[-1] + "\n" + ch
            else:
              last.children.append(ch)

    def endElement(self, tag):
      if len(self.stack) > 0:
        last = self.stack[-1]
      self.output = self.stack.pop()

    # If a processing instruction is outside the main <svg> tag, it will be lost.
    def processingInstruction(self, target, data):
      s = Instruction(target, data)
      if len(self.stack) > 0:
        last = self.stack[-1]
        last.children.append(s)
      self.output = s

    def comment(self, comment):
      s = Comment(re.sub("(^ | $)", "", comment))
      if len(self.stack) > 0:
        last = self.stack[-1]
        last.children.append(s)
      self.output = s

    def startCDATA(self):
      s = CDATA("")
      if len(self.stack) > 0:
        last = self.stack[-1]
        last.children.append(s)
      self.stack.append(s)
      self.CDATA = True

    def endCDATA(self):
      if len(self.stack) > 0:
        last = self.stack[-1]
      self.output = self.stack.pop()
      self.CDATA = False

    def startDTD(self, name, public_id, system_id): pass
    def endDTD(self): pass
    def startEntity(self, name): pass
    def endEntity(self, name): pass

  ch = ContentHandler()
  parser = make_parser()
  parser.setContentHandler(ch)
  parser.setProperty(handler.property_lexical_handler, ch)
  parser.setFeature(feature_namespaces, 0)
  parser.setFeature(feature_external_ges, 0)
  parser.parse(stream)
  return ch.output

############################### standard representation for transformations and parametric functions

def cannonical_transformation(expr):
  if expr is None:
    output = lambda x, y: (x, y)
    output.func_name = "identity"
    return output

  elif callable(expr):

    # 2 real -> 2 real
    if expr.func_code.co_argcount == 2:
      return expr

    # complex -> complex
    elif expr.func_code.co_argcount == 1:
      split = lambda z: (z.real, z.imag)
      output = lambda x, y: split(expr(complex(x, y)))
      output.func_name = expr.func_name
      return output

    else:
      raise TypeError, "Must be a 2 -> 2 real function or a complex -> complex function"

  else:
    compiled = compile(expr, expr, "eval")

    # 2 real -> 2 real
    if "x" in compiled.co_names and "y" in compiled.co_names:
      evalexpr = expr
      evalexpr = re.sub("x", "float(x)", evalexpr)
      evalexpr = re.sub("y", "float(y)", evalexpr)
      output = eval("lambda x,y: (%s)" % evalexpr, math.__dict__)
      output.func_name = "x, y -> %s" % expr
      return output

    # complex -> complex
    elif "z" in compiled.co_names:
      evalexpr = re.sub("z", "complex(x,y)", expr)
      output = eval("lambda x,y: ((%s).real, (%s).imag)" % (evalexpr, evalexpr), cmath.__dict__)
      output.func_name = "z -> %s" % expr
      return output

    else:
      raise TypeError, "Transformation string '%s' must contain real 'x' and 'y' or complex 'z'" % expr

def cannonical_parametric(expr):
  if callable(expr):

    # 1 real -> 2 real
    if expr.func_code.co_argcount == 1:
      return expr

    else:
      raise TypeError, "Must be a 1 -> 2 real function"

  else:
    compiled = compile(expr, expr, "eval")

    # 1 real -> 2 real
    if "t" in compiled.co_names:
      output = eval("lambda t: (%s)" % re.sub("t", "float(t)", expr), math.__dict__)
      output.func_name = "t -> %s" % expr
      return output

    # 1 real -> 1 real
    elif "x" in compiled.co_names:
      output = eval("lambda t: (t, %s)" % re.sub("x", "float(t)", expr), math.__dict__)
      output.func_name = "x -> %s" % expr
      return output

    # real (a complex number restricted to the real axis) -> complex
    elif "z" in compiled.co_names:
      evalexpr = re.sub("z", "complex(t,0)", expr)
      output = eval("lambda t: ((%s).real, (%s).imag)" % (evalexpr, evalexpr), cmath.__dict__)
      output.func_name = "z -> %s" % expr
      return output

    else:
      raise TypeError, "Parametric string '%s' must contain real 't', 'x', or 'z'" % expr

############################### make code objects pickleable (so that curves and transformations are pickleable)

def _code_constructor(code_args, python_version, svgfig_version):
  if python_version != sys.version_info or svgfig_version != defaults.version_info:
    warnings.warn("Function created in Python %s/SVGFig %s, but this is Python %s/SVGFig %s"
                  % (".".join(map(str, python_version)), ".".join(map(str, svgfig_version)), ".".join(map(str, sys.version_info)), ".".join(map(str, defaults.version_info))),
                  defaults.VersionWarning, 5)
  return new.code(*code_args)

def _code_serializer(code):
  if code.co_freevars or code.co_cellvars:
    raise ValueError, "Sorry, can't pickle code that depends on local variables %s %s" % (str(code.co_freevars), str(code.co_cellvars))
  return _code_constructor, ((code.co_argcount, code.co_nlocals, code.co_stacksize,
                              code.co_flags, code.co_code, code.co_consts, code.co_names,
                              code.co_varnames, code.co_filename, code.co_name,
                              code.co_firstlineno, code.co_lnotab),
                             sys.version_info, defaults.version_info)


copy_reg.pickle(types.CodeType, _code_serializer)
