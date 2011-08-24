import defaults

############################### convenient functions for making paths

def poly(*data, **kwds):
  errstring = "Arguments are: poly((x1,y1), (x2,y2), ..., loop=False)"
  loop = False
  if "loop" in kwds:
    loop = kwds["loop"]
    del kwds["loop"]
  if len(kwds) > 0: raise TypeError, errstring

  try:
    output = []
    for x, y in data:
      if output == []: output.append(("M", x, y))
      else: output.append(("L", x, y))
    if loop and len(data) > 0:
      output.append(("Z",))
    return output
  except (TypeError, ValueError): raise TypeError, errstring

def bezier(*data, **kwds):
  errstring = "Arguments are: bezier((x,y,c1x,c1y,c2x,c2y), ..., loop=False)"
  loop = False
  if "loop" in kwds:
    loop = kwds["loop"]
    del kwds["loop"]
  if len(kwds) > 0: raise TypeError, errstring

  try:
    output = []
    for x, y, c1x, c1y, c2x, c2y in data:
      if output == []: output.append(("M", x, y))
      else:
        output.append(("C", c1x, c1y, c2x, c2y, x, y))
    if loop and len(data) > 0:
      output.append(("Z",))
    return output
  except (TypeError, ValueError): raise TypeError, errstring

def velocity(*data, **kwds):
  errstring = "Arguments are: velocity((x,y,vx,vy), ..., loop=False)"
  loop = False
  if "loop" in kwds:
    loop = kwds["loop"]
    del kwds["loop"]
  if len(kwds) > 0: raise TypeError, errstring

  try:
    output = []
    indexes = range(len(data))
    if loop and len(data) > 0: indexes.append(0)

    for i in indexes:
      if output == []: output.append(("M", data[i][0], data[i][1]))
      else:
        inext = (i+1) % len(data)
        iprev = (i-1) % len(data)

        x, y = data[i][0], data[i][1]
        c1x, c1y = data[iprev][2]/3. + data[iprev][0], data[iprev][3]/3. + data[iprev][1]
        c2x, c2y = data[i][2]/-3. + x, data[i][3]/-3. + y

        output.append(("C", c1x, c1y, c2x, c2y, x, y))

    if loop and len(data) > 0:
      output.append(("Z",))
    return output
  except (TypeError, ValueError): raise TypeError, errstring

def foreback(*data, **kwds):
  errstring = "Arguments are: foreback((x,y,vfx,vfy,vbx,vby), ..., loop=False)"
  loop = False
  if "loop" in kwds:
    loop = kwds["loop"]
    del kwds["loop"]
  if len(kwds) > 0: raise TypeError, errstring

  try:
    output = []
    indexes = range(len(data))
    if loop and len(data) > 0: indexes.append(0)

    for i in indexes:
      if output == []: output.append(("M", data[i][0], data[i][1]))
      else:
        inext = (i+1) % len(data)
        iprev = (i-1) % len(data)

        x, y = data[i][0], data[i][1]
        c1x, c1y = data[iprev][4]/3. + data[iprev][0], data[iprev][5]/3. + data[iprev][1]
        c2x, c2y = data[i][2]/-3. + x, data[i][3]/-3. + y

        output.append(("C", c1x, c1y, c2x, c2y, x, y))

    if loop and len(data) > 0:
      output.append(("Z",))
    return output
  except (TypeError, ValueError): raise TypeError, errstring

def smooth(*data, **kwds):
  errstring = "Arguments are: smooth((x1,y1), (x2,y2), ..., loop=False)"

  loop = False
  if "loop" in kwds:
    loop = kwds["loop"]
    del kwds["loop"]
  if len(kwds) > 0: raise TypeError, errstring

  try:
    x, y = zip(*data)
    vx, vy = [0.]*len(data), [0.]*len(data)
    for i in xrange(len(data)):
      inext = (i+1) % len(data)
      iprev = (i-1) % len(data)

      vx[i] = (x[inext] - x[iprev])/2.
      vy[i] = (y[inext] - y[iprev])/2.
      if not loop and (i == 0 or i == len(data)-1):
        vx[i], vy[i] = 0., 0.

    return velocity(zip(x, y, vx, vy), loop)
  except (TypeError, ValueError): raise TypeError, errstring

############################### pathdata parsers

def parse_whitespace(index, pathdata):
  while index < len(pathdata) and pathdata[index] in (" ", "\t", "\r", "\n", ","): index += 1
  return index, pathdata

def parse_command(index, pathdata):
  index, pathdata = parse_whitespace(index, pathdata)

  if index >= len(pathdata): return None, index, pathdata
  command = pathdata[index]
  if "A" <= command <= "Z" or "a" <= command <= "z":
    index += 1
    return command, index, pathdata
  else: 
    return None, index, pathdata

def parse_number(index, pathdata):
  index, pathdata = parse_whitespace(index, pathdata)

  if index >= len(pathdata): return None, index, pathdata
  first_digit = pathdata[index]

  if "0" <= first_digit <= "9" or first_digit in ("-", "+", "."):
    start = index
    while index < len(pathdata) and ("0" <= pathdata[index] <= "9" or pathdata[index] in ("-", "+", ".", "e", "E")):
      index += 1
    end = index

    index = end
    return float(pathdata[start:end]), index, pathdata
  else: 
    return None, index, pathdata

def parse_boolean(index, pathdata):
  index, pathdata = parse_whitespace(index, pathdata)

  if index >= len(pathdata): return None, index, pathdata
  first_digit = pathdata[index]

  if first_digit in ("0", "1"):
    index += 1
    return int(first_digit), index, pathdata
  else:
    return None, index, pathdata

############################### main parsing function (keeps defaults from getting messy)

def parse(pathdata):
  if isinstance(pathdata, (list, tuple)): return pathdata

  output = []
  index = 0
  while True:
    command, index, pathdata = parse_command(index, pathdata)
    index, pathdata = parse_whitespace(index, pathdata)

    if command == None and index == len(pathdata): break  # this is the normal way out of the loop
    if command in ("Z", "z"):
      output.append((command,))

    ######################
    elif command in ("H", "h", "V", "v"):
      errstring = "Pathdata command \"%s\" requires a number at index %d" % (command, index)
      num1, index, pathdata = parse_number(index, pathdata)
      if num1 == None: raise ValueError, errstring

      while num1 != None:
        output.append((command, num1))
        num1, index, pathdata = parse_number(index, pathdata)

    ######################
    elif command in ("M", "m", "L", "l", "T", "t"):
      errstring = "Pathdata command \"%s\" requires an x,y pair at index %d" % (command, index)
      num1, index, pathdata = parse_number(index, pathdata)
      num2, index, pathdata = parse_number(index, pathdata)

      if num1 == None: raise ValueError, errstring

      while num1 != None:
        if num2 == None: raise ValueError, errstring
        output.append((command, num1, num2))

        num1, index, pathdata = parse_number(index, pathdata)
        num2, index, pathdata = parse_number(index, pathdata)

    ######################
    elif command in ("S", "s", "Q", "q"):
      errstring = "Pathdata command \"%s\" requires a cx,cy,x,y quadruplet at index %d" % (command, index)
      num1, index, pathdata = parse_number(index, pathdata)
      num2, index, pathdata = parse_number(index, pathdata)
      num3, index, pathdata = parse_number(index, pathdata)
      num4, index, pathdata = parse_number(index, pathdata)

      if num1 == None: raise ValueError, errstring

      while num1 != None:
        if num2 == None or num3 == None or num4 == None: raise ValueError, errstring
        output.append((command, num1, num2, num3, num4))

        num1, index, pathdata = parse_number(index, pathdata)
        num2, index, pathdata = parse_number(index, pathdata)
        num3, index, pathdata = parse_number(index, pathdata)
        num4, index, pathdata = parse_number(index, pathdata)

    ######################
    elif command in ("C", "c"):
      errstring = "Pathdata command \"%s\" requires a c1x,c1y,c2x,c2y,x,y sextuplet at index %d" % (command, index)
      num1, index, pathdata = parse_number(index, pathdata)
      num2, index, pathdata = parse_number(index, pathdata)
      num3, index, pathdata = parse_number(index, pathdata)
      num4, index, pathdata = parse_number(index, pathdata)
      num5, index, pathdata = parse_number(index, pathdata)
      num6, index, pathdata = parse_number(index, pathdata)

      if num1 == None: raise ValueError, errstring

      while num1 != None:
        if num2 == None or num3 == None or num4 == None or num5 == None or num6 == None: raise ValueError, errstring

        output.append((command, num1, num2, num3, num4, num5, num6))

        num1, index, pathdata = parse_number(index, pathdata)
        num2, index, pathdata = parse_number(index, pathdata)
        num3, index, pathdata = parse_number(index, pathdata)
        num4, index, pathdata = parse_number(index, pathdata)
        num5, index, pathdata = parse_number(index, pathdata)
        num6, index, pathdata = parse_number(index, pathdata)

    ######################
    elif command in ("A", "a"):
      errstring = "Pathdata command \"%s\" requires a rx,ry,angle,large-arc-flag,sweep-flag,x,y septuplet at index %d" % (command, index)
      num1, index, pathdata = parse_number(index, pathdata)
      num2, index, pathdata = parse_number(index, pathdata)
      num3, index, pathdata = parse_number(index, pathdata)
      num4, index, pathdata = parse_boolean(index, pathdata)
      num5, index, pathdata = parse_boolean(index, pathdata)
      num6, index, pathdata = parse_number(index, pathdata)
      num7, index, pathdata = parse_number(index, pathdata)

      if num1 == None: raise ValueError, errstring

      while num1 != None:
        if num2 == None or num3 == None or num4 == None or num5 == None or num6 == None or num7 == None: raise ValueError, errstring

        output.append((command, num1, num2, num3, num4, num5, num6, num7))

        num1, index, pathdata = parse_number(index, pathdata)
        num2, index, pathdata = parse_number(index, pathdata)
        num3, index, pathdata = parse_number(index, pathdata)
        num4, index, pathdata = parse_boolean(index, pathdata)
        num5, index, pathdata = parse_boolean(index, pathdata)
        num6, index, pathdata = parse_number(index, pathdata)
        num7, index, pathdata = parse_number(index, pathdata)

  return output

############################### transformation function (keeps defaults from getting messy)

def transform(func, pathdata):
  x, y, X, Y = None, None, None, None
  output = []
  for datum in pathdata:
    if not isinstance(datum, (tuple, list)):
      raise TypeError, "Pathdata elements must be lists/tuples"

    command = datum[0]
    args = datum[1:]

    ######################
    if command in ("Z", "z"):
      x, y, X, Y = None, None, None, None
      output.append(("Z",))

    ######################
    elif command in ("H", "h", "V", "v"):
      num1 = args[0]

      if command == "H" or (command == "h" and x == None): x = num1
      elif command == "h": x += num1
      elif command == "V" or (command == "v" and y == None): y = num1
      elif command == "v": y += num1

      X, Y = func(x, y)
      output.append(("L", X, Y))

    ######################
    elif command in ("M", "m", "L", "l", "T", "t"):
      num1, num2 = args

      if command.isupper() or x == None or y == None:
        x, y = num1, num2
      else:
        x += num1
        y += num2

      X, Y = func(x, y)
      output.append((command.capitalize(), X, Y))

    ######################
    elif command in ("S", "s", "Q", "q"):
      num1, num2, num3, num4 = args

      if command.isupper() or x == None or y == None:
        cx, cy = num1, num2
      else:
        cx = x + num1
        cy = y + num2

      if command.isupper() or x == None or y == None:
        x, y = num3, num4
      else:
        x += num3
        y += num4

      CX, CY = func(cx, cy)
      X, Y = func(x, y)
      output.append((command.capitalize(), CX, CY, X, Y))

    ######################
    elif command in ("C", "c"):
      num1, num2, num3, num4, num5, num6 = args

      if command.isupper() or x == None or y == None:
        c1x, c1y = num1, num2
      else:
        c1x = x + num1
        c1y = y + num2

      if command.isupper() or x == None or y == None:
        c2x, c2y = num3, num4
      else:
        c2x = x + num3
        c2y = y + num4

      if command.isupper() or x == None or y == None:
        x, y = num5, num6
      else:
        x += num5
        y += num6

      C1X, C1Y = func(c1x, c1y)
      C2X, C2Y = func(c2x, c2y)
      X, Y = func(x, y)
      output.append((command.capitalize(), C1X, C1Y, C2X, C2Y, X, Y))

    ######################
    elif command in ("A", "a"):
      num1, num2, angle, large_arc_flag, sweep_flag, num3, num4 = args

      oldx, oldy = x, y
      OLDX, OLDY = X, Y

      if command.isupper() or x == None or y == None:
        x, y = num3, num4
      else:
        x += num3
        y += num4
      X, Y = func(x, y)

      if x != None and y != None:
        centerx, centery = (x + oldx)/2., (y + oldy)/2.
      CENTERX, CENTERY = (X + OLDX)/2., (Y + OLDY)/2.

      rx = centerx + num1
      ry = centery + num2
      RX, RY = func(rx, ry)

      output.append((command.capitalize(), RX - CENTERX, RY - CENTERY, angle, large_arc_flag, sweep_flag, X, Y))

  return output

############################### bbox function (keeps defaults from getting messy)

def bbox(pathdata):
  x, y = None, None
  output = defaults.BBox(None, None, None, None)

  for datum in pathdata:
    if not isinstance(datum, (tuple, list)):
      raise TypeError, "Pathdata elements must be lists/tuples"

    command = datum[0]
    args = datum[1:]

    ######################
    if command in ("Z", "z"): pass

    ######################
    elif command in ("H", "h", "V", "v"):
      num1 = args[0]

      if command == "H" or (command == "h" and x == None): x = num1
      elif command == "h": x += num1
      elif command == "V" or (command == "v" and y == None): y = num1
      elif command == "v": y += num1

      output.insert(x, y)

    ######################
    elif command in ("M", "m", "L", "l", "T", "t"):
      num1, num2 = args

      if command.isupper() or x == None or y == None:
        x, y = num1, num2
      else:
        x += num1
        y += num2

      output.insert(x, y)

    ######################
    elif command in ("S", "s", "Q", "q"):
      num1, num2, num3, num4 = args

      if command.isupper() or x == None or y == None:
        cx, cy = num1, num2
      else:
        cx = x + num1
        cy = y + num2

      if command.isupper() or x == None or y == None:
        x, y = num3, num4
      else:
        x += num3
        y += num4

      output.insert(x, y)

    ######################
    elif command in ("C", "c"):
      num1, num2, num3, num4, num5, num6 = args

      if command.isupper() or x == None or y == None:
        c1x, c1y = num1, num2
      else:
        c1x = x + num1
        c1y = y + num2

      if command.isupper() or x == None or y == None:
        c2x, c2y = num3, num4
      else:
        c2x = x + num3
        c2y = y + num4

      if command.isupper() or x == None or y == None:
        x, y = num5, num6
      else:
        x += num5
        y += num6

      output.insert(x, y)

    ######################
    elif command in ("A", "a"):
      num1, num2, angle, large_arc_flag, sweep_flag, num3, num4 = args

      oldx, oldy = x, y
      OLDX, OLDY = X, Y

      if command.isupper() or x == None or y == None:
        x, y = num3, num4
      else:
        x += num3
        y += num4

      if x != None and y != None:
        centerx, centery = (x + oldx)/2., (y + oldy)/2.
      CENTERX, CENTERY = (X + OLDX)/2., (Y + OLDY)/2.

      output.insert(x, y)

  return output
  
