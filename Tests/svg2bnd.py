"""svg2bnd.py -- Converts a path in an SVG file to a boundary for Voro++."""

import sys,re
from svgfig import svg, pathdata
from math import sqrt

if len(sys.argv) < 2:
    print 'Usage: svg2bnd <file>'
    sys.exit(0)
stage = 0
count = 0
obj = svg.load(sys.argv[1])
semitotal = 0
for key, elem in obj.walk():
    if isinstance(elem, svg.SVG):
        if elem.tag == 'path':
            print '# Start'
            redo = 0
            commandStr = elem[u'd']
            commandStr = re.sub('([0-9])-', ' -', commandStr)
            # HACK: Add some spaces to ease along parsing.
            commands = pathdata.parse(commandStr)
            lx, ly = 0.0, 0.0
            x, y = 0.0, 0.0
            xs = []
            ys = []
            for c in commands:
                d = c[0]
                
                if d == u'L' or d == u'M': # Absolute position line
                    x, y = c[1], c[2]
                elif d == u'l' or d == u'm': # Relative position line
                    dx, dy = c[1], c[2]
                    x += dx
                    y += dy
                elif d == u'H': # Absolute position horizontal motion
                    x = c[1]
                elif d == u'h': # Relative position horizontal motion
                    x += c[1]
                elif d == u'V': # Absolute position vertical motion
                    y = c[1]
                elif d == u'v': # Relative position vertical motion
                    y += c[1]
                if (redo == 0 or (sqrt((x-xs[-1])**2+(y-ys[-1])**2) > 1e-5)):
                    xs.append(x)
                    ys.append(y)
                    count = count + 1
                    redo = redo + 1
            if sqrt((xs[-1]-xs[0])**2+(ys[-1]-ys[0])**2) < 1e-5:
                xs = xs[:-1]
                ys = ys[:-1]
                count = count - 1
            for i in xrange(len(xs)):
                print i+semitotal, xs[i], ys[i]

            print '# End'
            semitotal = count
        elif elem.tag == 'polygon':
            redo = 0
            print '# Start'
            pointsStr = elem['points']
            pointses = pointsStr.split()
            xs = []
            ys = []
            for p in pointses:
                x,y = p.split(u',')
                if (redo ==0 or (sqrt((x-xs[-1])**2+(y-ys[-1])**2) > 1e-5)):
                    xs.append(x)
                    ys.append(y)
                    count = count + 1
                    redo = redo + 1
            if sqrt((xs[-1]-xs[0])**2+(ys[-1]-ys[0])**2) < 1e-5:
                xs = xs[:-1]
                ys = ys[:-1]
                count = count - 1
            for i in xrange(len(xs)):
                print i + semitotal, xs[i], ys[i]
                
            print '# End'
            semittotal=count
#        else:
#            print 'Unsupported element tag: %s'%elem.tag
