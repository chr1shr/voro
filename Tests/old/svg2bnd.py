"""svg2bnd.py -- Converts a path in an SVG file to a boundary for Voro++."""

import sys
from svgfig import *

if len(sys.argv) < 2:
    print 'Usage: svg2bnd <file>'
    sys.exit(0)

obj = load(sys.argv[1])
for e in obj.values():
    if isinstance(e, SVG) and e.t == u'path':
        print '# Start'
        path = pathtoPath(e)
        commands = path.d
        count= 0 
        x,y = 0.0, 0.0
        print count,x,y
        for c in commands[1:]:
            if c[0] == u'm':
                count += 1
                lx=x
                ly=y
                x += c[1]
                y += c[2]
                if abs((x-lx)**2 + (y-ly)**2) > 1e-8:
                    print count,x,y
        print '# End'

