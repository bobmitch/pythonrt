
from math import sin, cos, pi, radians, degrees

def zrot (vertex, origin, degrees):
        r = radians(degrees)
        tx = vertex[0] - origin[0]
        ty = vertex[1] - origin[1]
        x = (tx * cos(r)) - (ty * sin(r))
        y = (tx * sin(r)) + (ty * cos(r))
        x = x + origin[0]
        y = y + origin[1]
        return [x,y,vertex[2]]
    
def yrot (vertex, origin, degrees):
        r = radians(degrees)
        tx = vertex[0] - origin[0]
        tz = vertex[2] - origin[2]
        x = (tz * sin(r)) + (tx * cos(r))
        z = (tz * cos(r)) - (tx * sin(r))
        x = x + origin[0]
        z = z + origin[2]
        return [x,vertex[1],z]

def xrot (vertex, origin, degrees):
        r = radians(degrees)
        ty = vertex[1] - origin[1]
        tz = vertex[2] - origin[2]
        z = (ty * sin(r)) + (tz * cos(r))
        y = (ty * cos(r)) - (tz * sin(r))
        y = y + origin[1]
        z = z + origin[2]
        return [vertex[0],y,z]
