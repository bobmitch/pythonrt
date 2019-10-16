""" Basic surfaces module
    Includes  Sphere class, Ray class, YPlane class
"""

from math import sqrt
from rotate import *
from vector import *
#import scipy


import copy

INFINITY = 999999999
EPSILON = 0.00001


    
class Ray():
    def __init__ (self, o=[0,0,0], d=[0,0,1]):
        self.o = o    # origin
        self.d = d    # direction (nominally a unit vector, but not always)
        
                    
    



class Material():
    def __init__ (self):
        self.colour=[1.0,0,0]
        self.reflect=0.0
        self.spec_level = 0.6
        self.spec_size = 16
        self.reflection_blur = 0.0
        self.refract = 0
        self.ior = 1.0
        self.send_caustics = False
        self.receive_caustics = False




class Entity ():
    def __init__ (self, center=[0,0,0], aabbmin=[-1,-1,-1], aabbmax=[1,1,1], mat=Material()):
        self.center=center
        self.aabbmin = aabbmin
        self.aabbmax = aabbmax
        self.mat = Material()
        self.photon_list = []
     
    def slow_ray_aabb_test (self, ray):
        tnear = -INFINITY
        tfar = INFINITY
        
        #TODO: UNROLL FOR LOOP AND REMOVED DOTS
        for p in xrange(3):
            # for each plane of bounding box
            if ray.dir.v[p] == 0.0: #if direction of ray is parallel to the axis
                if ray.orig.v[p] < self.aabbmin[p] or ray.orig.v[p] > self.aabbmax[p]:
                    return (False)
            else:
                t1 = (self.aabbmin[p] - ray.orig.v[p]) / ray.dir.v[p]
                t2 = (self.aabbmax[p] - ray.orig.v[p]) / ray.dir.v[p]
                if t1 > t2: #swap t1 and t2
                    temp = t1
                    t1 = t2
                    t2 = temp
                if t1 > tnear: tnear = t1
                if t2 < tfar: tfar = t2
                if tnear > tfar: return (False)
                if tfar < 0: return (False)
        #print "Hit something!"
        return (True) #Intersection with AABB occurs
    
       
    def fast_ray_aabb_test (self, ray):
        dx, dy, dz = ray.d[0], ray.d[1], ray.d[2]
        ox, oy, oz = ray.o[0], ray.o[1], ray.o[2]
        minx, miny, minz = self.aabbmin[0], self.aabbmin[1], self.aabbmin[2]
        maxx, maxy, maxz = self.aabbmax[0], self.aabbmax[1], self.aabbmax[2]
        
        if dx >= 0 and maxx < ox : return (False)
        if dx <= 0 and minx >= ox: return (False)
        if dy >= 0 and maxy < oy : return (False)
        if dy <= 0 and miny >= oy: return (False)
        if dz >= 0 and maxz < oz : return (False)
        if dz <= 0 and minz >= oz: return (False)
        
        #if this point reached, then ray is somewhat facing the aabb
        return (True)

    
    
       
        







class Plane():
    def __init__ (self, normal=[0,1,0], offset=-200):
        self.infinite=True
        self.offset=offset
        self.normal = normal
        self.mat = Material()
        self.photons=[]
        #self.p = vec_mul (self.normal , offset)
        
    def ray_intersect (self,ray):
        
        #v = self.normal.pEscalar(r.direccion)
        #if (v == 0): return False
        #
        #distChoque = -(self.normal.pEscalar(r.origen)+self.distancia)/v
        #if (distChoque<0): return False             # Direccion del rayo negativa
        #if (distChoque > r.disInter): return False  # No es el mas cercano
        #
        #r.disInter = distChoque
        
        v = vec_dot (self.normal, ray.d)
        if v == 0: return INFINITY, [-1,-1,-1], 0, False
        t =  -(vec_dot(self.normal,ray.o) + self.offset) / v
        if t < 0 : return INFINITY, [-1,-1,-1], 0, False     # no intersection (occurs behind ray origin)
        pt = vec_add ( ray.o , vec_mul (ray.d,t) )
        return t, pt, self.normal, False


class YPlane():
    def __init__ (self):
        self.infinite=True
        self.normal = [0,1,0]
        self.mat = Material()
        self.center = [0,0,0]
        self.photon_list = []
        
    def ray_intersect (self,ray):
        if (ray.o[1]>0 and ray.d[1]<0) or (ray.o[1]<0 and ray.d[1]<0):
            return INFINITY, [-1,-1,-1],0
        if ray.d[1]==0:
            # ray is parallel to Y plane
            return INFINITY, [-1,-1,-1], 0
        v = vec_dot (self.normal, ray.d)
        t = (-vec_dot (self.normal, ray.o))/v
        pt = vec_add ( ray.o , vec_mul (ray.d,t) )
        
        if t<0:
            return INFINITY,[-1,-1,-1],0
        
        # checkers!
        #r = (1+sin(pt[0]/21))/2
        #g = (1+sin(pt[2]/23))/2
        #b = (1+cos(pt[2]/17))/2
        #self.mat.colour = (r,g,b)
        #if pt[0]%100 < 50 and pt[2]%100 < 50:
        #    self.mat.colour=[0,0,0]
        #else:
        #    self.mat.colour=[1,1,1]
        #
        # kludge for double sided plane
        if v<0: return t,pt,[0,1,0], False
        return t, pt, [0,-1,0], False


class Tri (Entity):
    def __init__ (self, points=[[-10,-10,10],[0,10,10],[10,-10,10]]):
        self.p = points
        self.flat = True
        self.vertex_normal = [ [1,0,0],[1,0,0],[1,0,0]]
        # calc bounding box to pass to entity parent class
        minx, miny, minz = INFINITY, INFINITY, INFINITY
        maxx, maxy, maxz = -INFINITY, -INFINITY, -INFINITY     
        if self.p[0][0]<minx: minx = self.p[0][0]
        if self.p[0][1]<miny: miny = self.p[0][1]
        if self.p[0][2]<minz: minz = self.p[0][2]
        if self.p[1][0]<minx: minx = self.p[1][0]
        if self.p[1][1]<miny: miny = self.p[1][1]
        if self.p[1][2]<minz: minz = self.p[1][2]
        if self.p[2][0]<minx: minx = self.p[2][0]
        if self.p[2][1]<miny: miny = self.p[2][1]
        if self.p[2][2]<minz: minz = self.p[2][2]
        if self.p[0][0]>maxx: maxx = self.p[0][0]
        if self.p[0][1]>maxy: maxy = self.p[0][1]
        if self.p[0][2]>maxz: maxz = self.p[0][2]
        if self.p[1][0]>maxx: maxx = self.p[1][0]
        if self.p[1][1]>maxy: maxy = self.p[1][1]
        if self.p[1][2]>maxz: maxz = self.p[1][2]
        if self.p[2][0]>maxx: maxx = self.p[2][0]
        if self.p[2][1]>maxy: maxy = self.p[2][1]
        if self.p[2][2]>maxz: maxz = self.p[2][2]
        
        # TODO: Remove this ugly hack to avoid craziness when sorting in kdtree by making aabb little bigger than needed
        self.aabbmin = [minx-0.1, miny-0.1, minz-0.1]
        self.aabbmax = [maxx+0.1, maxy+0.1, maxz+0.1]
        
        center = vec_add (self.aabbmin , vec_sub (self.aabbmax,self.aabbmin) )
        
        Entity.__init__ (self, center, [minx-EPSILON, miny-EPSILON, minz-EPSILON], [maxx+EPSILON, maxy+EPSILON, maxz+EPSILON] )
        self.recalc ()
        
    def flip (self):
        """ reorded the vertices in the triangle - effectively reversing the normal """
        temp = copy.deepcopy (self.p[0])
        self.p[0] = copy.deepcopy (self.p[2])
        self.p[2] = copy.deepcopy (temp)
        
        
    def recalc (self):
        # calc bounding box to pass to entity parent class
        minx, miny, minz = INFINITY, INFINITY, INFINITY
        maxx, maxy, maxz = -INFINITY, -INFINITY, -INFINITY     
        if self.p[0][0]<minx: minx = self.p[0][0]
        if self.p[0][1]<miny: miny = self.p[0][1]
        if self.p[0][2]<minz: minz = self.p[0][2]
        if self.p[1][0]<minx: minx = self.p[1][0]
        if self.p[1][1]<miny: miny = self.p[1][1]
        if self.p[1][2]<minz: minz = self.p[1][2]
        if self.p[2][0]<minx: minx = self.p[2][0]
        if self.p[2][1]<miny: miny = self.p[2][1]
        if self.p[2][2]<minz: minz = self.p[2][2]
        if self.p[0][0]>maxx: maxx = self.p[0][0]
        if self.p[0][1]>maxy: maxy = self.p[0][1]
        if self.p[0][2]>maxz: maxz = self.p[0][2]
        if self.p[1][0]>maxx: maxx = self.p[1][0]
        if self.p[1][1]>maxy: maxy = self.p[1][1]
        if self.p[1][2]>maxz: maxz = self.p[1][2]
        if self.p[2][0]>maxx: maxx = self.p[2][0]
        if self.p[2][1]>maxy: maxy = self.p[2][1]
        if self.p[2][2]>maxz: maxz = self.p[2][2]
        
        # TODO: Remove this ugly hack to avoid craziness when sorting in kdtree by making aabb little bigger than needed
        self.aabbmin = [minx-0.1, miny-0.1, minz-0.1]
        self.aabbmax = [maxx+0.1, maxy+0.1, maxz+0.1]
        
        self.center = vec_add (self.aabbmin , vec_sub (self.aabbmax,self.aabbmin) )
        
        # triangle specific - very useful for lookups in intersection test
        self.edge1 = vec_sub ( self.p[1], self.p[0] )
        self.edge2 = vec_sub ( self.p[2], self.p[0] )
        self.normal = vec_cross ( self.edge1 , self.edge2 )
        self.normal = vec_norm (self.normal)
        
    def ray_intersect (self, ray):
        
        pt = [200,200,200] 
        pvec = vec_cross ( ray.d , self.edge2 )
        det = vec_dot ( self.edge1, pvec )
        
        if det < EPSILON and det > -EPSILON:
            # if det is near zero than plane of tri is same as ray
            return INFINITY, [-1,-1,-1], 0, False
        
        inv_det = 1.0 / det
        
        # calc distance from origin to v0
        tvec = vec_sub ( ray.o , self.p[0])
        
        # calc u and test bounds
        u = vec_dot ( tvec, pvec ) * inv_det
        
        if u < 0.0 or u > 1.0:
            return INFINITY, [-1,-1,-1], 0, False
        
        # calc v and test bounds
        qvec = vec_cross ( tvec, self.edge1 )
        v = vec_dot ( ray.d, qvec ) * inv_det
        if v < 0.0 or (u+v) > 1.0:
            return INFINITY, [-1,-1,-1], 0, False
    
        # intersection - calc t
        t = vec_dot ( self.edge2, qvec ) * inv_det
        if t < 0 : INFINITY, [-1,-1,-1], 0, False
        pt[0] = ray.o[0] + ( t * ray.d[0])
        pt[1] = ray.o[1] + ( t * ray.d[1])
        pt[2] = ray.o[2] + ( t * ray.d[2])
        
        #normal = vec_cross ( edge1 , edge2 )
        #normal = vec_norm (normal)
        
        # return positive normal, depending on which side got hit
        if det < 0:
            if self.flat:
                # return flat shaded normal
                return t,pt,(vec_mul ( self.normal, -1 )), False
            else:
                # calculate normal at intersection point
                N1 = self.vertex_normal[0]
                N2 = self.vertex_normal[1]
                N3 = self.vertex_normal[2]
                NU = vec_add (N1, vec_mul ( (vec_sub (N2,N1)) ,u))
                NV = vec_mul ((vec_sub (N3,N1)), v)
                N = vec_add (NU,NV)
                return t,pt,(vec_mul ( vec_norm(N), -1 )), False
        else:
            return INFINITY, [-1,-1,-1], 0, False
            #return t,pt,self.normal, False
            #return -1, [0,0,0], 0
        
        return t,pt, normal
    
        
    def zrot (self, o, d):
        for v in xrange(3):
            self.p[v] = zrot (self.p[v],o,d)
            self.recalc ()
            
    def yrot (self, o, d):
        for v in xrange(3):
            self.p[v] = yrot (self.p[v],o,d)
            self.recalc ()

   

class Torus(Entity):
    def __init__ (self, center=[0,0,0], r=28, R=40, facing=[0,1,0], mat=Material() ):
        self.center = center
        self.r = r
        self.R = R
        self.facing = facing
        self.mat = mat
        self.aabbmin = vec_add (center, [-R-r-r,-r-r,-R-r-r])
        self.aabbmax = vec_add (center, [R+r+r,r+r,R+r+r])
        
    def ray_intersect (self,ray):
        d = ray.d
        p = ray.o
        r = self.r
        R = self.R
        R2 = R*R
        adota = vec_dot (p,p)
        adotb = vec_dot (p,d)
        
        # define quartic
        K = adota - r*r - R2
        A = 4*adotb
        B = 2*(2*adotb*adotb + K + 2*R2*d[2]*d[2])
        C = 4*(K*adotb + 2*R2*p[2]*d[2])
        D = K*K + 4*R2*(p[2]*p[2] - r*r)
        
        polynomial = scipy.poly1d ( [A,B,C,D] )
        roots = list (polynomial.r)
       
        realroots=[]
        for root in roots:
            realroots.append (float(root) + float(root.imag))

        t = min(realroots)

        if t < 0:
            return INFINITY, [-1,-1,-1], [0,0,1], 0
        else:
            pt = [0,0,0]
            pt[0] = ray.o[0] + ( t * ray.d[0])
            pt[1] = ray.o[1] + ( t * ray.d[1])
            pt[2] = ray.o[2] + ( t * ray.d[2])
            return t, pt, [-1,0,0], 0
        
        
        #d = ray.d
        #p = ray.o
        #r = self.r
        #R = self.R
        #R2 = R*R
        #a = vec_dot (d, d)
        #b = 2 * vec_dot(p, d)
        #y = vec_dot(p, p) - r**2 - R**2
        #
        #a4 = a**2
        #a3 = 2*a*b
        #a2 = b**2 + 2*a*y + 4*(R**2)*(d[2]**2)
        #a1 = 2*b*y + 8*(R**2)*p[2]*d[2]
        #a0 = y**2 + 4*(R**2)*(p[2]**2) - 4*(R**2)*(r**2)
        #
        ## solve the following polynomial for t
        ## a4*t**4 + a3*t**3 + a2*t**2 + a1*t + a0 = 0
        #
        #polynomial = scipy.poly1d ( [a4, a3, a2, a1, a0] )
        #roots = list (polynomial.r)
        #print roots
        #realroots=[]
        #for root in roots:
        #    realroots.append (float(abs(root.imag)))
        #
        #t = min(realroots)
        #
        #if t < 0:
        #    return INFINITY, [-1,-1,-1], 0
        #else:
        #    pt = [0,0,0]
        #    pt[0] = ray.o[0] + ( t * ray.d[0])
        #    pt[1] = ray.o[1] + ( t * ray.d[1])
        #    pt[2] = ray.o[2] + ( t * ray.d[2])
        #    return t, pt, [-1,0,0], 0
       
        
class Sphere(Entity):
    def __init__ (self, center=[0,0,0], radius=10.0):
        Entity.__init__ (self, center, [center[0]-radius, center[1]-radius, center[2]-radius], [center[0]+radius, center[1]+radius, center[2]+radius])
        #self.center=center #deferred to entity class 
        self.radius=radius
        self.rsquared = self.radius**2
        
    def intersect_normal (self,point):
        """ Returns the normal between center of sphere and point"""
        return vec_norm (vec_sub (point,self.center))

    def ray_intersect (self, ray):
        """ Returns the distance to hitpoint as t, and intersect point as pt"""
        #if self.fast_ray_aabb_test ( ray ):
        if 1:
            v = vec_sub (self.center,ray.o) 
            b = vec_dot (v, ray.d)
            disc = (b*b) - vec_dot(v,v) + self.rsquared
            if disc<0:
                return INFINITY, [-1,-1,-1], 0
            d = sqrt(disc)
            t2= b + d
            if t2<0:
                return INFINITY, [-1,-1,-1], 0
            t1 = b - d
            if t1>0:
                t=t1
                inside = False
            else:
                inside = True
                t=t2
            pt=[-1,-1,-1]
            pt[0] = ray.o[0] + t * ray.d[0]
            pt[1] = ray.o[1] + t * ray.d[1]
            pt[2] = ray.o[2] + t * ray.d[2]
            
            # calc normal at intersection to pass back
            normal = vec_sub (pt, self.center)
            normal = vec_norm (normal)
            return t, pt, normal, inside
        else:
            return INFINITY, [-1,-1,-1], 0, False


class Light():
    def __init__(self, pos=[0,0,0], rgb=(1.0,1.0,1.0), power=1.0 ):
        self.pos=pos
        self.num_photons=0
        self.power = 1.0
        self.colour=rgb
        self.multiplier=1.0
        self.samples = 16
        self.radius=20
        self.rsquared = self.radius * self.radius
        
        
    def ray_intersect (self, ray):
        """ Returns the distance to hitpoint as t, and intersect point as pt"""
        #if self.fast_ray_aabb_test ( ray ):
        if 1:
            v = vec_sub (self.center,ray.o) 
            b = vec_dot (v, ray.d)
            disc = (b*b) - vec_dot(v,v) + self.rsquared
            if disc<0:
                return -1, [-1,-1,-1], 0
            d = sqrt(disc)
            t2= b + d
            if t2<0:
                return -1, [-1,-1,-1], 0
            t1 = b - d
            if t1>0:
                t=t1
            else:
                t=t2
            pt=[1,1,1]
            pt[0] = ray.o[0] + t * ray.d[0]
            pt[1] = ray.o[1] + t * ray.d[1]
            pt[2] = ray.o[2] + t * ray.d[2]
            
            # calc normal at intersection to pass back
            normal = vec_sub (pt, self.center)
            normal = vec_norm (normal)
            return t, pt, normal, False
        else:
            return -1, [-1,-1,-1], 0, False




class Box ():
    def __init__ (self, c=[0,0,300], w=20, mat=Material()):
        self.c = c
        self.w = w
        self.mat = mat
        self.objectlist = generate_box (c,w,mat)
        
    def zrot (self, d):
        for tri in self.objectlist:
            tri.zrot ( self.c, d)
            
    def yrot (self, d):
        for tri in self.objectlist:
            tri.yrot ( self.c, d)
            
    def recalc (self):
        self.objectlist = generate_box (self.c,self.w,self.mat)
                

def generate_box ( c=[0,0,300], w=20, mat=Material()):
    faces=[]
    
    # front face
    face = Tri ( [ [c[0]+w, c[1]+w, c[2]-w] , [c[0]-w, c[1]+w, c[2]-w] , [c[0]-w,c[1]-w,c[2]-w] ] )
    face.mat = mat
    faces.append (face)
    face = Tri ( [ [c[0]+w, c[1]+w, c[2]-w] , [c[0]-w, c[1]-w, c[2]-w] , [c[0]+w,c[1]-w,c[2]-w] ] )
    face.mat = mat
    faces.append (face)
    
    # back face
    face = Tri ( [ [c[0]+w, c[1]+w, c[2]+w] , [c[0]-w, c[1]+w, c[2]+w] , [c[0]-w,c[1]-w,c[2]+w] ] )
    face.mat = mat
    face.flip()
    faces.append (face)
    face = Tri ( [ [c[0]+w, c[1]+w, c[2]+w] , [c[0]-w, c[1]-w, c[2]+w] , [c[0]+w,c[1]-w,c[2]+w] ] )
    face.mat = mat
    face.flip()
    faces.append (face)
    
    # left face
    face = Tri ( [ [c[0]-w, c[1]+w, c[2]-w] , [c[0]-w, c[1]+w, c[2]+w] , [c[0]-w,c[1]-w,c[2]+w] ] )
    face.mat = mat
    faces.append (face)
    face = Tri ( [ [c[0]-w, c[1]+w, c[2]-w] , [c[0]-w, c[1]-w, c[2]+w] , [c[0]-w,c[1]-w,c[2]-w] ] )
    face.mat = mat
    faces.append (face)
    
    # right face
    face = Tri ( [ [c[0]+w, c[1]+w, c[2]-w] , [c[0]+w, c[1]+w, c[2]+w] , [c[0]+w,c[1]-w,c[2]+w] ] )
    face.mat = mat
    face.flip()
    faces.append (face)
    face = Tri ( [ [c[0]+w, c[1]+w, c[2]-w] , [c[0]+w, c[1]-w, c[2]+w] , [c[0]+w,c[1]-w,c[2]-w] ] )
    face.mat = mat
    face.flip()
    faces.append (face)
    
    # bottom face
    face = Tri ( [ [c[0]+w, c[1]+w, c[2]-w] , [c[0]-w, c[1]+w, c[2]-w] , [c[0]-w,c[1]+w,c[2]+w] ] )
    face.mat = Material()
    face.mat.colour=[1,1,1]
    face.flip()
    faces.append (face)
    face = Tri ( [ [c[0]+w, c[1]+w, c[2]-w] , [c[0]-w, c[1]+w, c[2]+w] , [c[0]+w,c[1]+w,c[2]+w] ] )
    face.mat = Material()
    face.mat.colour=[1,1,0]
    face.flip()
    faces.append (face)
    
    # top face
    face = Tri ( [ [c[0]+w, c[1]-w, c[2]-w] , [c[0]-w, c[1]-w, c[2]-w] , [c[0]-w,c[1]-w,c[2]+w] ] )
    face.mat = mat
    faces.append (face)
    face = Tri ( [ [c[0]+w, c[1]-w, c[2]-w] , [c[0]-w, c[1]-w, c[2]+w] , [c[0]+w,c[1]-w,c[2]+w] ] )
    face.mat = mat
    faces.append (face)
    
    return faces




    





    
    
        
        
