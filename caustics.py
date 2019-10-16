""" code for making caustics in a ray tracer """

import random
import copy
from vector import *
from surfaces import *
from colour import *
from tracer import MAX_PHOTON_RECURSION

MAX_PHOTON_TREE_DEPTH = 32
MAX_PHOTONS_PER_NODE = 1

def get_photon_cube_size (photon_list):
    # for a given list of objects, return the AABB of them all combined
    # nb. returned aabb is a tuple of aabbmin,aabbmax
    min = [INFINITY,INFINITY,INFINITY]
    max = [-INFINITY,-INFINITY,-INFINITY]
    for photon in photon_list:
        if photon.position[0] < min[0] : min[0]=photon.position[0]
        if photon.position[1] < min[1] : min[1]=photon.position[1]
        if photon.position[2] < min[2] : min[2]=photon.position[2]
        if photon.position[0] > max[0] : max[0]=photon.position[0]
        if photon.position[1] > max[1] : max[1]=photon.position[1]
        if photon.position[2] > max[2] : max[2]=photon.position[2]    
    return min, max

class Photon ():
    def __init__(self,position=[0,0,0],power=-5, d=[0,0,0],flag=5, c=[1,1,1]):
        self.position=position
        self.colour=c
        self.power=power
        self.d = d

class Photon_Tree_Node ():
    def __init__ (self):
        self.left = None
        self.right = None
        self.aabbmin = [-INFINITY, -INFINITY, -INFINITY]
        self.aabbmax = [INFINITY,INFINITY,INFINITY]
        self.photon = Photon()
        self.level = -1
        self.split_axis = -1
    
    def ray_aabb_test (self, ray):
        tnear = -INFINITY
        tfar = INFINITY
        for p in xrange(3):
            # for each plane of bounding box
            if ray.d[p] == 0.0: #if direction of ray is parallel to the axis
                if ray.o[p] < self.aabbmin[p] or ray.o[p] > self.aabbmax[p]:
                    return (False)
            else:
                t1 = (self.aabbmin[p] - ray.o[p]) / ray.d[p]
                t2 = (self.aabbmax[p] - ray.o[p]) / ray.d[p]
                if t1 > t2: #swap t1 and t2
                    temp = t1
                    t1 = t2
                    t2 = temp
                if t1 > tnear: tnear = t1
                if t2 < tfar: tfar = t2
                if tnear > tfar: return (False)
                if tfar < 0: return (False)
        return (True) #Intersection with AABB occurs
    
    def get_photon_list (self):
        """ return list of photons underneath this node"""
        templist=[]
        templist.append ( self.photon )
        if self.left is not None:
            templist = templist + self.left.get_photon_list ()
        if self.right is not None:
            templist = templist + self.right.get_photon_list ()
        return templist
    
    def point_near_node (self,p,d):
        """ for a given point p, check to see if it is within distance d of bounding box """
        if p[0] < self.aabbmin[0]-d : return False
        if p[1] < self.aabbmin[1]-d : return False
        if p[2] < self.aabbmin[2]-d : return False
        if p[0] > self.aabbmax[0]+d : return False
        if p[1] > self.aabbmax[1]+d : return False
        if p[2] > self.aabbmax[2]+d : return False
        return True
      
    def get_nearest_photons (self, p, d):
        """ gets points in photon map within radius d of point p """
        nearest_photons=[]
        
        if self.left is not None:
            if self.left.point_near_node (p,d):
                nearest_photons = nearest_photons + self.left.get_nearest_photons (p, d)
            
        if self.right is not None:
            if self.right.point_near_node (p,d):
                nearest_photons = nearest_photons + self.right.get_nearest_photons (p, d)
                        
        # calc square distance from current photon to point p
        dist = vec_distance (self.photon.position, p)
        if dist < d:
            # photon is close enough
            nearest_photons.append (self.photon)
            return nearest_photons
        return nearest_photons
    
    def populate (self, photon_list):
        if len (photon_list) > 1 and self.level < MAX_PHOTON_TREE_DEPTH:
            # need to split tree more
            # get the longest axis of cube surrounding photons
            dx=self.aabbmax[0] - self.aabbmin[0]
            dy=self.aabbmax[1] - self.aabbmin[1]
            dz=self.aabbmax[2] - self.aabbmin[2]
            axis = -1
            if dx > dy and dx > dz: axis = 0
            if dy > dx and dy > dz: axis = 1
            if dz > dx and dz > dy: axis = 2
            if axis == -1 : axis = 0
            self.split_axis = copy.deepcopy(axis)
            # find the median - sort list by axis, then pick half way point
            median = len(photon_list)/2 
            photon_list.sort(cmp=lambda x,y:cmp(x.position[axis],y.position[axis]))
            splitpos = copy.deepcopy (photon_list[median].position[axis])
            self.photon = copy.deepcopy(photon_list[median]) # put photon at median into current node
            photon_list.remove (photon_list[median])
            # get list of photons for each potential child node
            leftlist=[]
            rightlist=[]
            for photon in photon_list:
                if photon.position[axis] < splitpos:
                    leftlist.append (copy.deepcopy(photon))
                if photon.position[axis] > splitpos:
                    rightlist.append (copy.deepcopy(photon))
            if len (leftlist) > 0:
                self.left = Photon_Tree_Node ()
                self.left.level = self.level + 1
                self.left.aabbmin, self.left.aabbmax = get_photon_cube_size (leftlist)
                self.left.populate (leftlist)
            if len (rightlist) > 0:
                self.right = Photon_Tree_Node ()
                self.right.level = self.level + 1
                self.right.aabbmin, self.right.aabbmax = get_photon_cube_size (rightlist)
                self.right.populate (rightlist)
        else:
            # stick single photon left into current node
            self.photon = copy.deepcopy (photon_list[0])
        

class Photon_Tree ():
    def __init__ (self, photon_list=[], level=0):
        # create a root node
        print "Generating photon tree from list of photons...."
        self.node = Photon_Tree_Node ()
        self.level = level
        self.node.aabbmin, self.node.aabbmax = get_photon_cube_size (photon_list)
        self.node.populate (photon_list)
        print "Photon tree complete."
    
        

def trace_photon (scene_tree, photon, recursion_level=1):
    """ trace a photon defined by photon_ray into scene tree """
    """ return a list of photons - and update photon lists of surfaces intersected """
    photon_ray = Ray (photon.position, photon.d)
    object_list = scene_tree.get_object_list (photon_ray)
    photon_list = []
    # get nearest object
    nearest_object = None
    nearest_delta = INFINITY
    nearest_intersection = [0,0,INFINITY]
    for entity in object_list:
        result = entity.ray_intersect(photon_ray )
        delta = result[0]
        if delta > 0.0 and delta < INFINITY:
            #take note of distance
            if delta < nearest_delta:
                #current object is closest yet
                nearest_delta=delta
                nearest_object=entity
                nearest_intersection = copy.deepcopy(result[1])
                norm = copy.deepcopy(result[2])
                inside = copy.deepcopy(result[3])
                
    if nearest_object is not None and recursion_level < MAX_PHOTON_RECURSION:
        # photon hit an object
        
        # do caustic re-generation first
        if nearest_object.mat.send_caustics == True and nearest_object.mat.spec_level > 0:
            
            if nearest_object.mat.reflect > 0:
                # generate reflective caustic
                fudged_hitpoint = vec_add (nearest_intersection, vec_div (norm,1000000.0) )
                c = vec_dot (norm, photon.d)
                c = -c
                reflection_dir = vec_add (photon.d, vec_mul (vec_mul(norm,c),2) )
                reflection_dir = vec_norm (reflection_dir) # TO DO - remove, may not be necessary
                # now create new photon at intersection heading in new direction, reducing power correctly
                new_colour = colour_alpha (photon.colour, nearest_object.mat.colour, nearest_object.mat.reflect)
                new_photon = Photon (fudged_hitpoint, photon.power * nearest_object.mat.reflect, reflection_dir, 0, new_colour)
                photon_list = photon_list + trace_photon (scene_tree, new_photon, recursion_level + 1)
                
            if nearest_object.mat.refract > 0:
                # generate reflective caustic
                test = vec_dot (norm, vec_norm(photon.d))
                if test>0:
                    # then inside object, reverse n1 and n2
                    NR = vec_mul (norm, -1)
                    c1 = -vec_dot (NR, vec_norm(photon.d))
                    n2 = 1.0
                    n1 = nearest_object.mat.ior
                else:
                    NR = vec_norm(norm)
                    c1 = -vec_dot (NR, vec_norm(photon.d))
                    n2 = nearest_object.mat.ior
                    n1 = 1.0 
                n = n1 / n2
                c2 = sqrt( 1 - (n*n) * (1 - (c1*c1)) )
                a = vec_mul (photon.d,n)
                b = vec_mul (NR, n*c1-c2)
                Rr = vec_add (a,b)
                refraction_dir = vec_norm (Rr) # TODO - normalisation may not be needed here
                fudged_hitpoint = vec_add (nearest_intersection, vec_div (Rr,1000000.0))
                # now create new photon at intersection heading in new direction, reducing power correctly
                new_colour = colour_alpha (photon.colour, nearest_object.mat.colour, nearest_object.mat.refract)
                new_photon = Photon (fudged_hitpoint, photon.power * nearest_object.mat.refract, refraction_dir, 0, new_colour)
                photon_list = photon_list + trace_photon (scene_tree, new_photon, recursion_level + 1)
                
        # now do caustic receiving
        if nearest_object.mat.receive_caustics == True and recursion_level > 1:
            # photon has hit a surface to be displayed
            new_photon = Photon (nearest_intersection, photon.power, norm, 0, photon.colour)
            #nearest_object.photon_list.append (new_photon)
            photon_list.append (new_photon)
    return photon_list   

def generate_photon_list ( scene_tree, lights ):
    """ create photon kd tree """
    print "Generating photon list...."
    photon_tree = []
    count = 0
    for light in lights:
        if light.num_photons > 0:
            for ne in xrange (light.num_photons):
                #print "Photon:",ne,"/",light.num_photons
                photon_ok = False
                while photon_ok == False:
                    # TODO: get list of surfaces in scene which WANT caustics and make sure random rays go in that direction
                    # then I can get rid of ray checking a few lines down
                    x = -1 + random.random()*2
                    y = -1 + random.random()*2
                    z = -1 + random.random()*2
                    if (x**2 + y**2 + z**2 < 1):
                        ray_direction = vec_norm ([x,y,z])
                        #photon_ray = Ray (light.pos, ray_direction)
                        # check ray might intersect with a specular object
                        #object_list = scene_tree.get_object_list (photon_ray)
                        #for item in object_list:
                        #    if item.mat.spec_level > 0 and item.mat.send_caustics == True:
                        #        photon_ok = True
                        #        break
                        
                        #debug
                        photon_ok = True
                        break
                # assume point light for now
                photon = Photon (light.pos, light.power, ray_direction, 0, light.colour)
                traced_photon = trace_photon (scene_tree, photon)
                photon_tree = photon_tree + traced_photon
                count+=1
                if count%5000==0:print count,"/",light.num_photons
                    
    print "Photon list complete. Total of ",len (photon_tree),"photons visible out of",light.num_photons,"fired."
    return photon_tree
    
