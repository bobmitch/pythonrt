
# KD TREE python module
# Version 0.1

MAX_TREE_DEPTH = 17 #100,000
MAX_NODE_SIZE = 3
MAX_LEVEL = 0

from surfaces import INFINITY
from surfaces import EPSILON
from surfaces import Plane
import copy


def get_extents (objectlist):
    # for a given list of objects, return the AABB of them all combined
    # nb. returned aabb is a tuple of aabbmin,aabbmax
    min = [INFINITY,INFINITY,INFINITY]
    max = [-INFINITY,-INFINITY,-INFINITY]
    for item in objectlist:
        if item.aabbmin[0] < min[0]: min[0]=item.aabbmin[0]
        if item.aabbmax[0] > max[0]: max[0]=item.aabbmax[0]
        if item.aabbmin[1] < min[1]: min[1]=item.aabbmin[1]
        if item.aabbmax[1] > max[1]: max[1]=item.aabbmax[1]
        if item.aabbmin[2] < min[2]: min[2]=item.aabbmin[2]
        if item.aabbmax[2] > max[2]: max[2]=item.aabbmax[2]
        
    return min, max


def aabb_overlap (aabbmin1, aabbmax1, aabbmin2, aabbmax2):
    for axis in xrange(3):
        if aabbmin1[axis] >= aabbmax2[axis]: return False
        if aabbmax1[axis] <= aabbmin2[axis]: return False
    return True


class Node():
    def __init__(self):
        self.aabbmin = None
        self.aabbmin = None
        self.isleaf = False
        self.left = None
        self.right = None
        self.level = None
        self.objectlist = []
        
 
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
    
    def get_nodes_traversed (self,ray):
        count = 1
        if ray.o[0] > self.aabbmin[0] and ray.o[1] > self.aabbmin[1] and ray.o[2] > self.aabbmin[2] and ray.o[0] < self.aabbmax[0] and ray.o[1] < self.aabbmax[1] and ray.o[2] < self.aabbmax[2]:
            if self.left is None and self.right is None: return 1
            if self.left is not None:
                count = count + self.left.get_nodes_traversed (ray)
            if self.right is not None:
                count = count + self.right.get_nodes_traversed (ray)
            return count
        # ok, ray not within aabb, do proper check
        if self.ray_aabb_test (ray):
            if self.left is None and self.right is None: return count
            if self.left is not None:
                count = count + self.left.get_nodes_traversed (ray)
            if self.right is not None:
                count = count + self.right.get_nodes_traversed (ray)
            return count
        else:
            return count
    
    def get_object_list (self, ray, infinite_objects=[]):
        # first, check if ray origin is within aabb
        templist = infinite_objects
        if ray.o[0] > self.aabbmin[0] and ray.o[1] > self.aabbmin[1] and ray.o[2] > self.aabbmin[2] and ray.o[0] < self.aabbmax[0] and ray.o[1] < self.aabbmax[1] and ray.o[2] < self.aabbmax[2]:
            # ray origin within aabb - automatic hit
            if self.left is None and self.right is None: return templist + self.objectlist
            # else, traverse child nodes
            if self.left is not None:
                templist = templist + self.left.get_object_list (ray)
            if self.right is not None:
                templist = templist + self.right.get_object_list (ray)
            return templist
        # ok, ray not within aabb, do proper check
        if self.ray_aabb_test (ray):
            # ray hits current node, and no more nodes left, so return objects 
            if self.left is None and self.right is None: return templist + self.objectlist
            # else, traverse child nodes
            if self.left is not None:
                templist = templist + self.left.get_object_list (ray)
            if self.right is not None:
                templist = templist + self.right.get_object_list (ray)
            return templist
        else:
            return templist
        
    
    def populate_node_grid (self, scene):
        # SPLITS AT PARENT SPATIAL MEDIAN - WORKS OK
        if self.level > MAX_TREE_DEPTH or len(scene) <= MAX_NODE_SIZE:
            self.objectlist = scene
            self.isleaf = True
            return
    
        # rotate split axis according to level in tree
        axis = self.level % 3
        splitpos = self.aabbmin[axis] + ((self.aabbmax[axis] - self.aabbmin[axis])/2.0)
        
        # generate potential child node AABBs
        left_aabbmin = copy.deepcopy (self.aabbmin)
        left_aabbmax = copy.deepcopy (self.aabbmax)
        right_aabbmin = copy.deepcopy (self.aabbmin)
        right_aabbmax = copy.deepcopy (self.aabbmax)
        left_aabbmax[axis] = splitpos - 0.01
        right_aabbmin[axis] = splitpos
        
        # generate list of overlapping objects for each potential child node
        leftlist=[]
        rightlist=[]
        for item in scene:
            if aabb_overlap (left_aabbmin, left_aabbmax, item.aabbmin, item.aabbmax):
                leftlist.append (item)
        for item in scene:
            if aabb_overlap (right_aabbmin, right_aabbmax, item.aabbmin, item.aabbmax):
                rightlist.append (item)
        
        if len(leftlist) > 0:
            self.left = Node()
            self.left.level = self.level + 1
            self.left.aabbmin = left_aabbmin
            self.left.aabbmax = left_aabbmax
            self.left.populate_node_grid (leftlist)
        if len(rightlist) > 0:
            self.right = Node()
            self.right.level = self.level + 1
            self.right.aabbmin = right_aabbmin
            self.right.aabbmax = right_aabbmax
            self.right.populate_node_grid (rightlist)
    
        
    def populate_node_naive (self, scene):
        # USES OBJECT MEDIAN AS SPLITPOINT - DOESNT WORK RIGHT AT THE MOMENT
        if self.level > MAX_TREE_DEPTH or len(scene) <= MAX_NODE_SIZE:
            self.objectlist = scene
            self.isleaf = True
            return
    
        # rotate split axis according to level in tree
        axis = self.level % 3
        # sort scene
        median = len(scene)/2 # choose median
        scene.sort(cmp=lambda x,y:cmp(x.center[axis],y.center[axis]))
        splitpos = copy.deepcopy (scene[median].aabbmax[axis])
        
        # generate potential child node AABBs
        left_aabbmin = copy.deepcopy (self.aabbmin)
        left_aabbmax = copy.deepcopy (self.aabbmax)
        right_aabbmin = copy.deepcopy (self.aabbmin)
        right_aabbmax = copy.deepcopy (self.aabbmax)
        left_aabbmax[axis] = splitpos
        right_aabbmin[axis] = splitpos + EPSILON
        
        # generate list of overlapping objects for each potential child node
        leftlist=[]
        rightlist=[]
        for item in scene:
            if aabb_overlap (left_aabbmin, left_aabbmax, item.aabbmin, item.aabbmax):
                leftlist.append (item)
        for item in scene:
            if aabb_overlap (right_aabbmin, right_aabbmax, item.aabbmin, item.aabbmax):
                rightlist.append (item)
        
        if len(leftlist) > 0:
            self.left = Node()
            self.left.level = self.level + 1
            self.left.aabbmin = left_aabbmin
            self.left.aabbmax = left_aabbmax
            self.left.populate_node_naive (leftlist)
        if len(rightlist) > 0:
            self.right = Node()
            self.right.level = self.level + 1
            self.right.aabbmin = right_aabbmin
            self.right.aabbmax = right_aabbmax
            self.right.populate_node_naive (rightlist)
                
        

        
        
    
    def display_node (self):
        print "Level:",self.level
        print "AABBmin:",self.aabbmin
        print "AABBmax:",self.aabbmax
        print "Left:",type(self.left)
        print "Right:",type(self.right)
        print "Contents:",len(self.objectlist)
        if self.left is not None: self.left.display_node ()
        if self.right is not None: self.right.display_node ()
        
class Tree ():
    # constructor for root node of kd tree
    # takes a list of objects as defined in the surfaces.py module
    def __init__(self, scene=[]):
        # process infinite planes - place in special list, remove from scene
        # and dont add to kdtree
        self.infinite_objects = []
        for item in scene:
            if hasattr (item,'infinite'):
                self.infinite_objects.append (item)
        for item in self.infinite_objects:
            scene.remove (item)
        self.node = Node()
        self.node.level = 0
        self.node.aabbmin, self.node.aabbmax = get_extents (scene)
        #self.node.populate_node_naive (scene)
        self.node.populate_node_grid (scene)
        
    def display_tree (self):
        self.node.display_node ()
        
    def get_object_list (self, ray):
        # for a given ray return the object(s) potentially hit in the tree - including any infinite objects like planes
        return self.node.get_object_list (ray, self.infinite_objects)
    
    def get_nodes_traversed (self, ray):
        return self.node.get_nodes_traversed (ray)
        
    
