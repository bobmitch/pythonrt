
from surfaces import *
from kd_tree import *
from vector import *
from colour import *
from caustics import *
from objloader import *
import pygame
from pygame.locals import *
import copy
import psyco
import sys
import random
from rotate import *
from math import tan
from math import sqrt
import threading
psyco.full()

SCREENW=400.0
SCREENH=300.0
FOV=65

AIR_IOR = 1.0
AMBIENT = (0.1,0.1,0.1)
BACKGROUND=( 0.95, 0.92, 1.0 )
SKY_HORIZON = (0.89,0.92,1)
SKY_TOP = (0.65,0.63,0.9)
FLOOR = (0,0,0)
#SKY_HORIZON = (0,0,0)
#SKY_TOP = (0,0,0)
#FLOOR = (0,0,0)

PIXELSIZE=1
SCANLINES=1

SHADOWS=True
REFLECTIONS=1
REFRACTIONS=1
CAUSTICS=0
AA=0
MAX_RECURSION_LEVEL=3
MAX_PHOTON_RECURSION=5

CAUSTIC_RADIUS = 10.0
BLURRY_SAMPLES=32
SOFT_SHADOW_SAMPLES=32
SHOW_TREE=0
SHOW_SCENE=1



#rect=pygame.rect.Rect ( 10,10,2,2)
#clock=pygame.time.Clock()

def unique(s):
     """Return a list of the elements in s, but without duplicates.

     For example, unique([1,2,3,1,2,3]) is some permutation of [1,2,3],
     unique("abcabc") some permutation of ["a", "b", "c"], and
     unique(([1, 2], [2, 3], [1, 2])) some permutation of
     [[2, 3], [1, 2]].

     For best speed, all sequence elements should be hashable.  Then
     unique() will usually work in linear time.

     If not possible, the sequence elements should enjoy a total
     ordering, and if list(s).sort() doesn't raise TypeError it's
     assumed that they do enjoy a total ordering.  Then unique() will
     usually work in O(N*log2(N)) time.

     If that's not possible either, the sequence elements must support
     equality-testing.  Then unique() will usually work in quadratic
     time.
     """

     n = len(s)
     if n == 0:
         return []

     # Try using a dict first, as that's the fastest and will usually
     # work.  If it doesn't work, it will usually fail quickly, so it
     # usually doesn't cost much to *try* it.  It requires that all the
     # sequence elements be hashable, and support equality comparison.
     u = {}
     try:
         for x in s:
             u[x] = 1
     except TypeError:
         del u  # move on to the next method
     else:
         return u.keys()

     # We can't hash all the elements.  Second fastest is to sort,
     # which brings the equal elements together; then duplicates are
     # easy to weed out in a single pass.
     # NOTE:  Python's list.sort() was designed to be efficient in the
     # presence of many duplicate elements.  This isn't true of all
     # sort functions in all languages or libraries, so this approach
     # is more effective in Python than it may be elsewhere.
     try:
         t = list(s)
         t.sort()
     except TypeError:
         del t  # move on to the next method
     else:
         assert n > 0
         last = t[0]
         lasti = i = 1
         while i < n:
             if t[i] != last:
                 t[lasti] = last = t[i]
                 lasti += 1
             i += 1
         return t[:lasti]

     # Brute force is all that's left.
     u = []
     for x in s:
         if x not in u:
             u.append(x)
     return u



def check_for_exit():
    for event in pygame.event.get():
        if event.type==KEYDOWN:
            if event.key==K_ESCAPE:
                pygame.quit()
                sys.exit(0)
            


def get_colour ( ray, tree, lights, recursion_level=0, photon_tree=None ):
    """ for a given ray, and a given scene, and given lights, returns the colour calculated
        as viewed by the ray.  if no colour determined, returns -1 tupled"""
    nearestdelta = INFINITY
    nearesthitpoint = [0,0,INFINITY]
    count=0
    
    scene = tree.get_object_list (ray)
    #print "Get colour for ray returned",len(scene),"objects."
    for entity in scene:
        result = entity.ray_intersect(ray)
        delta = result[0]
        #print delta
        if delta > 0.0 and delta < INFINITY:
            #take note of distance
            if delta < nearestdelta:
                #current object is closest yet
                nearestdelta=delta
                nearestobject=entity
                nearesthitpoint = copy.deepcopy(result[1])
                norm = copy.deepcopy(result[2])
                inside = copy.deepcopy(result[3])
                
                
    if nearestdelta==INFINITY or nearestdelta<0:
        #no objects hit by ray
        # CALCULATE BACKGROUND GRADIENT
        angle = -ray.d[1]
        if angle > 1 : print "NOT NORMAL"
        if angle<0:
            angle=-angle
            mix1 = colour_scalar_mul ( FLOOR , angle )
            mix2 = colour_scalar_mul ( SKY_HORIZON, 1.0 - angle )
            skycolour = colour_add ( mix1, mix2 )
            return skycolour 
        else:
            mix1 = colour_scalar_mul ( SKY_TOP , angle )
            mix2 = colour_scalar_mul ( SKY_HORIZON, 1.0 - angle )
            skycolour = colour_add ( mix1, mix2 )
            return skycolour    
        return BACKGROUND      
    else:
        # ok, ray hit something ('nearestobject')
        final_colour = [0,0,0]
        
        # reset reflection and refraction
        refraction = [-1,-1,-1]
        reflection = [-1,-1,-1]
        
        if nearestobject.mat.refract < 1 and nearestobject.mat.reflect < 1:
            
            final_colour = colour_mul ( AMBIENT , nearestobject.mat.colour )  
            light_amount = 0.0 # assume full light
            for light in lights: 
                #create normalised vector from hitpoint to light
                hittolight = vec_sub (light.pos, nearesthitpoint)                    
                hittolight = vec_norm (hittolight)
                light_power = light.multiplier
                blocked = False
                shadow_ray = Ray (nearesthitpoint, hittolight)
                potential_shadow_objects=[]
                if light.radius==0.0:
                    # only one single possible group of shadow objects for point light with no radius
                    light.samples=1
                    potential_shadow_objects = tree.get_object_list (shadow_ray)
                else:
                    # figure out other objects depending on size of light - make 6 other rays at light extents
                    potential_shadow_objects = potential_shadow_objects + tree.get_object_list (shadow_ray) 
                    #original ray, do 6 others
                    test_pos = [light.pos[0]-light.radius, light.pos[1], light.pos[2]]
                    test_to_light = vec_norm (vec_sub (test_pos, nearesthitpoint))
                    test_ray = Ray (nearesthitpoint, test_to_light)
                    potential_shadow_objects = potential_shadow_objects + tree.get_object_list (test_ray)
                    test_pos = [light.pos[0]+light.radius, light.pos[1], light.pos[2]]
                    test_to_light = vec_norm (vec_sub (test_pos, nearesthitpoint))
                    test_ray = Ray (nearesthitpoint, test_to_light)
                    potential_shadow_objects = potential_shadow_objects + tree.get_object_list (test_ray)
                    test_pos = [light.pos[0], light.pos[1]+light.radius, light.pos[2]]
                    test_to_light = vec_norm (vec_sub (test_pos, nearesthitpoint))
                    test_ray = Ray (nearesthitpoint, test_to_light)
                    potential_shadow_objects = potential_shadow_objects + tree.get_object_list (test_ray)
                    test_pos = [light.pos[0], light.pos[1]-light.radius, light.pos[2]]
                    test_to_light = vec_norm (vec_sub (test_pos, nearesthitpoint))
                    test_ray = Ray (nearesthitpoint, test_to_light)
                    potential_shadow_objects = potential_shadow_objects + tree.get_object_list (test_ray)
                    test_pos = [light.pos[0], light.pos[1], light.pos[2]+light.radius]
                    test_to_light = vec_norm (vec_sub (test_pos, nearesthitpoint))
                    test_ray = Ray (nearesthitpoint, test_to_light)
                    potential_shadow_objects = potential_shadow_objects + tree.get_object_list (test_ray)
                    test_pos = [light.pos[0], light.pos[1], light.pos[2]-light.radius]
                    test_to_light = vec_norm (vec_sub (test_pos, nearesthitpoint))
                    test_ray = Ray (nearesthitpoint, test_to_light)
                    potential_shadow_objects = potential_shadow_objects + tree.get_object_list (test_ray)

                    potential_shadow_objects = unique (potential_shadow_objects)
                
                light_amount = 0.0 # assume full light
                for shadow_object in potential_shadow_objects:
                    if shadow_object != nearestobject:
                        # do light.samples number of samples
                        sample_list=[0]*light.samples
                        for s in xrange(light.samples):
                            dx = (-1 + random.random()*2) * (light.radius)
                            dy = (-1 + random.random()*2) * (light.radius)
                            dz = (-1 + random.random()*2) * (light.radius)
                            sample_pos = [light.pos[0]+dx,light.pos[1]+dy,light.pos[2]+dz]
                            sample_to_hitpoint = vec_sub (sample_pos, nearesthitpoint)
                            sample_to_hitpoint = vec_norm (sample_to_hitpoint)
                            sample_ray = Ray (nearesthitpoint, sample_to_hitpoint)
                            result = shadow_object.ray_intersect (sample_ray)
                            sample_list[s] = result[0]
                            if len(sample_list)==5:
                                # early out - if first so many rays are blocked, assume most will be :)
                                if sum(sample_list)<INFINITY:
                                    blocked = True
                                    break
                        sample_list.sort()
                        #print sample_list
                        if blocked==True or (sample_list[0] > 0 and sample_list[len(sample_list)-1] < INFINITY):
                            blocked = True
                            break

                        # reach this point if there is some soft shadow required
                        samples_missed = 0
                        for s in sample_list:
                            if s==INFINITY: samples_missed+=1
                        m = 1.0 / len(sample_list)
                        light_amount = light_amount + (1.0 - (samples_missed * m))
                        
                           

                        #result = shadow_object.ray_intersect (shadow_ray )
                        #if result[0] > 0 and result[0]!=INFINITY:
                        #    blocked = True
                        #    break
                        #else:
                        #    # do nothing
                        #    pass
                                         
                if blocked == False: # light has passed through some semi-transparent object(s)
                    # intersection point can see current light
                    light_normal = vec_dot (norm, hittolight)
                    if light_normal < 0 : light_normal = 0
                    if light_normal > 0:
                        #light_colour = light.colour
                        diffuse_colour = colour_scalar_mul (nearestobject.mat.colour, light_normal)
                        diffuse_colour = colour_mul (diffuse_colour, light.colour)
                        if light_amount < 0:light_amount=0.0
                        if light_amount > 1:light_amount=1.0
                        diffuse_colour = vec_mul ( diffuse_colour, 1.0 - light_amount )
                        final_colour = colour_add ( final_colour, diffuse_colour) 
                        #do specular phong - intensity = diffuse * (L.N) + specular * (V.R)n
                        if nearestobject.mat.spec_level > 0:
                            fudged_hitpoint = vec_add (nearesthitpoint, vec_div (norm,1000000.0) )
                            c = vec_dot (norm, ray.d)
                            c = -c
                            reflection_dir = vec_add (ray.d, vec_mul (vec_mul(norm,c),2) )
                            #reflection_ray = Ray (fudged_hitpoint, reflection_dir)
                            V = ray.d
                            R = reflection_dir
                            L = shadow_ray.d
                            vrdot = vec_dot ( V, R )
                            dot = vec_dot ( L, R )
                            if dot>0:
                                dot = dot ** nearestobject.mat.spec_size
                                spec_strength = dot * nearestobject.mat.spec_level * light_normal
                                spec = colour_scalar_mul ( light.colour, spec_strength )
                            else:
                                spec = (0,0,0)
                            final_colour = colour_add ( final_colour, spec )
                
        
        if REFRACTIONS and nearestobject.mat.refract > 0.0 and nearestobject.mat.reflect < 1 and recursion_level < MAX_RECURSION_LEVEL:
            test = vec_dot (norm, vec_norm(ray.d))
            if test>0:
                # then inside object, reverse n1 and n2
                NR = vec_mul (norm, -1)
                c1 = -vec_dot (NR, vec_norm(ray.d))
                n2 = 1.0
                n1 = nearestobject.mat.ior
            else:
                NR = vec_norm(norm)
                c1 = -vec_dot (NR, vec_norm(ray.d))
                n2 = nearestobject.mat.ior
                n1 = 1.0 
            n = n1 / n2
            c2 = sqrt( 1 - (n*n) * (1 - (c1*c1)) )
            a = vec_mul (ray.d,n)
            b = vec_mul (NR, n*c1-c2)
            Rr = vec_add (a,b)
            Rr = vec_norm (Rr)
            fudged_hitpoint = vec_add (nearesthitpoint, vec_div (Rr,1000000.0))
            refraction_ray = Ray (fudged_hitpoint, Rr)
            refraction = get_colour ( refraction_ray, tree, lights, recursion_level+1)
            refraction = colour_scalar_mul (refraction, nearestobject.mat.refract)
            final_colour = colour_add ( colour_scalar_mul (final_colour, 1-nearestobject.mat.refract) , refraction)
                      
        if nearestobject.mat.reflect > 0 and nearestobject.mat.refract < 1 and REFLECTIONS and recursion_level < MAX_RECURSION_LEVEL:
            # calulcate reflection vector using: Rl = V + (2 * N * c )
            #fudge the hitpoint out by a smidgeon along the normal so reflection ray doesnt hit the object itself
            fudged_hitpoint = vec_add (nearesthitpoint, vec_div (norm,1000000.0) )
            #fudged_hitpoint = nearesthitpoint
            c = vec_dot (norm, ray.d)
            c = -c
            reflection_dir = vec_add (ray.d, vec_mul (vec_mul(norm,c),2) )
            reflection_ray = Ray (fudged_hitpoint, reflection_dir)
            if recursion_level < MAX_RECURSION_LEVEL and REFLECTIONS:
                if nearestobject.mat.reflection_blur > 0:
                    # TODO: Implement blurry reflections
                    pass
                else:
                    reflection = get_colour ( reflection_ray, tree, lights, recursion_level+1)
                    reflection = colour_scalar_mul (reflection, nearestobject.mat.reflect)
                    final_colour = colour_add ( colour_scalar_mul (final_colour, 1-nearestobject.mat.reflect) , reflection)
        
        # process caustics
        if photon_tree is not None and CAUSTICS:
            photon_list = photon_tree.node.get_nearest_photons (nearesthitpoint, CAUSTIC_RADIUS)
            #photon_list = photon_tree.node.get_photon_list ()
            num_photons = len(photon_list)
            if num_photons > 0:
                #DEBUG
                #print num_photons,
                # get an average of all the photons within CAUSTIC_RADIUS of hitpoint (discard those facing the wrong way)
                avg_normal = 0
                avg_power = 0
                avg_colour = [0,0,0]
                total_nearby_photons = 0
                for photon in photon_list:
                    photon_normal = vec_dot (photon.d, norm)
                    if photon_normal > 0:
                        total_nearby_photons += 1
                        avg_normal = avg_normal + photon_normal
                        avg_power = avg_power + photon.power
                        avg_colour = vec_add (avg_colour, photon.colour)
                if total_nearby_photons > 0:
                    avg_colour = vec_div (avg_colour, total_nearby_photons)
                    avg_normal = (avg_normal / total_nearby_photons)
                    avg_power = avg_power / (num_photons/total_nearby_photons) #remove this line maybe
                    multiplier = (avg_normal * avg_power) / CAUSTIC_RADIUS * 30 
                    colour_multiplier = colour_scalar_mul (avg_colour, 0.01)
                    photon_colour = colour_scalar_mul (colour_multiplier,multiplier)
                    final_colour = colour_add (final_colour, photon_colour)
                    #debug
                    #final_colour = [1,1,1]
            
        if final_colour == [0,0,0] :
            #print "FINAL COLOUR IS 0,0,0"
            return AMBIENT 
        return final_colour


def rgb_clamp (r,g,b):
    if r>255:r=255
    if g>255:g=255
    if b>255:b=255
    if r<0:r=0
    if g<0:g=0
    if b<0:b=0
    return r,g,b



def generate_pyramid ( c=[0,0,300], size=20, h=40):
    faces=[]
    reflect=0.2
    face = Tri ( [ [c[0]-size, c[1], c[2]-size] , [c[0], c[1]-h, c[2]], [c[0]+size,c[1],c[2]-size] ] ) #front
    face.mat.reflect=reflect
    faces.append (face)
    face = Tri ( [ [c[0]-size, c[1], c[2]+size] , [c[0], c[1]-h, c[2]], [c[0]+size,c[1],c[2]+size] ] ) #back
    face.mat.reflect=reflect
    faces.append (face)
    face = Tri ( [ [c[0]-size, c[1], c[2]+size] , [c[0], c[1]-h, c[2]], [c[0]-size,c[1],c[2]-size] ] ) #left
    face.mat.reflect=reflect
    faces.append (face)
    face = Tri ( [ [c[0]+size, c[1], c[2]-size] , [c[0], c[1]-h, c[2]], [c[0]+size,c[1],c[2]+size] ] ) #right
    face.mat.reflect=reflect
    faces.append (face)
    
    return faces


def generate_sphere_flake ( pos=[-80,80,400], radius=70.0 , maxdepth=3, curdepth = 1 , flakeset=[]):
    """ takes position, radius and maxdepth, and returns a set of spheres """

    # make sphere at current location
    sphere = Sphere ( pos , radius )
    sphere.mat.reflect=0.2
    r=random.random()
    g=random.random()
    b=random.random()
    sphere.mat.spec_level = random.random()
    sphere.mat.spec_size = random.randrange(1,55)
    sphere.mat.colour=(1,1,0)
    flakeset.append ( sphere )
    if curdepth >= maxdepth:
        # at maximum depth, return set of spheres
        return flakeset
    else:
        # recursively create more flakes - only does vertical at the moment
        scale = 1.0 / 1.8
        newradius = radius * scale
        
        newpos = copy.deepcopy (pos)
        newpos[1] = pos[1] - radius - (newradius)
        flakeset =  ( generate_sphere_flake ( newpos, newradius, maxdepth, curdepth+1) )
        
        newpos = copy.deepcopy (pos)
        newpos[0] = pos[0] - radius - (newradius)
        flakeset =  ( generate_sphere_flake ( newpos, newradius, maxdepth, curdepth+1) )
        
        newpos = copy.deepcopy (pos)
        newpos[0] = pos[0] + radius + (newradius)
        flakeset =  ( generate_sphere_flake ( newpos, newradius, maxdepth, curdepth+1) )
        
        newpos = copy.deepcopy (pos)
        newpos[1] = pos[1] + radius + (newradius)
        flakeset =  ( generate_sphere_flake ( newpos, newradius, maxdepth, curdepth+1) )
        
        newpos = copy.deepcopy (pos)
        newpos[2] = pos[2] + radius + (newradius)
        flakeset =  ( generate_sphere_flake ( newpos, newradius, maxdepth, curdepth+1) )
        
        newpos = copy.deepcopy (pos)
        newpos[2] = pos[2] - radius - (newradius)
        flakeset =  ( generate_sphere_flake ( newpos, newradius, maxdepth, curdepth+1) )
        
        return flakeset
        
        
        
    
    


def main():
    global CAUSTIC_RADIUS
    scene=[]
    lights=[]

    light=Light ( [-100,-450,0] , (1, 1 ,1) )
    light.num_photons = 555000
    light.power = 0.4
    light.samples = 1
    light.radius = 0.0
    lights.append (light)

    
    posCamera = [0,0,-150]
    lookCamera = [0,0,0]
    upCamera = [0,1,0]
    #posCamera = [0,100,0]
    #lookCamera = [0,0,0]
    #upCamera = [0,0,-1]
    
    
    #box = Box([0,-30,100],30)
    #scene = scene + box.objectlist
    
    
    #sphere=Sphere([0,-40,0],20)
    #sphere.mat.refract=0.8
    #sphere.mat.ior = 1.09
    #sphere.mat.send_caustics = True
    #scene.append (sphere)
    
    #sphere=Sphere([10,-40,40],20)
    #sphere.mat.refract=0.0
    #sphere.mat.reflect = 0.3
    #sphere.mat.colour=[0.2,0.2,1]
    #scene.append (sphere)
    
    sphere=Sphere([0,0,50],30)
    sphere.mat.refract=0.0
    sphere.mat.reflect = 0.3
    sphere.mat.colour=[0.2,1.0,0.2]
    scene.append (sphere)
    
    #mat=Material()
    #mat.colour = [1.0,0.6,0.6]
    #mat.reflect = 0.0
    #mat.refract = 0.0
    #mat.ior = 1.3
    #mat.send_caustics = True
    #scale =332.2
    #scene = scene + (load_obj ("torus.obj", mat, [-2,-30,0], scale, True) )
    
    
    #plane=YPlane()
    #plane.mat.spec_level = 0.5
    #plane.mat.colour=[1,1,1]
    #plane.mat.receive_caustics = True
    #plane.mat.reflect=0.0
    #scene.append (plane)
    
    #torus = Torus()
    #scene.append (torus)
        
    
    # Construct kd-tree for scene
    print "Compiling surface kdtree with",len(scene),"objects...."
    tree = Tree(scene)
    print "Done."
    
    # Make photon tree
    # check needed first
    if CAUSTICS:
        photon_map_needed = False
        for item in scene:
            if item.mat.send_caustics == True and ( item.mat.reflect > 0 or item.mat.refract > 0):
                photon_map_needed = True
                break
        if photon_map_needed:
            photon_list = generate_photon_list (tree, lights) # doesnt make a tree yet - just a list of photons
            photon_tree = Photon_Tree (photon_list)
        else:
            print "No photon tree generated - no objects to interact with."
            photon_tree = None
    else:
        photon_tree = None
    
    
    
    eye=[0.0,0.0,-9.1] #eye coord
    vpw = (270.0/SCREENW) * PIXELSIZE  # pixel width on virtual viewing plane
    vph = (270.0/SCREENH) * PIXELSIZE  # pixel height on virtual viewing plane
    xaspect = float(SCREENW) / float (SCREENH)
    yaspect = float(SCREENH) / float (SCREENW)
    viewport = 10.0
    
    if SCREENW > SCREENH: yaspect = 1.0
    else: xaspect = 1.0
    
    vx_size = viewport * xaspect
    vy_size = viewport * yaspect
    vx_start = -vx_size / 2
    vy_start = -vy_size / 2
    
    screen=pygame.display.set_mode((SCREENW,SCREENH),DOUBLEBUF)
    pygame.display.set_caption ("Mitch Raytracer V2.1")
    pygame.key.set_repeat(500, 20)
    clock=pygame.time.Clock()
    
    while 1:
        # do camera calcs
        look = vec_sub (lookCamera, posCamera)
        #print "vecsub lookCamera, posCamera:",lookCamera,"-",posCamera,"=",look
        vh = vec_cross (look, upCamera)
        vh = vec_norm (vh)
        vv = vec_cross (look,vh)
        vv = vec_norm (vv)
        fl = SCREENW / (2 * tan((0.5 * FOV) * 0.017453292))
        Vp = vec_norm (look)
        Vp[0] = Vp[0] * fl - 0.5 * (SCREENW * vh[0] + SCREENH * vv[0])
        Vp[1] = Vp[1] * fl - 0.5 * (SCREENW * vh[1] + SCREENH * vv[1])
        Vp[2] = Vp[2] * fl - 0.5 * (SCREENW * vh[2] + SCREENH * vv[2])
        print "Time for last frame: ",clock.tick()
        for event in pygame.event.get():
            if event.type==KEYDOWN:
                    if event.key==K_ESCAPE:
                        return
                    if event.key==K_w:
                        posCamera[2]+=1
                        #lookCamera[2]+=10
                    if event.key==K_s:
                        posCamera[2]-=1
                        #lookCamera[2]-=10
                    if event.key==K_a:
                        posCamera[0] -=1
                        #lookCamera[0] -=10
                    if event.key==K_d:
                        posCamera[0] +=1
                        #lookCamera[0] +=10
                    if event.key==K_r:
                        posCamera[1] -=1
                        #lookCamera[1] -=10
                    if event.key==K_f:
                        posCamera[1] +=1
                        #lookCamera[1] +=10
                    if event.key==K_LEFT:
                        lookCamera = yrot (lookCamera ,posCamera, -4)
                    if event.key==K_RIGHT:
                        lookCamera = yrot (lookCamera ,posCamera, 4)
                    if event.key==K_q:
                        upCamera = zrot (upCamera, [0,0,0], -4)
                    if event.key==K_e:
                        upCamera = zrot (upCamera, [0,0,0], 4)
                    if event.key==K_UP:
                        lookCamera = xrot (lookCamera ,posCamera, 4)
                    if event.key==K_DOWN:
                        lookCamera = xrot (lookCamera ,posCamera , -4)
            if event.type==QUIT:
                return
        
        #vpw = PIXELSIZE * vx_size/SCREENW/0.04
        #vph = PIXELSIZE * vy_size/SCREENH/0.04
        for y in xrange(0,SCREENH, PIXELSIZE*SCANLINES):
            for x in range(0,SCREENW, PIXELSIZE):
                direc = [0,0,0]
                direc[0] = x * vh[0] + y * vv[0] + Vp[0]
                direc[1] = x * vh[1] + y * vv[1] + Vp[1]
                direc[2] = x * vh[2] + y * vv[2] + Vp[2]
                direc = vec_norm (direc)
                ray = Ray(posCamera, direc)
			
		#print "x:     ", x
		#print "y:     ", y
		#print "vh:    ", vh[0], vh[1], vh[2]
		#print "vv:    ", vv[0], vv[1], vv[2]
		#print "vp:    ", Vp[0], Vp[1], Vp[2]
		#print "look:  ", look[0], look[1], look[2]
		#print "DIREC: ", direc[0], direc[1], direc[2]
		#raw_input(".")


                f_colour = get_colour ( ray, tree, lights, 0, photon_tree )
                if f_colour[0] < 0:
                    # nothing hit by ray, draw background, or nothing
                    pass
                else:
                    f_colour = colour_clamp (f_colour)
                    if SHOW_TREE:
                        nodes = tree.get_nodes_traversed (ray)
                        nodes = nodes/255.0
                        k_colour = [nodes,nodes,nodes]
                        if SHOW_SCENE:
                            f_colour = colour_add ( k_colour, colour_scalar_mul (f_colour,0.25))
                            f_colour = colour_clamp (f_colour)
                        else:
                            f_colour = colour_clamp (k_colour)
                    rgb = float_to_rgb ( f_colour )
                    rect=pygame.rect.Rect (x, y, PIXELSIZE, PIXELSIZE * SCANLINES)
                    screen.fill ( rgb, rect )
                    
            
            for event in pygame.event.get():
                if event.type==KEYDOWN:
                    if event.key==K_ESCAPE:
                        return()
                    if event.key==K_UP:
                        CAUSTIC_RADIUS = CAUSTIC_RADIUS + 4
                    if event.key==K_DOWN:
                        CAUSTIC_RADIUS = CAUSTIC_RADIUS - 4
                    if event.key==K_SPACE:
                        clock.tick(1)
            pygame.display.update()
                             
        pygame.display.flip()
        #screen.fill(BACKGROUND)


if __name__=='__main__':
    main()
    pygame.quit()

    
