""" functions to load 3d objects """

from __future__ import with_statement
from surfaces import Tri
from surfaces import Material
from vector import vec_add, vec_norm
import copy

def load_obj (filename, mat=Material(), p=[0,0,0], s=1.0, SMOOTH=False):
    """ load list of vertices and faces from obj file """
    """ return a list of objects handled by ray tracer """
    """ p = position, s = scale """
    print "Loading object:",filename
    vertex_list=[]
    vertex_list.append ([-1,-1,-1])   # null vertex at start, because obj file starts face counting at 1 not 0
    vertex_normals=[]
    vertex_normals.append ([-1,-1,-1]) # see comment above
    tri_list=[]
    with open(filename) as f:
        for line in f:
            if len(line) > 0 and line.count(" ") > 0:
                if line.split()[0]=='v':
                    v1 = p[0] + (float(line.split()[1]) * s)
                    v2 = p[1] + (float(line.split()[2]) * -s)
                    v3 = p[2] + (float(line.split()[3]) * -s)
                    vertex_list.append ( [v1,v2,v3])
                if line.split()[0]=='f':
                    v1 = int(line.split()[1].split("/")[0])
                    v2 = int(line.split()[2].split("/")[0])
                    v3 = int(line.split()[3].split("/")[0])
                    #print "Make a face from verts: ",v1,v2,v3
                    vertex1 = vertex_list[v3]
                    vertex2 = vertex_list[v2]
                    vertex3 = vertex_list[v1]
                    tri = Tri ( [vertex1, vertex2, vertex3])
                    if len(vertex_normals) > 1 and SMOOTH:
                        tri.vertex_normal = [ vertex_normals[v3], vertex_normals[v2], vertex_normals[v1] ]
                    tri.mat = mat
                    if SMOOTH == False:
                        tri.flat = True
                    tri_list.append (tri)
                if line.split()[0]=='vn':
                    vn1 = float(line.split()[1])
                    vn2 = float(line.split()[2])
                    vn3 = float(line.split()[3])
                    vertex_normals.append ( [vn1,vn2,vn3] )
    
    if SMOOTH and len(vertex_normals) < 2:
        # calculate own normals
        print "Calculating normals..."
        count = 0
        # make a dictionary of point already calculated
        dic = {}
        for tri in tri_list:
            count = count + 1
            if count%100==0:
                print count,",",len(tri_list) 
            for p in tri.p:
                # for each triangle, and each point in triangle, check against all other points in all other
                # points in all other triangle and configure point normal as average of face normals of touching points
                # first check to see we haven`t calced this point already...
                try:
                    if dic[tuple(p)]:
                        tri.vertex_normal[tri.p.index(p)] = dic[tuple(p)] # point already calced by another face - use that result
                except:
                    # never checked this point before, so do normal calculations and store in dictionary
                    current_normal = copy.deepcopy ( tri.normal )
                    for other_tri in tri_list:
                        if other_tri is not tri:
                            for other_p in other_tri.p:
                                if other_p == p:
                                    # found another point which shares face
                                    current_normal = vec_add ( current_normal, other_tri.normal )
                                    break # dont need to check other points in other triangle :)
                    N = vec_norm (current_normal)
                    tri.vertex_normal[tri.p.index(p)] = N
                    dic[tuple(p)] = N
            tri.flat = False
    
    print "Loading complete."
    return tri_list