""" basic vector operations """

from math import sqrt

def vec_add (v1,v2):
    return [ v1[0]+v2[0],  v1[1]+v2[1],  v1[2]+v2[2] ]

def vec_sub (v1,v2):
    return [  v1[0]-v2[0],  v1[1]-v2[1],  v1[2]-v2[2] ]

def vec_mul (v1,s):
    """scalar vector multiplication"""
    return [ v1[0]*s, v1[1]*s,  v1[2]*s ]

def vec_div (v,s):
    """scalar vector division"""
    return [ v[0]/s, v[1]/s, v[2]/s ]

def vec_vector_mul (v1,v2):
    """vector multiplication"""
    return [ v1[0]*v2[0]. v1[1]*v2[1], v1[2]*v2[2] ]

def vec_cross (v1,v2):
    """ vector cross product """
    return [ (v1[1] * v2[2]) - ( v1[2] * v2[1]), (v1[2] * v2[0]) - ( v1[0] * v2[2]), (v1[0] * v2[1]) - ( v1[1] * v2[0])  ]

def vec_dot (v1,v2):
    return (v1[0]*v2[0]) + (v1[1]*v2[1]) + (v1[2]*v2[2]) 

def vec_norm (v1):
    modulo = sqrt ( vec_dot(v1,v1) )
    if modulo!=0:
        return  vec_mul (v1,( 1.0 / modulo))
    else:
        return [0,0,0]
    
def vec_distance (v1,v2):
    """ returns the distance between two 3d points """
    d = vec_sub (v1,v2)
    return sqrt ( d[0]**2 + d[1]**2 + d[2]**2)

def vec_square_distance (v1,v2):
    """ returns the distance between two 3d points """
    d = vec_sub (v1,v2)
    return ( d[0]**2 + d[1]**2 + d[2]**2)
    