""" basic operations on colour rgb truples """

def colour_mul ( rgb1, rgb2 ):
    return ( rgb1[0]*rgb2[0], rgb1[1]*rgb2[1], rgb1[2]*rgb2[2])

def colour_scalar_mul ( rgb1, s ):
    return ( rgb1[0]*s , rgb1[1]*s, rgb1[2]*s )
        

def float_to_rgb ( rgb ):
    return ( int(rgb[0]*255), int(rgb[1]*255), int(rgb[2]*255) )

def colour_add ( rgb1, rgb2 ):
    return [ rgb1[0]+rgb2[0], rgb1[1]+rgb2[1], rgb1[2]+rgb2[2] ]

def colour_subtract ( rgb1, rgb2 ):
    return ( rgb1[0]-rgb2[0], rgb1[1]-rgb2[1], rgb1[2]-rgb2[2] )

def colour_scalar_div ( rgb1, s ):
    return ( rgb1[0]/s , rgb1[1]/s, rgb1[2]/s )

def colour_clamp (c):
    t=[0,0,0]
    for n in xrange(3):
        if c[n]>1 or c[n]<0:
            if c[n]>1:t[n]=1
            if c[n]<0:t[n]=0
        else:
            t[n]=c[n]
    return t

def colour_average (rgb):
    return (rgb[0]+rgb[1]+rgb[2]/3.0)

def colour_alpha (rgb1, rgb2, alpha):
    a = colour_scalar_mul (rgb1,alpha)
    b = colour_scalar_mul (rgb2,1-alpha)
    return colour_add (a,b)