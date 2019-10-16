# Newton's method for polynomials

import copy

def poly_diff(poly):  
    """ Differentiate a polynomial. """  
  
    newlist = copy.deepcopy(poly)  
  
    for term in newlist:  
        term[0] *= term[1]  
        term[1] -= 1  
  
    return newlist  
  
def poly_apply(poly, x):  
    """ Apply a value to a polynomial. """  
  
    sum = 0.0 # force float  
  
    for term in poly:  
        sum += term[0] * (x ** term[1])  
  
    return sum  
  
def poly_root(poly, start, n, r_i):  
    """ Returns a root of the polynomial, with a starting value. 
    To get both roots in a quadratic, try using with n = 1 and n = -1."""  
  
    poly_d = poly_diff(poly)  
    x = start  
    counter = 0  
  
    while True:  
        if (n >= 0) and (n < 1):  
            break  
  
        x_n = x - (float(poly_apply(poly, x)) / poly_apply(poly_d, x))  
  
        if x_n == x:  
            break  
  
        x = x_n  
  
        n -= 1  
        counter += 1  
  
    if r_i:  
        return [x, counter]  
    else:  
        return x  