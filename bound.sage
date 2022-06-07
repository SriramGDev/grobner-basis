def asymptotic_bound(delta, n, l):
    # returns the asymptotic bound of Bardet, Faug\'ere, and Salvy for a homogeneous system of degree delta and dimension l in a polynomial ring with n variables
    var('x')
    eqn = 1/(1-delta*((x+1)^2-x^2)/((x+1)^3-x^3))-((x+1)/x)^(2*delta)
    lam_0 = eqn.find_root((delta-1)/2,delta-1)
    A = (1-delta^(-1))/(2*pi)*((1+lam_0^(-1))^3-1)/(1+lam_0)^(1+l)
    B = (((lam_0+1)/lam_0)^(2*delta)-1)/(1/lam_0^2-1/(lam_0+1)^2)
    return n*A*B^n




