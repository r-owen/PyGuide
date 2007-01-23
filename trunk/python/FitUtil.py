"""Some basic fitting routines from scipy 0.5.

Presently includes:
brent       ---     Use Brent's method to minimize a 1-d function
                    (does not need inital guess)
bracket     ---      Find a bracket containing the minimum.
"""
# ******NOTICE***************
# optimize.py module by Travis E. Oliphant
#
# You may copy and use this module as you see fit with no
# guarantee implied provided you keep this notice in all copies.
# *****END NOTICE************

# A collection of optimization algorithms.  Version 0.5
# CHANGES
#  Bug fixes: the two raised exceptions were misspelled
#  (caught by pychecker) (R. Owen. 2006-04-17)
#  Added fminbound (July 2001)
#  Added brute (Aug. 2002)
#  Finished line search satisfying strong Wolfe conditions (Mar. 2004)
#  Updated strong Wolfe conditions line search to use cubic-interpolation (Mar. 2004)

# Minimization routines
__all__ = ["bracket", "brent"]

def bracket(func, xa=0.0, xb=1.0, args=(), grow_limit=110.0):
    """Given a function and distinct initial points, search in the downhill
    direction (as defined by the initital points) and return new points
    xa, xb, xc that bracket the minimum of the function:
    f(xa) > f(xb) < f(xc)
    """
    _gold = 1.618034
    _verysmall_num = 1e-21
    fa = func(xa, *args)
    fb = func(xb, *args)
    if (fa < fb):                      # Switch so fa > fb 
        dum = xa; xa = xb; xb = dum
        dum = fa; fa = fb; fb = dum
    xc = xb + _gold*(xb-xa)
    fc = func(xc, *args)
    callNum = 3
    iterNum = 0
    while (fc < fb):
        tmp1 = (xb - xa)*(fb-fc)
        tmp2 = (xb - xc)*(fb-fa)
        val = tmp2-tmp1
        if abs(val) < _verysmall_num:
            denom = 2.0*_verysmall_num
        else:
            denom = 2.0*val
        w = xb - ((xb-xc)*tmp2-(xb-xa)*tmp1)/denom
        wlim = xb + grow_limit*(xc-xb)
        if iterNum > 1000:
            raise RuntimeError("Too many iterations")
        if (w-xc)*(xb-w) > 0.0:
            fw = func(w, *args)
            callNum += 1
            if (fw < fc):
                xa = xb; xb=w; fa=fb; fb=fw
                return xa, xb, xc, fa, fb, fc, callNum
            elif (fw > fb):
                xc = w; fc=fw
                return xa, xb, xc, fa, fb, fc, callNum
            w = xc + _gold*(xc-xb)
            fw = func(w, *args)
            callNum += 1
        elif (w-wlim)*(wlim-xc) >= 0.0:
            w = wlim
            fw = func(w, *args)
            callNum += 1
        elif (w-wlim)*(xc-w) > 0.0:
            fw = func(w, *args)
            callNum += 1
            if (fw < fc):
                xb=xc; xc=w; w=xc+_gold*(xc-xb)
                fb=fc; fc=fw; fw=func(w, *args)
                callNum += 1
        else:
            w = xc + _gold*(xc-xb)
            fw = func(w, *args)
            callNum += 1
        xa=xb; xb=xc; xc=w
        fa=fb; fb=fc; fc=fw
    return xa, xb, xc, fa, fb, fc, callNum

def brent(func, args=(), brack=None, tol=1.48e-8, full_output=False, maxiter=500):
    """Given a function of one-variable and a possible bracketing interval,
    return the minimum of the function isolated to a fractional precision of
    tol.

    Uses inverse parabolic interpolation when possible to speed up convergence
    of golden section method.
    
    Inputs:
    - func  the function to minimize;
            - the first argument is x, the value being found
            - additional arguments may be provided in the args argument
    - args  a sequence of one or more additional arguments for func
            (supplied after x, the first argument)
    - brack a triple (a,b,c) where (a<b<c) and func(b) < func(a),func(c).
            If bracket is two numbers then they are assumed to be
            a starting interval for a downhill bracket search
            (see bracket)
    - tol   final precision (if < _mintol, silently increased)
    - full_output   if True, return xmin, func(xmin), # iters, # func calls
            else return xmin
    - maxiter   maximum number of iterations
    """
    _mintol = 1.0e-11
    _cg = 0.3819660
    if brack is None:
        xa,xb,xc,fa,fb,fc,callNum = bracket(func, args=args)
    elif len(brack) == 2:
        xa,xb,xc,fa,fb,fc,callNum = bracket(func, xa=brack[0], xb=brack[1], args=args)
    elif len(brack) == 3:
        xa,xb,xc = brack
        if (xa > xc):  # swap so xa < xc can be assumed
            dum = xa; xa=xc; xc=dum
        assert ((xa < xb) and (xb < xc)), "Not a bracketing interval."
        fa = func(xa, *args)
        fb = func(xb, *args)
        fc = func(xc, *args)
        assert ((fb<fa) and (fb < fc)), "Not a bracketing interval."
        callNum = 3
    else:
        raise ValueError("Bracketing interval must be length 2 or 3 sequence")

    x=w=v=xb
    fw=fv=fx=func(x, *args)
    if (xa < xc):
        a = xa; b = xc
    else:
        a = xc; b = xa
    deltax= 0.0
    callNum = 1
    iterNum = 0
    while (iterNum < maxiter):
        tol1 = tol*abs(x) + _mintol
        tol2 = 2.0*tol1
        xmid = 0.5*(a+b)
        if abs(x-xmid) < (tol2-0.5*(b-a)):  # check for convergence
            xmin=x; fval=fx
            break
        if (abs(deltax) <= tol1):           
            if (x>=xmid): deltax=a-x       # do a golden section step
            else: deltax=b-x
            rat = _cg*deltax
        else:                              # do a parabolic step
            tmp1 = (x-w)*(fx-fv)
            tmp2 = (x-v)*(fx-fw)
            p = (x-v)*tmp2 - (x-w)*tmp1;
            tmp2 = 2.0*(tmp2-tmp1)
            if (tmp2 > 0.0): p = -p
            tmp2 = abs(tmp2)
            dx_temp = deltax
            deltax= rat
            # check parabolic fit
            if ((p > tmp2*(a-x)) and (p < tmp2*(b-x)) and (abs(p) < abs(0.5*tmp2*dx_temp))):
                rat = p*1.0/tmp2        # if parabolic step is useful.
                u = x + rat
                if ((u-a) < tol2 or (b-u) < tol2):
                    if xmid-x >= 0: rat = tol1
                    else: rat = -tol1
            else:
                if (x>=xmid): deltax=a-x # if it's not do a golden section step
                else: deltax=b-x
                rat = _cg*deltax

        if (abs(rat) < tol1):            # update by at least tol1
            if rat >= 0: u = x + tol1
            else: u = x - tol1
        else:
            u = x + rat
        fu = func(u, *args)      # calculate new output value
        callNum += 1

        if (fu > fx):                 # if it's bigger than current
            if (u<x): a=u
            else: b=u
            if (fu<=fw) or (w==x):
                v=w; w=u; fv=fw; fw=fu
            elif (fu<=fv) or (v==x) or (v==w):
                v=u; fv=fu
        else: 
            if (u >= x): a = x
            else: b = x
            v=w; w=x; x=u
            fv=fw; fw=fx; fx=fu
        
    xmin = x
    fval = fx
    if full_output:
        return xmin, fval, iterNum, callNum
    else:
        return xmin
