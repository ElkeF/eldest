##########################################################################
#                       COMPLEX INTEGRATION                              #
##########################################################################
# Purpose:                                                               #
#          - A python module to calculate integrals of complex numbers.  #
#                                                                        #
##########################################################################
# written by: Elke Fasshauer May 2018                                    #
##########################################################################

import numpy
import scipy.integrate as integrate

#-------------------------------------------------------------------------
#   integration


def complex_quadrature(func, a, b, **kwargs):
    def real_func(x, **kwargs):
        return numpy.real(func(x, **kwargs))
    def imag_func(x, **kwargs):
        return numpy.imag(func(x, **kwargs))
    real_integral = integrate.quad(real_func, a, b, **kwargs)
    imag_integral = integrate.quad(imag_func, a, b, **kwargs)
    return (real_integral[0] + 1j*imag_integral[0], real_integral[1:],
            imag_integral[1:])
# This function was taken from:
# https://stackoverflow.com/questions/5965583/use-scipy-integrate-quad-to-integrate-complex-numbers

def complex_romberg(func, a, b, **kwargs):
    def real_func(x):
        return numpy.real(func(x))
    def imag_func(x):
        return numpy.imag(func(x))
    real_integral = integrate.romberg(real_func, a, b, **kwargs)
    imag_integral = integrate.romberg(imag_func, a, b, **kwargs)
    return real_integral + 1j*imag_integral

def complex_double_quadrature(outer, inner, a, b, gfun, hfun, **kwargs):
    first_real = lambda y,x: numpy.real(outer(x)) * numpy.real(inner(y))
    sec_real   = lambda y,x: - numpy.imag(outer(x)) * numpy.imag(inner(y))
    first_imag = lambda y,x: numpy.imag(outer(x)) * numpy.real(inner(y))
    sec_imag   = lambda y,x: numpy.real(outer(x)) * numpy.imag(inner(y))

    def temp_ranges(*args):
        return [gfun(args[0]) if callable(gfun) else gfun,
                hfun(args[0]) if callable(hfun) else hfun]

    first_real_integral = integrate.nquad(first_real, [temp_ranges, [a,b]])
    sec_real_integral   = integrate.nquad(sec_real, [temp_ranges, [a,b]])
    first_imag_integral = integrate.nquad(first_imag, [temp_ranges, [a,b]])
    sec_imag_integral   = integrate.nquad(sec_imag, [temp_ranges, [a,b]])

#    print first_real_integral
#    print sec_real_integral
#    print first_imag_integral
#    print sec_real_integral

    return (first_real_integral[0] + sec_real_integral[0]
            + 1j*first_imag_integral[0] + 1j*sec_imag_integral[0],
            first_real_integral[1:])
