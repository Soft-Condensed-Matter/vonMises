von Mises Probability density function
=======

* [ `von Mises distribution` ](https://en.wikipedia.org/wiki/Von_Mises_distribution#cite_note-Mardia99-2):
  In probability theory and directional statistics, the von Mises distribution 
  (also known as the circular normal distribution or Tikhonov distribution) 
  is a continuous probability distribution on the circle.


  For an angle $x$ the function is given by:

  $$ f(x; \mu, \kappa) = \frac{\exp \left(\kappa \cos\left( x - \mu \right)\right)}{2 \pi I_{0}(\kappa)}$$,

  where $\mu $ is the location and $\kappa$ the concentration, while $I_{0}(\kappa)$ is the modified Bessel function of the first kind of order zero.


_Implementation_:

 The numerical implementation uses the modified Bessel function of first kind 
 that it is computed according to the numerical recipies in Fortran, code also 
 provided in a separate file (** ModifiedBessel.f90 **). 

 von Mises function for different $\kappa$ concentration values:

 ![Values](https://github.com/Soft-Condensed-Matter/vonMises/blob/master/vonMises.png)
 
