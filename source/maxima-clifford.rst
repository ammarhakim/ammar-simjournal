:Author: Ammar Hakim

On Using the Clifford Algebra Package in Maxima and Expression Manipulation
===========================================================================

.. content::

The Maxima Clifford Package
---------------------------

Geometric (Clifford) Algebra (GA) is a very powerful formalism that
extends standard vector calculus to arbitrary dimensions, including
manifolds of arbitrary signature. GA unifies vast areas of geometry,
from complex numbers, quaternions, exterior algebra, differential
forms, tensor algebra, projective geometry, and many other apparently
disparate mathematical systems. It should be apparent that it is not
trivial to quickly learn such a vast formalism (though the benefits of
doing so are huge), and certainly one can't learn it from this short
technical note. Here, my goal is to show how one can use the `Maxima
Clifford package <https://github.com/dprodanov/clifford>`_ to
manipulate systems of equations for use in various code generators in
the Gkeyll code chain.

To install this package simply clone the repo and add the path to it
in your Maxima load path. Please see the `Gkeyll documentation page
<https://gkeyll.readthedocs.io/en/latest/dev/onmaxima.html>`_ on how
to do this.


