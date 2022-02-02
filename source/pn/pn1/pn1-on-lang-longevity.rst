:Author: Ammar Hakim
:Date: Jan 14th 2022
:Completed: 
:Last Updated: 

JE1: Why Some Programming Languages Survive
===========================================

Why do certain programming languages have such longevity? Is it merely
the perverseness of some prominent people (bullying or coercion)? Or is
there a fundamental reason for this? Why has not 50 years of research by
“the smartest and brightest” people led to significant evolution of the
most widely used (but old) languages? I mean, after all, don’t the
brightest cream of humanity go into computer “science”? So what is the
explanation that C remains so widely used? (C and C++ combined are 2x as
popular as the next option, Python. Further, many other highly popular
languages like C# and Java are also very C-like in their basic
semantics).

The reasons for these are manifold. First, just because 50 years of
research have gone by it does not mean that any fundamental change has
occurred in the basis of computing. In fact, the situation should be
compared to that of calculus: centuries have gone by and we still use
the same notation as invented by Leibniz and concepts invented by him
and Newton. Similarly, 2300 years later, the essential ideas of
hydrostatic equilibrium are unchanged from what Archimedes laid down in
his treatises. Maxwell’s equations, once written down using vector
calculus 150 years ago, have remained unchanged. The reason for this (in
mathematics and physics) is obvious: it reflects, in a deep way, the way
our world works.

A similar situation exists in computing. The fundamental structure of
our hardware dictates what languages will survive in the long run.
Thousands of languages have come and gone, many invented by very famous
and bright people. However, the hardware architecture (von Neumann
architecture) has not changed and this means that languages that target
or mimic this architecture in software will naturally survive. The most
important language that was specifically invented to program machines
with the von Neumann architecture is C. The explicit goal of C design
was to be one level above assembly that the hardware itself understands.
This means that as long as the von Neumann architecture is around, C
will continue to dominate and survive. Even languages that wish to
supplant C will need to follow its lead in being compact, have few
orthogonal concepts that can be combined in non-trivial ways and have
the ability to talk to the hardware directly. This essentially explains
C’s longevity and continues wide-spread use and popularity. (Though it
would not *appear* C is that popular given all the FUD and lack of
online “hits” compared to JavaScript, Java and C++, for example).

Further, there is a significant cost to software abstractions
[1]_. The notion of “zero cost abstractions” is perversely and
criminally false.  This falsehood about abstractions having little
cost or that “computationally intensive” parts can be rewritten in the
future in a more optimal way, has led to a situation in which the
software bloat has out-paced hardware speedup. It is one thing to
demonstrate a small example of abstractions having little cost, but a
totally different thing to maintain that this small example actually
translates to production-level code. The former may be true, but the
latter is most certainly (and criminally) false.

Hence, modern software appears sluggish even on the fastest hardware.
Efficiency matters, and despite reassurances from the “brightest and
best of humanity”, efficiency can’t be easily (or at all) fixed later.
It is one thing for Knuth to say “premature optimization is the root of
all evil” as his genius is to write highly optimal solutions from the
get go. However, when this becomes a mantra in the hands of other bright
and best people (though say 1000x less bright than Knuth) it leads to a
disastrous mess of inefficient and resource hungry software that
consumes gargantuan amounts of energy and yet underperforms.

There are serious cost of abstractions (please see footnote) in
computing.

-  Abstractions in computing is hard. Tradeoffs are inevitable and not
   simply restricted to performance. For example, unlike mathematics in
   which integers can be unboundedly large and real numbers have
   infinite expansions, in a computer they must remain finite. Hence,
   knowledge of the hardware limitations must be accounted for when
   doing numeric-intensive work.
-  Memory is finite and linearly laid-out in a von Neumann architecture.
   Programs are treated as data (“stored program architecture”, first
   invented by Charles Babbage). Hence, the language must allow proper
   means of directly addressing the underlying memory and instructions
   efficiently via close-to-metal calls. This is particularly true when
   complex memory hierarchies, SIMD operations, pipelined CPUs and
   massively parallel threaded devices like GPUs are becoming
   commonplace. Abstractions here can be disastrously inefficient.
-  *Thin* abstract interfaces are good. C interfaces are ultimate in
   “thin-ness”: structs and functions combined together with pointers
   (direct memory manipulation). C++ and Java interfaces are often
   expressed as deep inheritance hierarchies. These are difficult to
   understand and highly brittle: modifying them is very difficult and
   probably impossible in a production system.
-  However, being one layer above assembly involves a cost in terms of
   safety, in that care must be taken to ensure one does not cross the
   process boundary. This requires discipline and consistent use of both
   static and dynamic analysis tools.
-  Safety in the fundamental sense of process-boundary violations is a
   false notion: errors can occur even in the most “safe” languages.
   Further, notional safety comes with large runtime as the concept of
   process-boundary must be abstracted (VMs). These runtimes hide more
   than you eventually want.
-  Language runtimes and especially garbage-collection (needed for
   safety and blithely using resources) adds unpredictable overhead.
   Further, reference leaking can still occur.
-  Error handling is also driven by the underlying von Neumann
   architecture. Obviously, hardware can’t throw exceptions! Hence, the
   hardware signals must be caught and propagated backwards via error
   codes, undoing changes as needed. Due to non-locality of exceptions,
   the only real option is often to merely abort.
-  Resource management in critical code must be done with care:
   delete/release when done, and don’t use after delete.
-  Abstraction layers in some languages can he hideous: innumerable
   constructors, complex operator overloading rules, polymorphism,
   implicit conversions, inheritance and virtual methods, generics,
   template meta-programming, exceptions, pure-interfaces, innumerable
   “decorators” etc, etc etc.

Another language that has survived, though not with the same popularity
as C, is Lisp. Here again, the reason for survival is Lisp reflects the
fundamental mathematical nature of computing. This theory was developed
in the early parts of the 20th century by Turning, Church and others
based on work emerging from the quest to solidify the foundations of
mathematics. Lisp is a concrete implementation of these fundamental
mathematical ideas. One reason, though, why Lisp did not become dominant
that the hardware to run Lisp efficiently never took off, and was soon
totally eclipsed by machines based on the von Neumann architecture.
Could things have turned out different, and Lisp become dominant
compared to C? I am not sure, as the fundamental data-structure in Lisp
(a cons-cell) requires indirection. So it is possible that the von
Neumann architecture is actually mathematically inevitable when
efficiency is accounted for. (Though Babbage invented the stored-program
concept in the 18th century, it appears that von Neumann did not know
about his work and independently rediscovered and fully developed it
right after WW-II).

Fortran and APL family of languages also have had significant longevity.
Again, this is as multi-dimensional array manipulations are fundamental
when mathematics is performed on a computer.

It hence appears that the reason the three language families (C, Lisp
and Fortran) have survived essentially unchanged over 50 years now is
that they are based on fundamental and, perhaps inevitable, underlying
mathematical structures. Hence, one can read a C or Lisp book from 1990s
and still find them refreshingly modern. Meanwhile, even important
landmark tomes like Knuth’s Art of Computer Programming are not studied
for the code they contain, but more for the mathematical (algorithmic)
analysis they develop and exploit. Thus, as Brother William would say,
as mathematics underlies Nature, it is natural and good that our
machines and language should also express this fundamental feature of
our Universe.

.. [1] I do not mean to imply *all* abstractions are bad!
   Mathematically well-founded abstractions like arrays and other
   data-structures, functions that act on them and the like are just
   fine. The issue comes about when the abstractions are too
   "high-level". The word "abstractions" perhaps is too ambiguous and
   means different thing to different people. At present I do not know
   a more precise term to use.
