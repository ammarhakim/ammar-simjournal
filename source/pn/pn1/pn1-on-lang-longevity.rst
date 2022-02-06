:Author: Ammar Hakim
:Date: Jan 14th 2022
:Completed: 
:Last Updated: Feb 6th 2022

PN1: Why Some Programming Languages Survive
===========================================

.. contents::

Why do certain programming languages have such longevity? Is there a
fundamental reason for this? Why has 50 years of research not led to
significant evolution of the most widely used (but old) languages? Of
course, many new languages, some with elegant concepts have appeared
over the years. So what is the explanation that C remains so widely
used? (C and C++ combined are 2x as popular as the next option,
Python. Further, many other popular languages like C# and Java are
also C-like in their basic semantics).

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
our world works, at least in some limits.

A similar situation exists in computing. The fundamental structure of
our hardware dictates what languages will survive in the long run.
Thousands of languages have come and gone. However, the hardware
architecture (von Neumann architecture) has not changed and this means
that languages that target or mimic this architecture in software will
naturally survive. The most important language that was specifically
invented to program machines with the von Neumann architecture
is C. The explicit goal of C design was to be one level above
assembly, a language that the hardware itself understands. This means
that as long as the von Neumann architecture is around, C will
continue to dominate and survive. Even languages that wish to supplant
C will need to follow its lead in being compact, have few orthogonal
concepts that can be combined in non-trivial ways and have the ability
to talk to the hardware directly. This essentially explains C’s
longevity and continues wide-spread use and popularity. (Though it
would not *appear* C is that popular given all the fear of C that is
cultivated, and lack of slick online marketing as compared to
JavaScript, Java and C++, for example).

On the cost of abstractions
---------------------------

Certain level of abstractions are essential for programming. However,
abstractions come with a cost. The more distant the abstraction takes
one from the machine, the more expensive it is. Depending on the
problem at hand, some abstractions can be very harmful. For example,
in high-performance computing or in core systems software, any
abstraction that hides the details of the machine is bad. For
computational physics the *correct* level of abstraction is that of
*mathematics*, i.e. arrays and other mathematically defined
data-structures (multi-dimensional index-sets, lists, graphs, trees,
etc) and functions/operators that act on them.

Some recent language evangelists have propagated the phrase "zero cost
abstractions", specially in the context of C++. Unfortunately, in most
production memory- and speed-critical software this concept is
completely false, perhaps dangerously so. This falsehood about
abstractions having little cost or that "computationally intensive"
parts can be rewritten in the future in a more optimal way, has led to
a situation in which the software bloat has out-paced hardware
speedup. It is one thing to demonstrate a small example of
abstractions having little cost, but a totally different thing to
maintain that this small example actually translates to
production-level code. The former may be true, but the latter is most
certainly false.

Hence, modern software appears sluggish even on the fastest hardware.
Efficiency matters, and despite reassurances from the language
evangelists, efficiency can’t be easily (or at all) fixed later. It is
one thing for Knuth to say "premature optimization is the root of all
evil" as his genius is to write highly optimal solutions from the get
go. However, when this becomes a mantra in the hands of others it
leads to a disastrous mess of inefficient and resource hungry software
that consumes gargantuan amounts of energy and yet underperforms.

- Abstractions in computing is hard. Tradeoffs are inevitable and not
  simply restricted to performance. For example, unlike mathematics in
  which integers are unboundedly large and real numbers have infinite
  expansions, in a computer they must remain finite. Hence, knowledge
  of the hardware limitations must be accounted for when doing
  numeric-intensive work.
- Memory is finite and linearly addressed in a von Neumann
  architecture.  Programs are treated as data ("stored program
  architecture", first invented by Charles Babbage). Hence, the
  language must allow efficiently addressing the underlying memory and
  instructions via close-to-metal calls. This is particularly true
  when complex memory hierarchies, SIMD operations, pipelined CPUs and
  massively parallel threaded devices like GPUs are becoming
  commonplace. Abstractions here can be disastrously inefficient.
- *Thin* abstract interfaces are good. C interfaces are ultimate in
  "thin-ness": structs and functions combined together with pointers
  (direct memory manipulation). These thin abstractions compile down
  to efficient assembly as the language allows expressing code in a
  manner that remains close to the machine.
- In languages like C++ and Java, abstractions are often expressed as
  deep inheritance hierarchies. These are difficult to understand and
  highly brittle: modifying them is difficult and probably impossible
  in a production system. Worse, these deep inheritance hierarchies
  can have significant cost, specially in performance-critical
  software.
- However, being one layer above assembly involves a cost in terms of
  safety, in that care must be taken to ensure one does not cross the
  process boundary. This requires discipline and consistent use of
  both static and dynamic analysis tools.
- Safety in the fundamental sense of process-boundary violations is a
  false notion: errors can occur even in the most "safe" languages.
  Further, notional safety comes with large runtime as the concept of
  process-boundary must be abstracted into virtual machines. These
  runtimes hide more than you eventually want.
- Language runtimes and especially garbage-collection (needed for
  safety and automatic management of resources) adds unpredictable
  overhead. Further, reference leaking can still occur.
- Error handling is also driven by the underlying von Neumann
  architecture. Hence, errors must be caught and propagated backwards
  via error codes, undoing changes as needed.
- C++ exceptions create *non-locality of the control flow*. In
  general, it is not obvious where and why an exception was thrown and
  the only real option is often to merely abort.
- Further, using exceptions to implement control flow for errors is
  inherently clumsy. It is better to design state machines that
  properly handle recoverable errors and gracefully exit when the
  error is not recoverable.
- Resource management in critical code must be done with care:
  delete/release when done, and don’t use after delete. In C one can
  do this cleanly by using gotos and compiler-extensions that allow
  you to call cleanup code when an object goes out of scope.
- Abstraction layers in some languages can he hideous: innumerable
  constructors, complex operator overloading rules, polymorphism,
  implicit conversions, inheritance and virtual methods, generics,
  template meta-programming, exceptions, pure-interfaces, innumerable
  “decorators” etc, etc etc. Besides putting serious linguistic burden
  on the programmer, these promote programming practices that lead to
  very inefficient code [#]_.

On the elegance of Lisp. Fortran and APL
----------------------------------------
  
Another language that has survived, though not with the same
popularity as C, is Lisp. Here again, the reason for survival is Lisp
reflects the fundamental mathematical nature of computing. This theory
was developed in the early parts of the 20th century by Turning,
Church and others based on work emerging from the quest to solidify
the foundations of mathematics. Lisp is a concrete implementation of
these fundamental mathematical ideas. One reason, though, why Lisp did
not become dominant is that the hardware to run Lisp efficiently never
become mainstream, and was soon totally eclipsed by machines based on
the von Neumann architecture. Could things have turned out different,
and Lisp become dominant compared to C? I am not sure, as the
fundamental data-structure in Lisp (a cons-cell) requires
indirection. So it is possible that the von Neumann architecture is
actually mathematically inevitable when efficiency is accounted
for. (Though Babbage invented the stored-program concept in the 18th
century, it appears that von Neumann did not know about his work and
independently rediscovered and fully developed it right after WW-II).

However, despite run-time inefficiencies, Lisp and other functional
languages are mathematically elegant. Lexical binding is magic and
allows a very powerful way to program. Re-entrant procedures (like
coroutines and continuations) can be used to implement complicated
control structures and iterators over data-structures efficiently.

Fortran and APL family of languages also have had significant
longevity. Again, this is as multi-dimensional array manipulations,
that lie at the heart of Fortran and APL, are fundamental when
mathematics is performed on a computer. APL family of languages (and
its descendants like `J <https://www.jsoftware.com/#/>`_), in
particular, are close both to *mathematics* as well as the *machine*
(as the data-structures and operators have direct representation in C)
and hence continue to have usage in performance-critical applications
decades after they were created.

In conclusion
-------------

It hence appears that the reason the three language families (C, Lisp
and Fortran) have survived essentially unchanged over 50 years now is
that they are based on fundamental and, perhaps inevitable, underlying
mathematical structures. Hence, one can read a C or Lisp book from
1990s and still find them refreshingly modern. Meanwhile, even
important landmark tomes like Knuth’s Art of Computer Programming are
not studied for the code they contain, but more for the mathematical
(algorithmic) analysis Knuth develops and exploits. Thus, as Brother
William would say, as mathematics underlies Nature, it is natural and
good that our machines and language should also express this
fundamental feature of our Universe.

Footnotes
---------

.. [#] C++ is particularly bad at the linguistic overload it has
   introduced on programmers. If one uses C++, it is not enough to
   focus on the problem at hand, but one ends up spending countless
   hours navel-gazing obscure language features and quirks. In fact,
   one ends fighting the language far more than focusing on the
   problem at hand. Sadly, the C++ committee has gone on an offensive,
   releasing new extensions to the language every 3 years, making it
   impossible for an ordinary programmer, with real problems to solve,
   to ever master it fully. C++ is an overweight and unwieldy beast
   that is best avoided, unless one wants an ego boost from being
   called an "C++ expert".
