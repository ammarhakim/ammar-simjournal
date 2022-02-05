:Author: Ammar Hakim
:Date: March 25th 2021
:Completed: 
:Last Updated: Jan 30th 2022

PN0: A Minimalist Approach to Software
======================================

.. epigraph::

   The fate of a writer is strange. He begins his career by being a
   baroque writer, pompously baroque, and after many years, he might
   attain if the stars are favorable, not simplicity, which is
   nothing, but rather a modest and secret complexity.

   -- Jorge Luis Borges, Prologue to the "Self and the Other"

This foundational document defines my approach to programming. This is
a minimalist approach: programs are built with a tight and elegant
code structure, with an focus on performance. This forces a choice of
the *core* programming language (the C programming language) and
programming model (procedural with a mathematically clean separation
of data and operations). Despite the 50 year history of C (or perhaps
because of it), it is the only serious programming language ever
invented and will continue to be so as long as the von Neumann
architecture, to which it is closely tied, remains dominant.

The document is divided into sections, with a few short aphorisms
making up each section. The style is deliberately terse and
prescriptive.

.. contents::

Minimalism defined
------------------

#. What is minimalism? Minimalism does not mean writing less code, but
   writing code that counts. Minimalist programs are elegant and have a
   tight code structure, and do one thing well.

#. *Minimalism* means removing all superfluous features and code not
   needed to achieve the minimum viable program (MVP). One may term
   this approach as *brutal minimalism*. Some features are "good to
   have" but if they are not "must have" they should be eliminated
   when aiming for a minimalistic design.

#. *Minimalism* goes hand-in-hand with simplicity. A program should be
   as simple as possible. Over-engineering is a serious issue in
   software development, leading to code bloat, unnecessary layering
   that eventually leads to an unmaintainable mess.

#. Minimalist programs do not require complex processes, described in
   long, never-to-be-read documents. Processes documents become stale
   quickly. Processes also encourage a bureaucratic system, which is
   completely antithetical to minimalist programming. Hence, minimize
   process.

#. Programming is never about lines of code or less typing or other
   such superficial measures. (Though all these are typically the
   outcome of minimalist design). Programming is about expressing
   executable ideas cleanly. It is a very difficult art and requires
   removing the fear of (full or partial) rewrites and a sharp,
   mathematical and axiomatic focus on minimal concepts required to
   implement features efficiently.

Separate data and functions
---------------------------

#. In our notation an object is a unit in which *related data* are
   kept together. Example: instances of a C struct containing
   plain-old-data (POD).

#. To achieve a clean design, *data* and *operators* on those data
   must be kept separate. This is the *mathematically correct* thing
   to do as it allows constructing different systems of functions to
   manipulate the same data in an independent and non-intrusive
   manner.

#. Data, once created, should be treated as read-only and not directly
   modified. Data modification should only happen via functions.

#. Separation of data and operators allows dispatch on multiple object
   types. That is, functions can be written that take two or more
   objects to perform an action. This removes the incestuous state
   sharing that occurs when data and operators are mixed.   

Extensible code
---------------

#. Very often extensibility is not important. In such situations only
   the special case should be handled, but handled well.

#. Extensibility should not be based on the existence of common terms
   in ordinary language to describe two otherwise disparate
   systems. Ordinary language is not precise enough to express
   commonality and only through very careful analysis one discovers
   commonality (or lack thereof).

#. Different systems should not be shoehorned into one without
   significant analysis. In fact, extensibility usually is *increased*
   when systems are cleanly separated, but allow structured
   communication between them. Consider Unix command-line tools and
   their simple and elegant chaining mechanisms via pipes and
   output/input redirection.

#. When extensibility is required implement it with a minimalist
   design without complicated class hierarchies and fat interfaces. A
   hierarchical class structure that first suggests itself usually
   does not work cleanly in practice. Data nesting is fine, class
   inheritance is not as it leads to incestuously shared state.

#. Use code layering as indirections and not bandages. Do not add yet
   more layers to hide bad code. Rewrite it.   

Minimize dependencies
---------------------
   
#. Dependencies should be minimized. Avoid dependencies you do not
   understand completely.

#. Avoid a dependency if it needs two or more dependencies of its own.

#. Avoid dependency management tools that install dozens of packages
   from scratch.

#. Do not use popularity as a metric to understand an existing
   software library/framework. Some popular libraries may have some
   high-quality code but more often popularity is simply an indicator
   of good marketing (funding pressures or corporate greed to
   establish platform tie-in).

#. Minimalist programs do not require including
   everything-under-the-sun frameworks. In fact, it is a good idea to
   avoid anything that has the word "framework" or other such
   buzzwords in them.

Use C and if needed, add a scripting layer
------------------------------------------

#. Stay close to the machine. Learn the details of the von Neumann
   architecture and keep the structure of memory and the execution
   engine (CPU/GPU) in mind while designing code.

#. At first one wants results but very soon one wants control. To
   achieve results *and* control, use the C programming language. It
   gives you complete control, though at the price of great
   discipline. Learn to be disciplined.

#. C is the *only* serious programming language ever invented. A
   simple proof of this statement is to look at the core compute
   kernel of *any* critical software: it will be written in C.

#. Proper use of C structs and function pointers can lead to
   surprisingly elegant designs and clean, tight code.

#. Managing memory yourself is not a burden. Highly robust and
   reliable software like the Linux kernel, `redis
   <https://redis.io/>`_, `haproxy <https://www.haproxy.org/>`_,
   `sqlite <https://sqlite.org/index.html>`_ etc are written in C and
   all of them manage memory manually.

#. Pay attention to all compiler warnings and use a static analysis
   tool. Learn to use `valgrind <https://valgrind.org/>`_ and all the
   tools it provides. Ensure all code is "valgrind clean".

#. C code and C APIs are very easy to bind in multiple
   languages. Hence a good architectural motif (used in redis, haprox
   and elsewhere) is to write the low-level performance critical code
   in C and use scripting to provide higher level control.

#. Modern scripting languages are very flexible and powerful. Some
   like `Lua <https://www.lua.org/>`_ are specially designed for
   embedding in larger applications. Lua has a very tiny footprint, is
   written in portable C, making it universally usable on all types
   of systems, however constrained.

#. Defer complex control to the scripting layer. Higher-level
   scripting languages allow more complex and elegant control
   structures even when they are missing from the low-level language
   used to implement the performance critical aspects of the code.

#. The API exposed to the scripting language should be fine-grained
   enough to allow them to be used from of complex control structures
   like lexical closure, coroutines and iterators.

#. Allow users the ability to pass structured data between the script
   and compiled layer.

Focus on performance; measure it carefully
------------------------------------------

#. Focus on writing performant code. Inefficiencies can’t be easily (or
   at all) fixed later.

#. It is one thing for Knuth to say "premature optimization is the
   root of all evil" as his genius is to write highly optimal
   solutions from the get go. You are not Knuth.

#. Use the Linux perf tool to measure performance. Record the number
   of instructions run, the instructions-per-cycle and the chip speed
   for the run.

#. Aim to minimize instructions run (better algorithms) while
   maximizing instructions-per-cycle (unrolling code and avoid
   anything that interrupts the CPU like thread context
   switches). Minimize cache load misses, especially for L1-dcache (do
   as much work as possible with loaded memory and avoid indirections
   in inner loops).

Use a simple build system and automate testing
----------------------------------------------

#. Minimalist programs should be quick to build. Incremental builds
   should not take more than a few seconds and clean rebuild should
   not take more than a few minutes

#. Untested code might as well not exist. Maximize code coverage using
   unit (individual functions and structures) and regression (whole
   system) testing. Ensure there is a "make check" target that runs
   all unit tests.

In summary: creating efficient and innovative software requires a
minimalist approach. The goal should be to construct one or more
minimalist programs that have structured data exchange protocols
instead of giant monolithic programs. Frequent rewrites and
refactoring may be needed before one discovers the correct
design. Monolithic programs and over-engineered systems are almost
invariably slower, harder to maintain (despite their developers having
used the latest object-oriented-programming and "Agile" fads to make
them extensible) and difficult to understand.
