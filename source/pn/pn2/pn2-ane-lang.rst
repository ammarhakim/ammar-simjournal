:Author: Ammar Hakim
:Date: Oct 23rd 2022
:Completed: 
:Last Updated: Feb 6th 2022

PN2: Concatenative Programming for Data Analysis
================================================

Consider the following function written in a typical programming
language

.. code::

   sq(x) : return x*x

This function takes an input number and multiplies it by itself. The
key observation here is that inside the function body the variable "x"
is bound to the input value of "x". So when calling sq(5), x inside
the function body is bound to the value 5. Although this is a
perfectly good way to design programing languages it requires careful
binding of names to values. This is far from a trivial task, requiring
each function to have a private environment that must be managed in a
way that it does not interact with the calling environment.

An alternative to this approach is to eliminate the need for names
completely, instead using a stack to hold values. In such programming
languages, instead we define the sq function as follows

.. code::

   sq : dup *

To call this function one would write the following program:

.. code::

   5 sq

To run this program, we scan from left to right. When we encounter a
number or a literal (strings, arrays etc) we push it on the
stack. When we encounter a function, we pop the appropriate number of
top values from the stack, run the function with those values and push
the result back in the stack. So, after substituting the definition of
the sq function. the stack looks like

.. code::

   5 dup *
   5 5 *
   25

Note the stack modifier "dup" duplicated the value on top of the stack
and pushed it on the stack. The multiply function, as it is a binary
function, then popped two values, multiplied them, and pushed the
answer back on the stack again. When the function completes, the
answer is on top of the stack.

Note that in this approach we never introduced any variables like "x"
above, hence eliminating the need for binding variables to
values. Instead, we introduced two features: a stack to hold values
and a set of "stack rewriting rules" that modify the stack based on
the function one encounters.

The stack rewrite rule for dup are

.. code::

   dup : (a -- a a)

The left of the "--" is the state of the stack before encountering the
dup function, and the right of "--" the state when the dup function is
complete. The stack rewrite rule of * is

.. code::

   * : ( a -- b)

Note the rule does not say what "b" is. It merely describes the state
of the stack. In this case, the top of the stack is replaced by
another value.     

In this stack-based language, we never need to pass an environment to
a function, nor do we need to worry about binding variables to
values. Everything occurs on the stack, including intermediate
calculations.

Now consider the following function in a typical language

.. code::

   double(x) : return 2*x

To double the square of 5, say, we would need to do double(sq(5)). The
nesting of calls is needed to pass the output of sq to double. In a
concatenative language, on the other hand, the double function is
defined as

.. code::

   double : 2 *

that is, we push two on the stack and then the * function pops two
values, multiplying them and pushing the result back on the stack. So,
to double the sqaure of 5 we would write the following program

.. code::

   5 sq double

When this program completes, the value of 50 will be on top of the
stack. This way of combining functions by merely placing them one
after another is the reason why *concatenative languages* are called
such. Concatenative progamming allows chaining a series of functions,
all working on the present state of the stack, to perform a
computation. All this is achieved without introducing any variables in
the program at all.
