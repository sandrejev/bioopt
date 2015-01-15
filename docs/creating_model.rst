================
 Building model
================

There is two ways to create models. Shorthand notation is useful for quickly prototyping models. It takes less space
and syntax is designed to resemble original BioOpt format. However for people unfamiliar with shorthand notation it
might prove difficult to understand your code. On the other hand explicit notation is unambiguous but longer. It is
also faster if you are trying to improve code performance for big models. The choice is yours but in long-term using
explicit notation is advised.


Explicit notation
==================

First let's try to create a model using explicit notation

.. code:: python

    from bioopt.model import *

    # Create left-hand-side of reaction
    reactants = ReactionMemberList()
    reactants.append(ReactionMember(Metabolite("Na"), coefficient=2))
    reactants.append(ReactionMember(Metabolite("H2O"), coefficient=2))

    # Create right-hand-side of reaction
    products = ReactionMemberList()
    products.append(ReactionMember(Metabolite("NaOH"), coefficient=2))
    products.append(ReactionMember(Metabolite("H2"), coefficient=1))

    # Combine previously created objects into complete reaction definition.Note that direction and constraints are not
    # required. By default reaction would be unconstrained and allow to flow both directions. However if direction is
    # specified default constraints will reflect reaction directionality.
    r1 = Reaction("Burn", reactants, products, Direction.forward(), Bounds(-1, Bounds.inf())) # Lower bound will be reset to 0 because reaction is not reversible
    r2 = Reaction("Break", products, reactants)

    # Let's now create model with this reaction
    m = Model()
    m.reactions = [r1, r2]
    print(m)

.. parsed-literal::

   -REACTIONS
   Burn[0, inf]     : 2 Na + 2 H2O   -> 2 NaOH + 1 H2
   Break[-inf, inf] : 2 NaOH + 1 H2 <-> 2 Na + 2 H2O

   -CONSTRAINTS
   Burn    [0, inf]
   Break   [-inf, inf]


Short-hand notation
====================

Now let's create same model using short-hand notation. Notice how more convenient is notation for reactants and
products

.. code:: python

    reactants = 2*Metabolite("Na") + 2*Metabolite("H2O")
    products = 2*Metabolite("NaOH") + 2*Metabolite("H2")
    r1 = Reaction("Burn", reactants, products, Direction.forward(), Bounds(0, Bounds.inf()))
    r2 = Reaction("Break", products, reactants)
    m = Model()
    m.reactions = [r1, r2]


Write BioOpt file
=================

Once you have a created BioOpt model it can be saved on disc by calling `save()` member function.

.. code:: python

    m.save("export.bioopt")