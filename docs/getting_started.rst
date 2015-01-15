================
 Getting started
================

Open BioOpt file
----------------
BioOpt library comes with a toy model (toy.bioopt) which will be used through the rest of this tutorial.

.. code:: python

    from bioopt.bioopt_parser import BiooptParser
    model = BiooptParser().parse_file("toy.bioopt")
    print(model)


Bioopt model contains a list of reactions reactions and a description of objective function (like biomass)

.. code:: python

    print("Total number of reactions: {0}".format(len(model.reactions)))
    print("Objective: {0}".format(model.objective))

.. parsed-literal::

    Total number of reactions: 11
    Objective: R3 * 1.0 * 1.0

Retrieve reactions
------------------
Because reactions in BioOpt library are regular lists they can be accessed by index starting at 0. Additionally a
a model can be searched for a specific reaction using reaction name or regular expression. Finally, there is a
special function to search all export reactions


.. code:: python

    print(model.find_reaction("R1"))
    print(model.find_reactions(["R1", "R2"]))
    print(model.find_reactions(r".*xtI$", regex=True))
    print("Model have {0} export/import reactions".format(len(model.find_boundary_reactions()))) # search for all export reactions

.. parsed-literal::

    R1[-inf, inf]: 1 A <-> 2 B
    [R1[-inf, inf]: 1 A <-> 2 B, R2[0, inf]: 1 A -> 1 C]
    [A_xtI[0.0, 1.0]: 1 A_xtX* -> 1 A, B_xtI[0.0, 1.0]: 1 B_xtX* -> 1 B, C_xtI[0.0, 1.0]: 1 C_xtX* -> 1 C, D_xtI[0, inf]: 1 D_xtX* -> 1 D]
    Model have 8 export/import reaction


Retrieve metabolites
--------------------
Metabolites can not be accessed as a list but similar to reactions a model can be searched for specific metabolites using
metabolite name or regular expression. And if you won't provide any pattern for metabolite name all metabolites will be
returned. Additionally, a model can be searched for a specific reaction using reaction name or regular expression.
Finally, there is a special function to search all export reactions

.. code:: python

    print("Model contains {0} unique metabolites".format(len(model.find_metabolites())))
    print(model.find_metabolite("A"))
    print(model.find_metabolites(["B", "B_xtX"]))
    print(model.find_metabolites(r"B.*", regex=True))

.. parsed-literal::

    Model contains 8 unique metabolites
    A
    [B, B_xtX*]
    [B, B_xtX*]


Access reaction information
---------------------------
Every reaction include flux constraints, directionality and list of reactants and products.

.. code:: python

    r = model.find_reaction("A_xtI")
    print("Reaction '{react}' has direction {dir} and constraints from {lb} to {ub}".format(react=r.name, dir=r.direction, lb=r.bounds.lb, ub=r.bounds.ub))
    print("Reactants: {0}".format(','.join(rs.metabolite.name for rs in r.reactants)))
    print("Products: {0}".format(','.join(ps.metabolite.name for ps in r.products)))


.. parsed-literal::

    Reaction 'A_xtI' has direction -> and constraints from 0.0 to 1.0
    Reactants: A_xtX
    Products: A

Access metabolite information
-----------------------------

Every metabolite have information about its name and whether it is boundary (exported/imported) or balanced metabolite

.. code:: python

    m = model.find_metabolite("A_xtX")
    print("Metabolite '{0}' is {1}boundary".format(m.name, "not " if not m.boundary else ""))

.. parsed-literal::

    Metabolite 'A_xtX' is boundary

In context of reaction metabolites are wrapped inside `ReactionMember` class which allows access to a coefficient before
metabolite in this reaction

.. code:: python

    r = model.find_reaction("A_xtI")
    print("Metabolite {0} have coefficient {1:.0f}".format(r.reactants[0].metabolite, r.reactants[0].coefficient))

.. parsed-literal::

    Metabolite A_xtX* have coefficient 1