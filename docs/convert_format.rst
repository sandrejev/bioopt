=========================
 Converting BioOpt model
=========================

Convert to libSBML object
==========================

Toy model from BioOpt package will be used to demonstrate how to convert BioOpt model to libSBML model object.
Consult `libSBML documentation <http://sbml.org/Software/libSBML/docs/python-api/index.html>`_ for additional
information on handling libSBML objects.

.. code:: python

     from bioopt.bioopt_parser import BiooptParser
     from bioopt.converter import Bioopt2SbmlConverter

     model = BiooptParser().parse_file("toy.bioopt")
     converter = Bioopt2SbmlConverter(compartment_pattern=None)
     sbml = converter.convert(model)

     print(sbml.model.toSBML())





