from model import *
import re


# class Bioopt2CobraPy(object):
#     def parse_file(self, filename = None, name = 'here sould be name of your model'):
#
#         from bioopt.bioopt_parser import BiooptParser
#         import cobra
#         model = BiooptParser().parse_file(filename)
#
#         cobra_model = cobra.Model(name)
#         metabolites = {str(m.name): cobra.Metabolite(m.name) for m in model.find_metabolites()}
#
#         #metabols = dict((str(m.name), cobra.Metabolite(m.name)) for m in model.find_metabolites() if not m.boundary)
#         for r in model.reactions:
#             reaction = cobra.Reaction(r.name)
#             bounds = r.find_effective_bounds()
#             reaction.lower_buond = bounds.lb
#             reaction.upper_bound = bounds.ub
#             reaction.add_metabolites({metabolites[' '.join(str(rr).split(' ')[1::])] : -rr.coefficient*((rr in r.reactants)*2-1) for rr in [i for i in r.reactants + r.products if not i.boundary] })
#             cobra_model.add_reaction(reaction)
#             print reaction.reaction
#         print('%i reaction in model' % len(cobra_model.reactions))
#         print('%i metabolites in model' % len(cobra_model.metabolites))
#         print('%i genes in model' % len(cobra_model.genes))
#         print cobra_model.description
#         return cobra_model

class Bioopt2CobraPyConverter:
    """
    :param inf: A number to substitute __infinity__ values. SBML doesn't have number __infinity__  this number
    (default: 1000) will be used instead
    :rtype: :class:`Bioopt2COBRApyConverter`
    """

    def __init__(self, inf=1000, name = 'model'):
        if isinstance(inf, (float, int)):
            self.inf = inf
        else:
            raise ValueError("Infinity value '{0}' is not a number".format(int))

        if isinstance(name, str):
            self.name = name
        else:
            raise ValueError("Name '{0}' is not a string".format(name))

    def parse_file(self, filename = None):
        from bioopt.bioopt_parser import BiooptParser
        import cobra
        bioopt_model = BiooptParser().parse_file(filename)
        cobra_model = cobra.Model(self.name)
        metabolites = {str(m.name): cobra.Metabolite(m.name) for m in bioopt_model.find_metabolites()}
        #metabols = dict((str(m.name), cobra.Metabolite(m.name)) for m in model.find_metabolites() if not m.boundary)
        for bioopt_reaction in bioopt_model.reactions:
            cobra_reaction = cobra.Reaction(bioopt_reaction.name)
            bounds = str(bioopt_reaction.find_effective_bounds())[1:-1].split(', ')
            cobra_reaction.lower_buond = bounds[0]
            cobra_reaction.upper_bound = bounds[1]
            #cobra_reaction.add_metabolites({metabolites[' '.join(str(rr).split(' ')[1::])] :
             #                                   -rr.coefficient*((rr in bioopt_reaction.reactants)*2-1)
              #                                  for rr in [i for i in bioopt_reaction.reactants + bioopt_reaction.products if not i.boundary] })
            cobra_reaction.add_metabolites(self.__parse_bioopt_reaction(bioopt_reaction))
            cobra_model.add_reaction(cobra_reaction)
        return cobra_model

    def __parse_bioopt_reaction(self, bioopt_reaction):
        import cobra
        #boundary reagents are not taken into account
        meta_dict = {}
        #print bioopt_reaction.boundary()
        for participant in bioopt_reaction.participants:
            if not participant.metabolite.boundary:
                name = participant.metabolite
                if participant in bioopt_reaction.products:
                    coefficient = participant.coefficient
                else:
                    coefficient = -participant.coefficient
                meta_dict[cobra.Metabolite(name)] = coefficient
        return meta_dict

class Bioopt2SbmlConverter:
    """
    Bioopt model object to libSBML document object converter. More details on SBML specification can be found
    `online <http://sbml.org/Documents/Specifications>`_

    :param level: SBML level from 1 to 3 (default: 2)
    :param version: SBML version number (default: 3). Check SBML specification for details
    :param inf: A number to substitute __infinity__ values. SBML doesn't have number __infinity__ so everywhere infinity
            is used in BioOpt it will be replaced with this number (default: 1000)
    :param reaction_id: Reaction id pattern. **"name"** - reaction name will be used as id. **"auto"** - an
            auto-incremental "R_XXXX" value will be used for reaction id (default: **"auto"**)
    :param metabolite_id: Metabolite id pattern. **"name"** - metabolite name will be used as id. **"auto"** - an
            auto-incremental "M_XXXX" value will be used for metabolite id (default: **"auto"**)
    :param compartment_id: Compartment id pattern. **"name"** - compartment name will be used as id. **"auto"** - an
            auto-incremental "C_XXXX" value will be used for compartment id (default: **"auto"**)
    :param compartment_suffix: Specifies whether metabolite names should be appended with compartment name
    :param compartment_pattern: Regular expression pattern describing how to extract compartment information from
            metabolite name. By default everything after last underscore ``_`` character is considered a compartment.
            i.e., metabolite_e - is from compartment "e"
    :rtype: :class:`Bioopt2SbmlConverter`
    """

    def __init__(self, level=2, version=3, inf=1000, reaction_id="auto",
             metabolite_id="auto", compartment_id="auto", compartment_suffix=True, compartment_pattern=r"_(\w+)$"):

        if isinstance(level, int):
            self.level = level
        else:
            raise ValueError("SBML level '{0}' is not an integer".format(level))

        if isinstance(version, int):
            self.version = version
        else:
            raise ValueError("SBML version '{0}' is not an integer".format(version))

        if not compartment_pattern:
            self.compartment_pattern = None
        elif isinstance(compartment_pattern, str):
            self.compartment_pattern = re.compile(compartment_pattern)
        elif hasattr(compartment_pattern, "search"):
            self.compartment_pattern = compartment_pattern
        else:
            raise ValueError("Compartment pattern '{0}' is not regular expression pattern".format(compartment_pattern))

        if isinstance(inf, (float, int)):
            self.inf = inf
        else:
            raise ValueError("Infinity value '{0}' is not a number".format(inf))

        if reaction_id in ["auto", "name"]:
            self.reaction_id = reaction_id
        else:
            raise ValueError("Reaction id pattern can be either 'auto' or 'name' ('{0}' was specified)".format(reaction_id))

        if metabolite_id in ["auto", "name"]:
            self.metabolite_id = metabolite_id
        else:
            raise ValueError("Metabolite id pattern can be either 'auto' or 'name' ('{0}' was specified)".format(metabolite_id))

        if compartment_id in ["auto", "name"]:
            self.compartment_id = compartment_id
        else:
            raise ValueError("Compartment id pattern can be either 'auto' or 'name' ('{0}' was specified)".format(compartment_id))


        if isinstance(compartment_suffix, bool):
            self.compartment_suffix = compartment_suffix
        else:
            raise ValueError("Expected boolean value for compartment_suffix ('{0}' was specified)".format(compartment_suffix))

    def __to_sbml_id(self, name):
        idStream = []
        count = 0
        end = len(name)

        if '0' <= name[count] <= '9':
            idStream.append('_')

        last_character = ''
        for count in range(0, end):
            if '0' <= name[count] <= '9' or 'a' <= name[count] <= 'z' or 'A' <= name[count] <= 'Z':
                last_character = name[count]
                idStream.append(last_character)
            elif last_character != '_':
                idStream.append('_')
                last_character = '_'

        id = ''.join(idStream)
        if id[len(id) - 1] != '_':
            return id

        return id[:-1]

    def __get_valid_sbml_id(self, name, existingIds=()):
        baseString = self.__to_sbml_id(name)
        id = baseString

        count = 1
        while id in existingIds:
            id = "{0}_{1}".format(baseString, count)
            count += 1

        return id

    def convert(self, bmodel):
        """
        Convert BioOpt model to libSBML document instance

        :param bmodel: BioOpt model of type :class;`model.Model`
        :rtype: :class:`libsbml.SBMLDocument`
        """
        import libsbml

        doc = libsbml.SBMLDocument(self.level, self.version)
        model = doc.createModel()

        unit_definition = model.createUnitDefinition()
        unit_definition.setId("mmol_per_gDW_per_hr")
        mole = unit_definition.createUnit()
        mole.setKind(libsbml.UNIT_KIND_MOLE)
        mole.setScale(-3)

        second = unit_definition.createUnit()
        second.setKind(libsbml.UNIT_KIND_SECOND)
        second.setMultiplier(0.00027778)
        second.setExponent(-1)

        c_dict = {}
        mc_dict = {}
        m_dict = {}
        r_dict = {}

        class IdMap:
            def __init__(self, id, name, short=None):
                self.id = id
                self.name = name
                self.short = short if short else name

        # Find all compartments
        metabolites = bmodel.find_metabolites()
        for i, m in enumerate(metabolites, start=1):
            if self.compartment_pattern:
                c_pattern_res = self.compartment_pattern.search(m.name)
                if c_pattern_res:
                    c_name = c_pattern_res.group(1)
                else:
                    raise ValueError("Metabolite '{0}' doesn't match compartment pattern".format(m.name))
            else:
                c_name = "cell"

            if c_name not in c_dict:
                c_id = self.__get_valid_sbml_id(
                    "C_" + ("{0:04d}".format(len(c_dict) + 1) if self.compartment_id == "auto" else c_name), c_dict.keys())
                c_dict[c_name] = IdMap(c_id, c_name)
                compartment = model.createCompartment()
                compartment.setId(c_dict[c_name].id)
                compartment.setName(c_dict[c_name].name)

            mc_dict[m.name] = c_dict[c_name]

        if "boundary" not in c_dict:
            c_dict["boundary"] = IdMap("boundary", "boundary")

        # Assign abbreviations to compartments
        for c_name, c in c_dict.iteritems():
            for i in xrange(1, (len(c_name)-1)):
                short = c_name[0:(0+i)]
                if short not in (v.short for v in c_dict.itervalues()):
                    c.short = short
                    break

        mids_set = set()
        for i, m in enumerate(metabolites, start=1):
            m_id = "M_" + ("{0:04d}".format(i) if self.metabolite_id == "auto" else m.name)
            m_suffix = "_" + (c_dict["boundary"].short if m.boundary else mc_dict[m.name].short)
            abbr = self.compartment_suffix and not m_id.endswith(m_suffix)
            if abbr:
                m_id += m_suffix

            m_id = self.__get_valid_sbml_id(m_id, mids_set)
            mids_set.add(m_id)

            m_dict[m.name] = IdMap(m_id, m.name)
            species = model.createSpecies()
            species.setId(m_dict[m.name].id)
            species.setName(m_dict[m.name].name)
            species.setBoundaryCondition(m.boundary)
            species.setCompartment(mc_dict[m.name].id)
            species.setInitialAmount(0)

        for i, r in enumerate(bmodel.reactions, start=1):
            r_id = self.__get_valid_sbml_id("R_" + ("{0:04d}".format(i) if self.reaction_id == "auto" else r.name),
                                         r_dict.keys())
            r_dict[r.name] = IdMap(r_id, r.name)

            reaction = model.createReaction()
            reaction.setId(r_dict[r.name].id)
            reaction.setName(r.name)
            reaction.setCompartment(compartment.getId())
            reaction.setReversible(r.direction == Direction.reversible())

            for rm in r.reactants:
                reactant = reaction.createReactant()
                reactant.setSpecies(m_dict[rm.metabolite.name].id)
                reactant.setStoichiometry(rm.coefficient)

            for rm in r.products:
                product = reaction.createProduct()
                product.setSpecies(m_dict[rm.metabolite.name].id)
                product.setStoichiometry(rm.coefficient)

            law = reaction.createKineticLaw()

            ast_flux = libsbml.ASTNode(libsbml.AST_NAME)
            ast_flux.setName("FLUX_VALUE")
            law.setMath(ast_flux)

            r_lb = r.bounds.lb if abs(r.bounds.lb) != Bounds.inf() else math.copysign(self.inf, r.bounds.lb)
            r_ub = r.bounds.ub if abs(r.bounds.ub) != Bounds.inf() else math.copysign(self.inf, r.bounds.ub)
            # print "{0}: [{1}, {2}]".format(r_id, r_lb, r_ub)

            lower_bound = law.createParameter()
            lower_bound.setId("LOWER_BOUND")
            lower_bound.setUnits("mmol_per_gDW_per_hr")
            lower_bound.setValue(r_lb)

            upper_bound = law.createParameter()
            upper_bound.setId("UPPER_BOUND")
            upper_bound.setUnits("mmol_per_gDW_per_hr")
            upper_bound.setValue(r_ub)

            objective = law.createParameter()
            objective.setId("OBJECTIVE_COEFFICIENT")
            objective.setUnits("dimensionless")
            objective.setValue(0)
            o = bmodel.objective is not None and bmodel.objective.operands is not None and len(bmodel.objective.operands) and \
                r in bmodel.objective.operands
            objective.setValue(int(o))

            flux = law.createParameter()
            flux.setId("FLUX_VALUE")
            flux.setUnits("mmol_per_gDW_per_hr")
            flux.setValue(0)

        return doc