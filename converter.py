from model import *
import re
from operator import itemgetter as _

class Bioopt2CobraPyConverter:
    """
    Bioopt model object to COBRApy model object converter. More details on COBRApy library can be found
    `online <https://github.com/opencobra/cobrapy>`_

    :rtype: :class:`Bioopt2COBRApyConverter`
    """

    def __init__(self, description=None):
        self.inf = 1000

        if not description or isinstance(description, str):
            self.description = description
        else:
            raise ValueError("Name '{0}' is not a string or None".format(description))

    def convert(self, bioopt_model):
        """
        Convert BioOpt model to COBRApy model instance

        :param bioopt_model: BioOpt model of type :class;`model.Model`
        :rtype: :class:`cobra.Model`
        """
        import cobra
        from cobra.core import DictList

        cobra_model = cobra.Model(self.description)

        # Precache converter metabolites for performance reasons !!!
        metabolites = {}
        for bm in bioopt_model.find_metabolites():
            if bm.boundary:
                continue

            cobra_metabolite = cobra.Metabolite(bm.name)
            cobra_metabolite._model = cobra_model
            metabolites[bm.name] = cobra_metabolite

        cobra_model.metabolites = DictList()
        cobra_model.metabolites.extend(metabolites.values())

        cobra_reactions = []
        for bioopt_reaction in bioopt_model.reactions:
            bounds = bioopt_reaction.find_effective_bounds()

            inf = Bounds.inf()
            cobra_reaction = cobra.Reaction(bioopt_reaction.name)
            cobra_reaction.lower_bound = bounds.lb if abs(bounds.lb) != inf else math.copysign(self.inf, bounds.lb)
            cobra_reaction.upper_bound = bounds.ub if abs(bounds.ub) != inf else math.copysign(self.inf, bounds.ub)
            cobra_reaction.objective_coefficient = bioopt_model.objective_dict.get(bioopt_reaction.name, 0)

            # Conversion is split into two loops for performance reasons !!!
            meta_dict = {}
            for member in bioopt_reaction.reactants:
                if member.metabolite.boundary:
                    continue

                cobra_metabolite = metabolites[member.metabolite.name]
                cobra_metabolite._reaction.add(cobra_reaction)
                meta_dict[cobra_metabolite] = -member.coefficient

            for member in bioopt_reaction.products:
                if member.metabolite.boundary:
                    continue

                cobra_metabolite = metabolites[member.metabolite.name]
                cobra_metabolite._reaction.add(cobra_reaction)
                meta_dict[cobra_metabolite] = member.coefficient

            cobra_reaction.add_metabolites(meta_dict, combine=False, add_to_container_model=False)
            cobra_reactions.append(cobra_reaction)

        #cobra_model.add_reactions(cobra_reactions)
        cobra_model.reactions = DictList()
        cobra_model.reactions.extend(cobra_reactions)

        return cobra_model


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
             metabolite_id="auto", compartment_id="auto", compartment_suffix=True, compartment_pattern=r"_(\w+)$",
             simulate_boundary=False, metabolite_map={}):

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

        if isinstance(metabolite_map, (dict)):
            self.metabolite_map = metabolite_map
        else:
            raise ValueError("Metabolite map is not a dictionary")

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

        self.simulate_boundary = simulate_boundary

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

    def convert(self, bioopt_model):
        """
        Convert BioOpt model to libSBML document instance

        :param bioopt_model: BioOpt model of type :class;`model.Model`
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


        compartment_names = {'e': "Extracellular", 'p': "Peroxisome", 'm': "Mitochondria", 'c': "Cytosol",
                             'l': "Lysosome", 'r': "Endoplasmic reticulum", 'g': "Golgi apparatus",
                             'n': "Nucleus", 'x': "Boundary"}

        # Find all compartments
        metabolites = sorted(bioopt_model.find_metabolites(), key=lambda x: x.name)
        for i, m in enumerate(metabolites, start=1):
            if self.compartment_pattern:
                c_pattern_res = self.compartment_pattern.search(m.name)
                if c_pattern_res:
                    c_name = c_pattern_res.group(1)
                else:
                    raise ValueError("Metabolite '{0}' doesn't match compartment pattern".format(m.name))
            else:
                c_name = "cell"

            if not c_name:
                pass

            if c_name not in c_dict:
                c_id = self.__get_valid_sbml_id(
                    "C_" + ("{0:04d}".format(len(c_dict) + 1) if self.compartment_id == "auto" else c_name), c_dict.keys())


                c_dict[c_name] = IdMap(c_id, compartment_names.get(c_name, c_name))
                compartment = model.createCompartment()
                compartment.setId(c_dict[c_name].id)
                compartment.setName(c_dict[c_name].name)

            mc_dict[m.name] = c_dict[c_name]

        if self.simulate_boundary and "boundary" not in c_dict:
            boundary_c = IdMap("boundary", "boundary")

        # Assign abbreviations to compartments
        for c_name, c in sorted(c_dict.iteritems(), key=_(0)):
            for i in xrange(1, (len(c_name)-1)):
                short = c_name[0]
                if short not in (v.short for v in c_dict.itervalues()):
                    c.short = short
                    break

        mids_set = set()
        for i, m in enumerate(metabolites, start=1):
            m_id = ("{0:04d}".format(i) if self.metabolite_id == "auto" else m.name)
            if not m_id.startswith("M_"):
                m_id = "M_" + m_id

            m_suffix = "_" + (boundary_c.short if m.boundary and self.simulate_boundary else mc_dict[m.name].short)
            abbr = self.compartment_suffix and not m_id.endswith(m_suffix)
            if abbr:
                m_id += m_suffix

            m_id = self.__get_valid_sbml_id(m_id, mids_set)
            mids_set.add(m_id)

            m_name2 = re.sub("xtX$", "", m.name)
            if m.name in self.metabolite_map:
                m_dict[m.name] = IdMap(m_id, self.metabolite_map[m.name])
            elif m_name2 in self.metabolite_map:
                m_dict[m.name] = IdMap(m_id, self.metabolite_map[m_name2])
            else:
                if self.metabolite_map:
                    print "{} couldn't be mapped to metabolite name...".format(m.name)
                m_dict[m.name] = IdMap(m_id, m.name)

            species = model.createSpecies()
            species.setId(m_dict[m.name].id)
            species.setName(m_dict[m.name].name)
            species.setBoundaryCondition(m.boundary)
            species.setCompartment(mc_dict[m.name].id)
            species.setInitialAmount(0)

        for i, r in enumerate(bioopt_model.reactions, start=1):
            r_id = self.__get_valid_sbml_id(("{0:04d}".format(i) if self.reaction_id == "auto" else r.name), r_dict.keys())
            if not r_id.startswith("R_"):
                r_id = "R_" + r_id

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
            o = bioopt_model.objective is not None and bioopt_model.objective.operands is not None \
                and len(bioopt_model.objective.operands) and r in bioopt_model.objective.operands
            objective.setValue(int(o))

            flux = law.createParameter()
            flux.setId("FLUX_VALUE")
            flux.setUnits("mmol_per_gDW_per_hr")
            flux.setValue(0)

        return doc

