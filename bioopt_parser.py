import re
from model import *
import warnings

class BiooptParseWarning(Warning):
    pass

class BiooptParser(object):
    def __init__(self):
        # TODO: replace number with float() for performance reasons
        self.re_number = r"(?:[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)|(?:[-+]?(?:[0-9]*\.[0-9]+|[0-9]+))"
        self.re_member = r"(\(?(" + self.re_number + r") *\)? +)?(.*)"

    def parse_file(self, path):
        f = open(path, "r")
        return self.__parse(f.read(), filename=path)

    def parse_reactions_section(self, section_text):
        return list(r for r, i in self.__parse_reactions_section(section_text))

    def __parse_reactions_section(self, section_text, filename=None, section_start=0):
        if not isinstance(section_text, str):
            raise TypeError("Reactions section text is not of type string")

        return self.__parse_section(section_text, self.parse_reaction, filename=filename, section_start=section_start)

    def parse_reaction_member(self, member_str):
        if not isinstance(member_str, str):
            raise TypeError("Reaction member string is not of type string")

        member_str, multiline_comment = self.strip_comments(member_str, False)
        if not len(member_str):
            raise ValueError("Reaction member string is empty")

        m = re.match(self.re_member, member_str)
        if not m:
            raise SyntaxError("Could not parse reaction member: {0}".format(member_str))

        tmp, coef, name = m.groups()
        coef = float(coef) if coef else 1

        return ReactionMember(Metabolite(name), coef)

    def parse_reaction_member_list(self, list_str):
        if not isinstance(list_str, str):
            raise TypeError("Reaction member list string is not of type string")

        list_str, multiline_comment = self.strip_comments(list_str, False)
        if not len(list_str):
            raise ValueError("Reaction member list string is empty")

        parts = re.split(r"\s+\+\s+", list_str)
        members = ReactionMemberList([self.parse_reaction_member(s) for s in parts])

        return members

    def parse_reaction(self, line):
        if not isinstance(line, str):
            raise TypeError("Reaction line was not a string")

        line, multiline_comment = self.strip_comments(line, False)
        if not len(line):
            raise ValueError("Reaction string is empty")

        sep = ":"
        parts = line.split(sep)

        if len(parts) > 2:
            raise SyntaxError("{0} separator split reaction line into more than two parts [{1}]".format(sep, line))
        if len(parts) < 2:
            raise SyntaxError("Could not split reaction line using {0} separator [{1}]".format(sep, line))

        reaction_name = parts[0].strip()

        d = re.search("(\s+<\->|<\-|\->\s+)", line)
        direction = d.groups()[0]
        parts = parts[1].split(direction)
        direction = direction.strip() #removing white space after splitting into parts

        if len(parts) != 2:
            raise SyntaxError("Reaction doesn't consist of exactly two parts (reactants & products)")

        reactants = self.parse_reaction_member_list(parts[0])
        products = self.parse_reaction_member_list(parts[1])

        if direction == "<-":
            return Reaction(reaction_name, products, reactants, Direction.forward())
        elif direction == "->":
            return Reaction(reaction_name, reactants, products, Direction.forward())
        elif direction == "<->":
            return Reaction(reaction_name, reactants, products, Direction.reversible())
        else:
            raise Exception("Unknown direction ({0})".format(direction))

    # TODO: Add "%" commens functionality
    def strip_comments(self, line, multiline_comment=False):
        short_comment = False
        output_line = ""
        for l in line:
            if multiline_comment and l == "%":
                multiline_comment = False
                continue

            if multiline_comment or short_comment:
                continue

            if l == "%":
                multiline_comment = True
                continue

            if l == "#":
                short_comment = True
                continue

            output_line += l

        output_line = output_line.strip()

        return output_line, multiline_comment

    def parse_constraint(self, constraint_text):
        constraint_text, multiline_comment = self.strip_comments(constraint_text, False)

        re_bounds = "\[\s*(" + self.re_number + r")\s*,\s*(" + self.re_number + r")\s*\]"
        re_constraint = r"(.*)\s*" + re_bounds
        m = re.match(re_constraint, constraint_text)

        if m is None:
            raise SyntaxError("Could parse reaction constraint: {0}".format(constraint_text))

        reaction_name, lb, ub = m.groups()
        reaction_name = reaction_name.strip()
        bounds = Bounds(float(lb), float(ub))

        return Reaction(reaction_name, ReactionMemberList(), ReactionMemberList(), bounds.direction, bounds)

    def parse_constraints_section(self, section_text):
        return list(c for c, i in self.__parse_constraints_section(section_text))

    def __parse_constraints_section(self, section_text, filename=None, section_start=0):
        if not isinstance(section_text, str):
            raise TypeError("External metabolites section text is not of type string")

        return self.__parse_section(section_text, self.parse_constraint, filename=filename, section_start=section_start)

    def parse_objective_section_line(self, objective_line):
        parts = re.split(r"\s+", objective_line)
        if not parts:
            raise SyntaxError("Could not parse objective section line: {0}".format(objective_line))

        parsed_parts = []
        for i, p in enumerate(parts):
            try:
                p = float(p)
            except ValueError:
                p = Reaction(p)

            parsed_parts.append(p)

        if len(parsed_parts) == 1:
            return MathExpression(None, [parsed_parts[0]])
        else:
            return MathExpression(Operation.multiplication(), parsed_parts)

    def parse_objective_section(self, section_text):
        return self.__parse_objective_section(section_text)

    def __parse_objective_section(self, section_text, filename=None, section_start=0, reactions=None, section_name="-DESIGN OBJECTIVE/-OBJECTIVE", reactions_section_name="-REACTIONS"):
        if not isinstance(section_text, str):
            raise TypeError("Objective section text is not of type string")

        add_operands = list(self.__parse_section(section_text, self.parse_objective_section_line, filename=filename, section_start=section_start))

        if not reactions is None and len(reactions) > 0 and len(add_operands) > 0:
            line_reactions = ((i, r) for expression, i in add_operands for r in expression.find_variables() if isinstance(r, Reaction))
            for i, r in line_reactions:
                if r.name not in reactions:
                    warnings.warn_explicit(
                        "Reaction '{0}' from '{1}' section is not present in '{2}' section".format(r.name, section_name, reactions_section_name),
                        BiooptParseWarning, filename=filename, lineno=section_start+i+1)

        if len(add_operands) == 1:
            return add_operands[0][0]
        elif len(add_operands) > 1:
            return MathExpression(Operation.addition(), list(e[0] for e in add_operands))


    def parse_external_metabolites_section(self, section_text):
        return list(e for e, i in self.__parse_external_metabolites_section(section_text))

    def __parse_external_metabolites_section(self, section_text, filename=None, section_start=0):
        if not isinstance(section_text, str):
            raise TypeError("External metabolites section text is not of type string")

        return self.__parse_section(section_text, Metabolite, filename=filename, section_start=section_start)

    def find_sections(self, text):
        s = re.compile(r"^-[\w ]+$", re.MULTILINE)
        sections1 = [(m.group(), m.start(), m.end()) for m in re.finditer(s, text)]

        sections2 = dict()
        for i, s in enumerate(sections1):
            name = s[0]
            start = s[2]+1
            end = sections1[i+1][1]-1 if i+1 < len(sections1) else len(text)
            line = text[0:start].count("\n")

            sections2[name] = (start, end, line)

        return sections2

    def __parse_section(self, section_text, method, filename=None, section_start=0):
        nl = re.compile("\n\r|\r\n|\n")
        lines = nl.split(section_text)
        comment = False
        for i, line in enumerate(lines):
            line, multiline_comment = self.strip_comments(line, comment)
            # TODO: Return None to have errors with line information anntached
            if not len(line):
                continue

            saved_warnings = []
            with warnings.catch_warnings(record=True) as ws:
                warnings.simplefilter("always")
                result = method(line)

                for w in ws:
                    w_outer = warnings.WarningMessage(message=w.message, category=BiooptParseWarning, filename=filename, lineno=section_start+i+1, line=line)
                    saved_warnings.append(w_outer)

            for w in saved_warnings:
                warnings.warn_explicit(message=w.message, category=w.category, filename=w.filename, lineno=w.lineno)

            yield result, i

    def __find_section(self, text, sections, fun):
        for name in sections.keys():
            if fun(name):
                start, end, line = sections[name]
                return name, text[start:end], line

        return None, None, None

    def parse(self, text):
        return self.__parse(text)

    def __parse(self, text, filename=None):
        sections = self.find_sections(text)

        model = Model()
        react_name, react_text, react_line = self.__find_section(text, sections, lambda x: re.search(r"reac", x, re.I))
        const_name, const_text, const_line = self.__find_section(text, sections, lambda x: re.search(r"cons", x, re.I))
        ext_m_name, ext_m_text, ext_m_line = self.__find_section(text, sections, lambda x: re.search(r"ext", x, re.I))
        obj_name, obj_text, obj_line       = self.__find_section(text, sections, lambda x: re.search(r"obj", x, re.I) and not re.search("des", x, re.I))
        dobj_name, dobj_text, line         = self.__find_section(text, sections, lambda x: re.search(r"obj", x, re.I) and re.search("des", x, re.I))

        if react_text:
            model.reactions = list(r for r, i in self.__parse_reactions_section(react_text, filename=filename, section_start=react_line))
            model.unify_metabolite_references()
        else:
            warnings.warn("Could not find '-REACTIONS' section", BiooptParseWarning)

        reactions = dict()
        if model.reactions:
            reactions = dict((r.name, r) for r in model.reactions)

        if const_text:
            for c, i in self.__parse_constraints_section(const_text, filename=filename, section_start=const_line):
                if c.name in reactions and c.bounds.lb < 0 and reactions[c.name].direction != Direction.reversible():
                    warnings.warn_explicit(
                        "Reaction '{0}' from '{1}' has effective bounds not compatible with reaction direction in '{2}' section ({3} : {4})".format(c.name, const_name, react_name, reactions[c.name].direction, c.bounds),
                        BiooptParseWarning, filename=filename, lineno=const_line+i+1)

                if c.name in reactions:
                    reactions[c.name].bounds = c.bounds
                elif react_text:
                    warnings.warn_explicit(
                        "Reaction '{0}' from '{1}' section is not present in '{2}' section".format(c.name, const_name, react_name),
                        BiooptParseWarning, filename=filename, lineno=const_line+i+1)
        else:
            warnings.warn("Could not find '-CONSTRAINS' section", BiooptParseWarning)

        if ext_m_text:
            metabolites = dict((m.name, m) for m in model.find_metabolites())
            for m, i in self.__parse_external_metabolites_section(ext_m_text, filename=filename, section_start=ext_m_line):
                if m.name in metabolites:
                    metabolites[m.name].boundary = True
                elif react_text:
                    warnings.warn_explicit(
                        "Metabolite '{0}' from '{1}' section is not present in any reaction from '{2}' section".format(m.name, ext_m_name, react_name),
                        BiooptParseWarning, filename=filename, lineno=ext_m_line+i+1)
        else:
            warnings.warn("Could not find '-EXTERNAL METABOLITES' section", BiooptParseWarning)

        if obj_text:
            objective = self.__parse_objective_section(obj_text, section_name=obj_name, reactions_section_name=react_name, filename=filename, section_start=obj_line, reactions=reactions)
            model.objective = objective
        else:
            warnings.warn("Could not find '-OBJECTIVE' section", BiooptParseWarning)

        if dobj_text:
            design_objective = self.__parse_objective_section(dobj_text, section_name=dobj_name, filename=filename, section_start=line, reactions=reactions)
            model.design_objective = design_objective
        else:
            warnings.warn("Could not find '-DESIGN OBJECTIVE' section", BiooptParseWarning)


        model.unify_reaction_references()

        return model

if __name__ == "__main__":
    parser = BiooptParser()
    models = ["Respiro_ferm_prep.bioopt", "iAG612.bioopt", "X:/Joana/Bioopt models/Rich20/rich20MbarkeriiAF692.bioopt",
              "X:/Joana/Bioopt models/Rich20/richaeroEcoliWmodified.bioopt",
              "X:/Joana/Bioopt models/Rich20/richbacillus_Oh.bioopt",
              "X:/Joana/Bioopt models/Rich20/richchloroflexi.bioopt",
              "X:/Joana/Bioopt models/Rich20/richclostridiumbeijerinckiiANAER.bioopt",
              "X:/Joana/Bioopt models/Rich20/richEcoliiAF1260_ANAEROBIC.bioopt",
              "X:/Joana/Bioopt models/Rich20/richKLEBSIELLA.bioopt",
              "X:/Joana/Bioopt models/Rich20/richpseudomonasaerob.bioopt",
              "X:/Joana/Bioopt models/Rich20/richpyloriIT341aerobe.bioopt",
              "X:/Joana/Bioopt models/Rich20/richsalmonella.bioopt",
              "X:/Joana/Bioopt models/Rich20/richSaureusaerobiciSB619.bioopt",
              "X:/Joana/Bioopt models/Rich20/richshewanella.bioopt",
              "X:/Joana/Bioopt models/Rich20/richsynecho.bioopt",
              "X:/Joana/Bioopt models/Rich20/richthermotogaanaer.bioopt",
              "X:/Joana/Bioopt models/Rich20/richtuberculosisiNJ661aerobe.bioopt"]

    for path in models:
        print "Parsing '{0}'".format(path)
        model = parser.parse_file(path)
