import re
from model import *

class BiooptParser(object):
    def __init__(self):
        # TODO: replace number with float() for performance reasons
        self.re_number = r"(?:[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)|(?:[-+]?(?:[0-9]*\.[0-9]+|[0-9]+))"
        self.re_member = r"(\(?(" + self.re_number + r") *\)? +)?(.*)"

    def parse_file(self, path):
        f = open(path, "r")
        self.parse(f.read())

    def parse_reactions_section(self, section_text):
        if not isinstance(section_text, str):
            raise TypeError("Reactions section text is not of type string")

        return list(self.__parse_section(section_text, self.parse_reaction))

    def parse_reaction_member(self, member_str):
        if not isinstance(member_str, str):
            raise TypeError("Reaction member string is not of type string")

        member_str = self.strip_comments(member_str)
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

        list_str = self.strip_comments(list_str)
        if not len(list_str):
            raise ValueError("Reaction member list string is empty")

        parts = re.split(r"\s+\+\s+", list_str)
        members = ReactionMemberList([self.parse_reaction_member(s) for s in parts])

        return members

    def parse_reaction(self, line):
        if not isinstance(line, str):
            raise TypeError("Reaction line was not a string")

        line = self.strip_comments(line)
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

    def strip_comments(self, line):
        line = line.strip()
        line = re.sub("\s*[#\n\r].*", "", line)

        return line

    def parse_constraint(self, constraint_text):
        constraint_text = self.strip_comments(constraint_text)

        re_bounds = "\[\s*(" + self.re_number + r")\s*,\s*(" + self.re_number + r")\s*\]"
        re_constraint = r"(.*)\s*" + re_bounds
        m = re.match(re_constraint, constraint_text)

        if m is None:
            raise SyntaxError("Could parse reaction constraint: {0}".format(constraint_text))

        reaction_name, lb, ub =  m.groups()
        reaction_name = reaction_name.strip()
        bounds = Bounds(float(lb), float(ub))

        return Reaction(reaction_name, ReactionMemberList(), ReactionMemberList(), bounds.direction, bounds)

    def parse_constraints_section(self, section_text):
        if not isinstance(section_text, str):
            raise TypeError("External metabolites section text is not of type string")

        return list(self.__parse_section(section_text, self.parse_constraint))

    def __parse_objective_section_line(self, objective_line):
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
            # TODO: test for operation None
            return MathExpression(None, [parsed_parts[0]])
        else:
            return MathExpression(Operation.multiplication(), parsed_parts)

    def parse_objective_section(self, section_text):
        if not isinstance(section_text, str):
            raise TypeError("Objective section text is not of type string")

        expression = None
        for exp in self.__parse_section(section_text, self.__parse_objective_section_line):
            if expression is None:
                expression = exp
            else:
                expression = MathExpression(expression, exp, Operation.addition())

        return expression

    def parse_external_metabolites_section(self, section_text):
        if not isinstance(section_text, str):
            raise TypeError("External metabolites section text is not of type string")

        return list(self.__parse_section(section_text, Metabolite))

    def find_sections(self, text):
        s = re.compile(r"^-[\w ]+$", re.MULTILINE)
        sections1 = [(m.group(), m.start(), m.end()) for m in re.finditer(s, text)]

        sections2 = dict()
        for i, s in enumerate(sections1):
            name = s[0]
            start = s[2]+1
            end = sections1[i+1][1]-1 if i+1 < len(sections1) else len(text)

            sections2[name] = (start, end)

        return sections2

    def __parse_section(self, section_text, method):
        nl = re.compile("\n\r|\r\n|\n")
        lines = nl.split(section_text)
        for i, line in enumerate(lines):
            line = self.strip_comments(line)
            if not len(line):
                continue

            yield method(line)

    def parse(self, text):
        sections = self.find_sections(text)

        # TODO: print line where error appeared
        model = Model()

        react_section = sections["-REACTIONS"]
        model.reactions = self.parse_reactions_section(text[react_section[0]:react_section[1]])

        # TODO: test for duplicate reactions
        reactions = dict((r.name, r) for r in model.reactions)

        const_section = sections["-CONSTRAINTS"]
        constrains = dict((reaction.name, reaction) for reaction in self.parse_constraints_section(text[const_section[0]:const_section[1]]))
        for r in model.reactions:
            r.bounds = constrains[r.name].bounds

        model.unify_metabolite_references()

        ext_m_section = sections["-EXTERNAL METABOLITES"]
        external_metabolites = [m.name for m in self.parse_external_metabolites_section(text[ext_m_section[0]:ext_m_section[1]])]
        for m in model.find_metabolites():
            m.boundary = m.name in external_metabolites

        obj_section = sections["-OBJ"]
        objective = self.parse_objective_section(text[obj_section[0]:obj_section[1]])
        model.objective = objective

        dobj_section = sections["-DESIGNOBJ"]
        design_objective = self.parse_objective_section(text[dobj_section[0]:dobj_section[1]])
        model.design_objective = design_objective

        model.unify_reaction_references()

        return model

if __name__ == "__main__":
    parser = BiooptParser()
    parser.parse_file("Respiro_ferm_prep.bioopt")
