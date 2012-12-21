import re
from model import *
class Test(object):
    def run(self):
        return True

class BiooptParser(object):
    def parse_file(self, path, force = 1):
        f = open(path, "r")
        self.__parse(f.read(), force)

    def parse_reactions_section(self, section_text):
        print "Parse reactions"

        nl = re.compile("[\n\r]")
        rstrings = nl.split(section_text)
        for line in rstrings:
            yield self.parse_reaction(line)

    def parse_reaction_member(self, member_str):
        if not isinstance(member_str, str):
            raise TypeError("Reaction member string is not of type string")

        member_str = self.strip_comments(member_str)
        if not len(member_str):
            raise ValueError("Reaction member string is empty")

        re_number = r"(?:[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)|(?:[-+]?(?:[0-9]*\.[0-9]+|[0-9]+))"
        re_member = r"(\(?(" + re_number + r") *\)? +)?(.*)"

        m = re.match(re_member, member_str)
        if m is None:
            raise SyntaxError("Could not parse reaction member: {0}".format(member_str))

        tmp, coef, name = m.groups()
        coef = float(coef) if coef else 1

        return ReactionMember(Metabolite(name), coef)

    def parse_reaction_member_list(self, list_str):
        list_str_split = re.split(r"\s+\+\s+", list_str)
        members = ReactionMemberList()
        for member_str in list_str_split:
            members.append(self.parse_reaction_member(member_str))

        return members

    def parse_reaction(self, line):
        if line is None:
            return None

        if not isinstance(line, str):
            raise TypeError("Reaction line was not a string")

        line = self.strip_comments(line)
        if not len(line):
            return None

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
        line = re.sub("[#\n\r].*", "", line)

        return line

    def parse_constraints_section(self, section_text):
        print "Parse constraints"
        pass

    def parse_objective_section(self, section_text):
        print "Parse objective"
        pass

    def parse_external_metabolites_section(self, section_text):
        print "Parse external metabolites"
        pass

    def __find_sections(self, text):
        s = re.compile(r"^-[\w ]+$", re.MULTILINE)
        #return dict((m.group(), text.count("\n", 0, m.start())) for m in re.finditer(s, text))
        sections1 = [(m.group(), m.start(), m.end()) for m in re.finditer(s, text)]

        sections2 = dict()
        for i, s in enumerate(sections1):
            name = s[0]
            start = s[2]+1
            end = sections1[i+1][1]-1 if i+1 < len(sections1) else len(text)

            sections2[name] = (start, end)

        return sections2



    def __parse(self, text, force=1):
        #self.__remove_comments(text)
        sections = self.__find_sections(text)

        model = Model()

        react_section = sections["-REACTIONS"]
        reactions = self.parse_reactions_section(text[react_section[0]:react_section[1]])

        model.reactions.extend(reactions)

#        ext_m_section = sections["-EXTERNAL METABOLITES"]
#        self.parse_external_metabolites(text[ext_m_section[0]:ext_m_section[1]])
#
#        const_section = sections["-CONSTRAINTS"]
#        self.parse_constraints(text[const_section[0]:const_section[1]])
#
#        obj_section = sections["-OBJ"]
#        self.parse_objective(text[obj_section[0]:obj_section[1]])
