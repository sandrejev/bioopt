#!/usr/bin/env python
import sys,re
from sets import Set
import dis



path="./Scere_master_Kiran.txt"

def error(code):
    error_map = {1: "Check tags order: -REACTIONS -CONSTRAINTS -EXTERNAL METABOLITES -OBJ -DESIGNOBJ",
                2 : "Check format, seems that -REACTIONS is not found of after -CONSTRAINTS",
                3 : "Error3",
                4 : "Direction can be only f, b , r" }
    if code in error_map.keys():
        return error_map[code]
    else:
        return "wrong code"



class Metabolite:
    def __init__(self, name, coefficient):
        if not isinstance(name, basestring):
            raise TypeError("Name is not a string")

        if not isinstance(coefficient, (int, float)):
            raise TypeError("Coefficient is not a number")

        coefficient = float(coefficient)

        if coefficient <= 0:
            raise ValueError("Coefficient value is zero or negative")

        self.__name = name
        self.__coefficient = coefficient

    def __eq__(self, other):
        return self.get_name() == other.get_name()

    def get_name(self):
        return self.__name
    def get_coefficient(self):
        return self.__coefficient

class Direction:
    forward = 0
    reversible = 1

    def is_valid(self, d):
        return d in [self.forward, self.reversible]

class Bounds:
    def ___init__
    lb =


class Reaction:
    def __init__(self, name = None, SubstrateList = None, ProductList = None, direction = None, bounds = None):
        """Creates object representing a reaction

        Args:
            name (str): Reaction name
            Substrates (Metabolite): List of metabolites with coeficients
        """
        if not isinstance(name, basestring):
            raise TypeError("Name is not a string")

        if not isinstance(SubstrateList, list):
            raise TypeError("SubstrateList is not a list")

        if not isinstance(SubstrateList, list):
            raise TypeError("SubstrateList is not a list")

        if not Direction.is_valid(direction):
            raise ValueError("Wrong direction were specified, can only be f,r")

        if not isinstance(bounds, tuple):
            raise TypeError("Bounds is not a tuple")

        self.__name       = name          if name          else ""
        self.__substrates = SubstrateList if SubstrateList else []
        self.__products   = ProductList   if ProductList   else []
        self.__direction  = direction     if direction     else Direction.reversible

        if bounds is None:
            if self.__direction == Direction.forward:
                self.__lb = 0
                self.__ub = float("inf")
            elif self.__direction == Direction.reversible:
                self.__lb = float("-inf")
                self.__ub = float("inf")
        else:



        self.check();

    def check(self):
        if bool(self.__substrates.get_all()) and bool(self.__products.get_all()):
            a = Set(self.__substrates.get_all().keys())
            b = Set(self.__products.get_all().keys())
            if a.intersection(b):
                print "Error: cannot be same metabolites in substrates/products, reaction: {reaction_name}".format(reaction_name=self.__name)


    #SET methods
    def set_name(self, name):
        self.__name = name
    def set_direction(self, direction):
        dirs = ['f','b','r']
        if direction in dirs:
            self.__direction = direction
        else:
            print "Check direction!","current:", direction
            exit(-1)
    #GET methods
    def get_name(self):
        return self.__name

    def get_substrates(self):
        return self.__substrates.get_all()

    def get_products(self):
        return self.__products.get_all()
    def get_direction(self):
        return self.__direction




def make_reaction (line, sep = ":"):
    line = line.strip()
    fails = 0

    parts = line.split(sep)

    if len(parts) != 2:
        print make_reaction.__name__ , "error: delimiter",sep, line
        fails += 1
        exit(-1)

    reaction_name = parts[0].strip()
    d = re.search("(\s+<\->|<\-|\->\s+)", line)
    direction = d.groups()[0]
    parts = parts[1].split(direction)
    #removing white space after splitting into parts
    direction = direction.strip()

    if len(parts) != 2:
        print make_reaction.__name__ , "error: direction",direction, line
        fails += 1
        exit(-1)

    left_part = parts[0]
    right_part = parts[1]

    substrates_part = re.split(r"\s+\+\s+",left_part)
    products_part   = re.split(r"\s+\+\s+",right_part)

    re_number = r"(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?"

    number = re.compile("[\(]?\s?"+"("+re_number+")"+"\s?[\)]?\s+")
    substrates = {}
    for metabolite in substrates_part:
        metabolite = metabolite.strip();
        coefficient = 1
        n = re.search(number, metabolite)
        if n:
            coefficient = float(n.groups()[0])
        metabolite = re.sub(number, "", metabolite)
        substrates[metabolite] = coefficient

    products = {}
    for metabolite in products_part:
        metabolite = metabolite.strip();
        coefficient = 1
        n = re.search(number, metabolite)
        if n:
            coefficient = float(n.groups()[0])
        metabolite = re.sub(number, "", metabolite)
        products[metabolite] = coefficient

    if direction == "->":
        direction = "f"
    elif direction == "<->":
        direction = "r"

    if direction == "<-":
        return reaction_name, products, substrates, "f"
    else:
        return reaction_name, substrates, products, direction

def parse_file(path, force = 1):
    try:
        f = open(path)
    except IOError as e:
        print "I/O error({0}): {1}".format(e.errno, e.strerror)
        exit(-1)
    parse(f, force)


def parse(f, force=1):

    skip_flag     = 0
    reactions_flag = 0
    constraints_flag = 0
    external_metabolites_flag = 0
    obj_flag = 0
    designobj_flag = 0

    fails = 0
    minor_fails = 0

    reactions = {}
    constraints = {}
    external_metabolites = {}
    obj = {}
    design_obj = {}


    for nr, line in enumerate(f):
        line = line.strip()
        c = nr + 1

        #checking flags
        if not line:
            continue
        if re.match('^#',line):
            continue
        if re.match('^%',line) and skip_flag == 0:
            skip_flag = 1
        elif re.match('^%',line) and skip_flag == 1:
            skip_flag = 0
            continue
        if skip_flag == 1:
            continue

        if re.match('^-REACTIONS', line):
            reactions_flag = 1
            continue
        if re.match('^-CONSTRAINTS', line) and reactions_flag == 1:
            reactions_flag = 0
            constraints_flag = 1
            continue

        elif re.match('^-CONSTRAINTS', line) and reactions_flag == 0:
            print parse.__name__, "at line {}, check if -REACTIONS flag was defined before. {}".format(c, error(1))
            fails += 1
            continue

        if re.match('^-EXTERNAL METABOLITES', line) and constraints_flag == 1:
            constraints_flag = 0
            external_metabolites_flag = 1
            continue
        elif re.match('^-EXTERNAL METABOLITES', line) and constraints_flag == 0:
            print parse.__name__, "at line {}, check if -CONSTRAINTS flag was defined before. {}".format(c, error(1))
            fails += 1
            continue

        if re.match('^-OBJ', line) and external_metabolites_flag == 1:
            external_metabolites_flag = 0
            obj_flag = 1
            continue
        elif re.match('^-OBJ', line) and external_metabolites_flag == 0:
            print parse.__name__, "at line {}, check if -EXTERNAL METABOLITES flag was defined before. {}".format(c, error(1))
            fails += 1
            continue

        if re.match('^-DESIGNOBJ', line) and obj_flag == 1:
            obj_flag = 0
            designobj_flag = 1
            continue
        elif re.match('^-OBJ', line) and external_metabolites_flag == 0:
            print parse.__name__, "at line {}, check if -OBJ flag was defined before. {}".format(c, error(1))
            fails += 1
            continue

        #colecting data
        if reactions_flag == 1:
            #print c, line
            name, substrates, products, direction = make_reaction(line)
            if name in reactions:
                print "Such reaction {} exists at {}".format(name,c)
                fails += 1
            else:
                reactions[name] = [substrates,products,direction]
            continue

        if constraints_flag == 1:
            re_number = r"(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?"
            m = re.match(r'(\w+)\s?\[\s?'+'('+re_number+')'+'\s?,\s?'+'('+re_number+')'+'\s?\]', line)
            if m:
                #re_number if matched returns 3 elements tuple
                #m will contain total 7 elements in 0th element will be name +3 for lb and 3 for ub
                name = m.groups()[0]
                lb   = m.groups()[1]
                ub   = m.groups()[4]
                if name in reactions:
                    constraints[name] = [lb,ub]
                else:
                    print "No such reaction {} for constraint {} at line {}".format(name,line,c)
            else:
                print "Wrong format at {}, expected constraints format: reaction[number, number]".format(c)
                fails +=1
            continue

        if external_metabolites_flag == 1:
            if line in external_metabolites:
                print "Duplicate found of external metabolite at {}".format(c)
                minor_fails += 1
            else:
                external_metabolites[line] = 1
            continue

        if obj_flag == 1:
            re_number = r"(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?"
            m = re.match(r'(\w+)\s+'+'('+re_number+')'+'\s+'+'('+re_number+')', line)
            if m:
                name = m.groups()[0]
                c1   = m.groups()[1]
                c2  = m.groups()[4]
                obj[name] = [c1,c2]
            else:
                print "Wrong format at {}, expected objective format: reaction number number".format(c)
                fails +=1
            continue
        if designobj_flag == 1:
            re_number = r"(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?"
            m = re.match(r'(\w+)\s+'+'('+re_number+')'+'\s+'+'('+re_number+')', line)
            if m:
                name = m.groups()[0]
                c1   = m.groups()[1]
                c2  = m.groups()[4]
                design_obj[name] = [c1,c2]
            else:
                print "Wrong format at {}, expected design objective format: reaction number number".format(c)
                fails +=1
            continue

    if force == 1:

        if fails:
            print "There were critical {} fails found and {}".format(fails, minor_fails)
            quit (-1)
        else:

            return reactions, constraints, external_metabolites, obj, design_obj
    else:
        return reactions, constraints, external_metabolites, obj, design_obj







if __name__ == "__main__":
    print "asdasdasd"






"""



a = parse_file(path)

print a


name, substrates, products, direction = make_reaction(r"reaction: 0.25 2'-bgdfg + 3445 z <-> 0.34354E10 2'-bgdfdfg ")

r = Reaction(name, Metabolites(substrates), Metabolites(products), direction)
"""