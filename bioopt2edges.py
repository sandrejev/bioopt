from bioopt_parser import *
import argparse
import cplex_utils
import cplex

def blocked(prob, reactions):
    """
    Find blocked reactions in constraint based model
    :param prob: Problem instance
    :param reactions: List of reactions to be checked for "blocked" condition
    :return: List of blocked reactions
    """
    m = prob.model

    # Constraint all reactions with [-1, 1] or [0, 1] range (reversible, irreversible)
    reactions_set = set(reactions)
    m.variables.set_lower_bounds([(r_i, -1 if prob.rxn2bounds[rxn].lb < 0 else 0) for rxn, r_i in prob.rxn2i.iteritems() if rxn in reactions_set])
    m.variables.set_upper_bounds([(r_i, 1 if prob.rxn2bounds[rxn].ub > 0 else 0) for rxn, r_i in prob.rxn2i.iteritems() if rxn in reactions_set])

    # Optimize sum of all reactions
    m.objective.set_linear([(r_i, 1) for r_i in xrange(len(prob.rxn2i))])

    # Create a list of candidates for being called blocked
    reactions_set = set(reactions)
    blocked_candidates_set = set(r_i for rxn, r_i in prob.rxn2i.iteritems() if rxn in reactions_set)
    blocked_candidates = list(blocked_candidates_set)

    # Maximize (and minimize) sum of all reactions to see which reactions can possibly have flux.
    # The ones that will have non-zero flux can not possibly be blocked. Repeat N times.
    for p in range(50):
        # MAXIMIZE CANDIDATES
        m.objective.set_linear([(r_i, 0) for r_i in xrange(len(prob.rxn2i))])
        m.objective.set_linear([(r_i, 1) for r_i in blocked_candidates])
        m.objective.set_sense(m.objective.sense.maximize)
        m.solve()
        active_fwd_set = set(r_i for r_i, flux in zip(blocked_candidates, m.solution.get_values(blocked_candidates)) if abs(flux) > 1e-10)

        # MINIMIZE CANDIDATES
        active_rev_set = set(r_i for r_i, flux in zip(blocked_candidates, m.solution.get_values(blocked_candidates)) if abs(flux) > 1e-10)
        m.objective.set_sense(m.objective.sense.minimize)
        m.solve()

        blocked_candidates_set = blocked_candidates_set - active_fwd_set.union(active_rev_set)
        blocked_candidates = list(blocked_candidates_set)
        print "Suspected blocked reactions: {}".format(len(blocked_candidates_set))

    # Now check blocked reactions one by one by setting objective coefficient for each reaction to 1 and then
    # maximizing and minimizing that reaction
    blocked = []
    m.objective.set_linear([(r_i, 0) for r_i in xrange(len(prob.rxn2i))])
    for i, r_i in enumerate(blocked_candidates):
        print "{}/{}: {}".format(i, len(blocked_candidates), prob.i2rxn[r_i]),

        v_min, v_max = 0.0, 0.0

        # Set objective coefficient to 1
        m.objective.set_linear(r_i, 1)
        m.objective.set_sense(m.objective.sense.maximize)
        prob.model.solve()
        if cplex_utils.is_optimal(prob.model):
            v_max = m.solution.get_objective_value()

        # If reaction can't have positive flux try to minimize the same reaction
        if v_max == 0:
            m.objective.set_sense(m.objective.sense.minimize)
            m.solve()
            if cplex_utils.is_optimal(prob.model):
                v_min = m.solution.get_objective_value()


        # If reaction don't have flux it is blocked!
        if abs(v_min) < 1e-10 and abs(v_max) < 1e-10:
            blocked.append(r_i)
            print ""
        else:
            print "[{} {}]".format(v_min, v_max)

        # Reset objective equation
        m.objective.set_linear(r_i, 0)

    return blocked


def bin(val):
    """
    Convert boolean value to 0/1 string
    :param val: Boolean input value
    :return: String representing input value
    """
    return "1" if val else "0"

class Edge:
    def __init__(self, source, destination):
        self.source = source
        self.destination = destination
        self.constrained = False
        self.blocked = False # Is reaction in this edge blocked
        self.removed = False # Was any metabolite/reaction in this edge specifically removed by regex filter
        self.genesis = False # Can reaction in this edge have flux even with [0,0] uptake constraints

class Node:
    def __init__(self, name, type):
        self.name = name
        self.type = type # 'reaction' or 'metabolite'
        self.constrained = False # Is this reaction blocked (only for type='reaction)
        self.blocked = False # Is this reaction blocked (only for type='reaction)
        self.removed = False # Was this metabolite/reaction specificaly removed by regex filter
        self.genesis = False # Can this reaction (only for type='reaction) have flux even with [0,0] uptake constraints

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Converts moel from bioopt format to optflux format')
    parser.add_argument('bioopt', action='store', help='File containing bioopt model')
    parser.add_argument('output', action='store', help='Output file prefix (.edges and .nodes files are created')
    parser.add_argument('--metabolite-map', dest="metabolite_map", action='store', help="Map metabolite identifiers to names")
    parser.add_argument('--remove-nodes', dest="remove_nodes", action='store', help="Comma separated list of nodes not to be included in the final edges list")

    args = parser.parse_args()

    #
    # Read metabolite_id => metabolite_name conversion file. The file must have at least two columns (first: id, second: name)
    #
    metabolite_map = {}
    if args.metabolite_map:
        for i, line in enumerate(open(args.metabolite_map)):
            l = re.split("\t", line.rstrip())
            if len(l) >= 2:
                metabolite_map[l[0]] = l[1]
            else:
                raise RuntimeError("Error on line {} in '{}' file. Map file should have 2 or more columns".format(i, args.metabolite_map))

    remove_nodes = re.compile(args.remove_nodes) if args.remove_nodes else None

    # Read bioopt model
    parser = BiooptParser()
    model = parser.parse_file(args.bioopt)

    # Find boundary reactions
    model_sinks = set(r.name for r in model.reactions if r.products[0].metabolite.boundary or r.reactants[0].metabolite.boundary)

    edges = []
    nodes = {}
    for r in model.find_reactions():
        r_node = Node(r.name, "reaction")
        r_node.constrained = r.bounds.lb == 0 and r.bounds.ub == 0
        r_node.blocked = r_node.constrained
        r_node.removed = remove_nodes and remove_nodes.match(r_node.name)
        nodes[r.name] = r_node

        # Create Metabolite->Reaction edges
        for reactant in r.reactants:
            m_name = reactant.metabolite.name
            m_name2 = metabolite_map.get(m_name, m_name)
            m_node = nodes.get(m_name, Node(m_name2, "metabolite"))
            if m_name not in nodes:
                nodes[m_name] = m_node

            edge = Edge(m_node, r_node)
            edge.constrained = r_node.constrained
            edge.blocked = r_node.constrained
            edge.removed = remove_nodes and (remove_nodes.match(m_name2) or remove_nodes.match(m_name) or remove_nodes.match(r.name))
            edges.append(edge)

        # Create Reaction->Metabolite edges
        for product in r.products:
            m_name = product.metabolite.name
            m_name2 = metabolite_map.get(m_name, m_name)
            m_node = nodes.get(m_name, Node(m_name2, "metabolite"))
            if m_name not in nodes:
                nodes[m_name] = m_node

            edge = Edge(r_node, m_node)
            edge.constrained = r_node.constrained
            edge.blocked = r_node.constrained
            edge.removed = remove_nodes and (remove_nodes.match(m_name2) or remove_nodes.match(m_name) or remove_nodes.match(r.name))
            edges.append(edge)

    prob = cplex_utils.bioopt2cplex(model)

    #
    # Find blocked reactions
    #
    prob.model.variables.set_lower_bounds([(r_i, -1 if prob.rxn2bounds[rxn].lb < 0 else 0) for rxn, r_i in prob.rxn2i.iteritems() if rxn in model_sinks])
    prob.model.variables.set_upper_bounds([(r_i, 1 if prob.rxn2bounds[rxn].ub > 0 else 0) for rxn, r_i in prob.rxn2i.iteritems() if rxn in model_sinks])
    potential_blocked_reactions = [n.name for n in nodes.itervalues() if n.type == "reaction" and not n.blocked and not n.constrained]
    blocked_reactions = set(prob.i2rxn[r_i] for r_i in blocked(prob, potential_blocked_reactions))
    #blocked_reactions = set(r.strip() for r in open("/g/patil/Sergej/CancerHeterogeneity/data/model/blocked_reactions.txt"))
    for n in nodes.itervalues():
        if n.type == "reaction":
            n.blocked = n.name in blocked_reactions

    for e in edges:
        e.blocked = e.source.blocked or e.destination.blocked


    #
    # Find infinite flux
    #
    prob.model.variables.set_lower_bounds([(r_i, 0) for rxn, r_i in prob.rxn2i.iteritems() if rxn in model_sinks])
    prob.model.variables.set_upper_bounds([(r_i, 0) for rxn, r_i in prob.rxn2i.iteritems() if rxn in model_sinks])
    potential_gen_reactions = set(n.name for n in nodes.itervalues() if n.type == "reaction" and not n.blocked and not n.constrained)
    genesis_reactions = potential_gen_reactions - set(prob.i2rxn[r_i] for r_i in blocked(prob, list(potential_gen_reactions)))
    #genesis_reactions = set(r.strip() for r in open("/g/patil/Sergej/CancerHeterogeneity/data/model/genesis_reactions.txt"))
    for n in nodes.itervalues():
        if n.type == "reaction":
            n.genesis = n.name in genesis_reactions

    for e in edges:
        e.genesis = e.source.genesis or e.destination.genesis


    #
    # Find infinie metabolites (unfinished)
    #
    if False:
        prob.model.variables.set_lower_bounds([(prob.rxn2i[r.name], -1 if r.bounds.lb < 0 else 0) for r in model.reactions])
        prob.model.variables.set_upper_bounds([(prob.rxn2i[r.name], 1 if r.bounds.ub > 0 else 0) for r in model.reactions])
        prob.model.variables.set_lower_bounds([(r_i, 0) for rxn, r_i in prob.rxn2i.iteritems() if rxn in model_sinks])
        prob.model.variables.set_upper_bounds([(r_i, 0) for rxn, r_i in prob.rxn2i.iteritems() if rxn in model_sinks])
        prob.model.objective.set_linear([(r_i, 0) for r_i in prob.i2rxn])

        for c in model.find_metabolites():
            c_i = prob.cpd2i[c.name]
            prob.model.variables.add(lb=[0], ub=[1000], names=["MET_EXPORT_TEST"], columns=[cplex.SparsePair([c_i], [-1])], obj=[1])
            prob.model.solve()

            status = prob.model.solution.get_status_string()
            obj = prob.model.solution.get_objective_value() if status.startswith("optimal") else 0
            if obj > 1e-10:
                print "{} {}".format(c.name, obj)

            prob.model.variables.delete(prob.rxnnum)




    # Write edges file
    with open(args.output + ".edges", 'w') as f_output:
        f_output.writelines("source\tdestination\tconstrained\tremoved\tblocked\tinfinite_flux\n")
        for e in edges:
            f_output.writelines("{}\t{}\t{}\t{}\t{}\t{}\n".format(e.source.name, e.destination.name, bin(e.constrained), bin(e.removed), bin(e.blocked), bin(n.genesis)))

        print "Created {}.edges file with {} edges".format(args.output + ".edges", len(edges))

    # Write nodes file
    with open(args.output + ".nodes", 'w') as f_output:
        f_output.writelines("node\ttype\tid\tconstrained\tremoved\tblocked\tinfinite_flux\n")
        for id, n in nodes.iteritems():
            f_output.writelines("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(n.name, n.type, id, bin(n.constrained), bin(n.removed), bin(n.blocked), bin(n.genesis)))

        print "Created {}.nodes file with {} nodes".format(args.output + ".nodes", len(nodes))

