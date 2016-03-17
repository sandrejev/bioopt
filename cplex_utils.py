import cplex
import model
from tempfile import NamedTemporaryFile
import random
from itertools import chain


dir_fwd = model.Direction.forward()
dir_rev = model.Direction.reversible()


def is_positive(model):
    obj = model.solution.get_objective_value()
    return is_optimal(model) and obj > 1e-6


def is_nonzero(model):
    return is_optimal(model) and abs(model.solution.get_objective_value()) > 1e-6


def is_optimal(model):
    status = model.solution.get_status_string()
    return status == 'optimal'


def copy(model):
    with NamedTemporaryFile(delete=False) as f:
        params_path = f.name

    with NamedTemporaryFile(delete=False) as f:
        problem_path = f.name

    model.parameters.write_file(params_path + ".prm")
    model.write(problem_path + ".sav")

    model_copy = cplex.Cplex()
    model_copy.parameters.read_file(params_path)
    model_copy.read(problem_path + ".sav")

    return model_copy



class FbaBounds:
    def __init__(self, lb, ub):
        self.lb = lb
        self.ub = ub

class FbaProblem:
    def __init__(self, model, rxn2i, cpd2i, rxn2ext, rxn2cpd, rxn2dir, rxn2bounds, obj):
        self.model = model
        self.rxn2i = rxn2i
        self.rxn2i_external = {r_id: r_i for r_id, r_i in self.rxn2i.iteritems() if rxn2ext[r_id]}
        self.rxn2i_internal = {r_id: r_i for r_id, r_i in self.rxn2i.iteritems() if not rxn2ext[r_id]}
        self.rxn2ext = rxn2ext
        self.rxn2cpd = rxn2cpd
        self.cpd2i = cpd2i
        self.i2rxn = [x[0] for x in sorted(rxn2i.iteritems(), key=lambda x: x[1])]
        self.i2cpd = [x[0] for x in sorted(cpd2i.iteritems(), key=lambda x: x[1])]
        self.rxn2bounds = rxn2bounds
        self.i2bounds = {self.rxn2i[rxn]: bounds for rxn, bounds in self.rxn2bounds.iteritems()}
        self.rxnnum = len(self.rxn2i)
        self.cpdnum = len(self.cpd2i)
        self.obj = obj
        self.is_reversible = len(rxn2dir) == 0
        self.rxn2dir = rxn2dir
        self.rxn2fwd = {rxn: r[dir_fwd] for rxn, r in rxn2dir.iteritems()}
        self.rxn2rev = {rxn: r[dir_rev] for rxn, r in rxn2dir.iteritems() if "rev" in r}
        self.fwd2rxn = {fwd: rxn for rxn, fwd in self.rxn2fwd.iteritems()}
        self.rev2rxn = {rev: rxn for rxn, rev in self.rxn2rev.iteritems()}
        self.dir2rxn = {fwd: rxn for rxn, fwd in self.rxn2fwd.iteritems()}
        self.dir2rxn.update({rev: rxn for rxn, rev in self.rxn2rev.iteritems()})
        self.rxnnum_external = len(self.rxn2i_external)
        self.rxnnum_internal = len(self.rxn2i_internal)
        self.coupling = {'allow': {}, 'depend': {}}

    def set_coupling(self, coupling_table):
        for c in coupling_table:
            if c['REACTION1'] not in self.coupling['depend']: self.coupling['depend'][c['REACTION1']] = {dir_fwd: set(), dir_rev: set()}
            if c['REACTION2'] not in self.coupling['allow']: self.coupling['allow'][c['REACTION2']] = {dir_fwd: set(), dir_rev: set()}

            self.coupling['depend'][c['REACTION1']][c['DIRECTION1']].add(c['REACTION2'])
            self.coupling['allow'][c['REACTION2']][c['DIRECTION2']].add(c['REACTION1'])

    def coupling_dependencies(self, rxn, direction=dir_fwd):
        """
        list of reactions which RXN depends upon
        """
        return self.coupling['depend'].get(rxn, {dir_fwd: set(), dir_rev: set()})[direction]

    def coupling_dependants(self, rxn, direction=dir_fwd):
        """
        list of reactions which depend on reaction RXN
        """
        return self.coupling['allow'].get(rxn, {dir_fwd: set(), dir_rev: set()})[direction]

    def set_media(self, media):
        media_dict = {r['COMPOUND']: FbaBounds(r['LB'], r['UB']) for r in media}
        cplex_lb = []
        cplex_ub = []
        default_bounds = FbaBounds(-cplex.infinity, 0)

        for rxn, i in self.rxn2i_external.iteritems():
            if rxn == self.obj:
                continue

            cpd = self.rxn2cpd[rxn]
            lb = media_dict.get(cpd, default_bounds).lb
            ub = media_dict.get(cpd, default_bounds).ub
            cplex_lb.append((i, -cplex.infinity if lb <= 100 else lb))
            cplex_ub.append((i, cplex.infinity if ub >= 100 else ub))

        self.model.variables.set_lower_bounds(cplex_lb)
        self.model.variables.set_upper_bounds(cplex_ub)

    def switch_reactions(self, indices, values=None):
        """
        Set problem variable bounds. If bounds are not specified reactions are blocked
        """
        indices_range = range(len(indices))

        raw_lb = self.model.variables.get_lower_bounds(indices_range)
        raw_ub = self.model.variables.get_upper_bounds(indices_range)

        inf = float("Inf")

        backup_lb, backup_ub = {}, {}
        new_lb, new_ub = [], []
        for i in indices_range:
            ub_bin = int(raw_ub[i] != 0)
            ub = values[i] if values else self.i2bounds[i].ub
            ub = ub if ub < inf else cplex.infinity
            ub_changed = ub_bin != indices[i] or (values and raw_ub[i] != ub)
            if ub_changed:
                backup_ub[i] = raw_ub[i]
                new_ub.append((i, (ub if indices[i] else 0.0), ))

            lb_bin = int(raw_lb[i] != 0)
            lb = -values[i] if values and self.i2bounds[i].lb < 0 and -values[i] else self.i2bounds[i].lb
            lb = lb if lb > -inf else -cplex.infinity
            if lb > ub: lb = self.i2bounds[i].lb
            lb_changed = lb_bin != indices[i] or (values and raw_lb[i] != lb)
            if lb_changed:
                backup_lb[i] = raw_lb[i]
                new_lb.append((i, (lb if indices[i] else 0.0), ))

        if len(new_lb):
            self.model.variables.set_lower_bounds(new_lb)
        if len(new_ub):
            self.model.variables.set_upper_bounds(new_ub)

        return backup_lb, backup_ub

    def active_reactions2(self, tested_reactions_init):
        m = self.model
        reversible_reactions = set(r_i for r_i, b in self.i2bounds.iteritems() if b.lb < 0)

        active_reactions2 = {}
        m.objective.set_linear([(r_i, 0) for r_i in xrange(self.rxnnum)])
        for r_i in tested_reactions_init:
            m.objective.set_linear([(r_ii, 0) for r_ii in xrange(self.rxnnum)])
            m.objective.set_linear(r_i, 1)
            m.solve()
            status = m.solution.get_status_string()
            obj = m.solution.get_objective_value() if status == "optimal" else 0

            if status == "optimal" and abs(obj) > 1e-9:
                active_reactions2[r_i] = obj
                continue

                backup_lb = m.variables.get_lower_bounds(r_i)
                backup_ub = m.variables.get_upper_bounds(r_i)
                m.variables.set_lower_bounds(r_i, obj)
                m.variables.set_upper_bounds(r_i, obj)
                m.objective.set_linear([(r_ii, -1) for r_ii in xrange(self.rxnnum)])
                m.solve()
                active_reactions2[r_i] = sum([v > 0 for v in m.solution.get_values()])
                m.variables.set_lower_bounds(r_i, backup_lb)
                m.variables.set_upper_bounds(r_i, backup_ub)

            if r_i in reversible_reactions and obj == 0:
                m.objective.set_linear([(r_ii, 0) for r_ii in xrange(self.rxnnum)])
                m.objective.set_linear(r_i, -1)
                m.solve()
                status_rev = m.solution.get_status_string()
                obj_rev = m.solution.get_objective_value() if status_rev == "optimal" else 0

                if abs(obj_rev) > 1e-9:
                    active_reactions2[r_i] = obj_rev

                continue

                if status_rev == "optimal" and obj_rev != 0:
                    backup_lb = m.variables.get_lower_bounds(r_i)
                    backup_ub = m.variables.get_upper_bounds(r_i)
                    m.variables.set_lower_bounds(r_i, -obj_rev)
                    m.variables.set_upper_bounds(r_i, -obj_rev)
                    m.objective.set_linear([(r_ii, -1) for r_ii in xrange(self.rxnnum)])
                    m.solve()
                    active_reactions2[r_i] = sum([v > 0 for v in m.solution.get_values()])
                    m.variables.set_lower_bounds(r_i, backup_lb)
                    m.variables.set_upper_bounds(r_i, backup_ub)
            # if r_i in active_reactions2:
            #     print active_reactions2[r_i]

            m.objective.set_linear([(r_ii, 0) for r_ii in xrange(self.rxnnum)])

        return active_reactions2

    def active_reactions(self, individual):
        """
        This is not exactly the same as testing all reactions one by one, but it gives 99.7% correct answer on average.
        Doing a complete scan would be too slow.
        It is expected to run 15-20 simulations
        """
        tested_reactions_init = set(r_i for r_i, s in enumerate(individual) if s > 0)
        individual_vals = [s*(10000 if self.i2rxn[r_i] in self.rxn2i_external else 1) for r_i, s in enumerate(individual)]
        backup_lb, backup_ub = self.switch_reactions(individual, values=individual_vals)
        backup_obj = self.model.objective.get_linear()

        m = self.model

        # Find potentially active reactions
        active_reactions = {}
        tested_reactions = tested_reactions_init
        reversible_reactions = set(r_i for r_i, b in self.i2bounds.iteritems() if b.lb < 0)
        added_reactions = []

        dir = "fwd"
        force_dir = ""
        sample_n = 1
        chunked = 0
        for i in range(200):
            dir = "rev" if added_reactions and not added_reactions[-1] and dir == "fwd" else "fwd"
            if force_dir:
                dir = force_dir
                force_dir = ""

            sign = 1
            if dir == "rev":
                tested_reactions_this = reversible_reactions.intersection(tested_reactions)
                sign = -1
            else:
                tested_reactions_this = tested_reactions

            tested_reactions_this_list = list(tested_reactions_this)

            if sample_n > 1:
                random.shuffle(tested_reactions_this_list)
                chunked += 1

            iactive = {}
            for tested_reactions_chunk in list_utils.chunkify(tested_reactions_this_list, sample_n, jumpy=True):
                m.objective.set_linear([(r_i, int(sign*r_i in tested_reactions_chunk)) for r_i in xrange(self.rxnnum)])
                m.solve()

                fluxes = zip(tested_reactions_chunk, m.solution.get_values(tested_reactions_chunk))
                iactive.update({r_i: val for r_i, val in fluxes if abs(val) > 1e-9})

            # iactive_coupled = {}
            # for r_i in iactive:
            #     iactive_coupled.update({self.rxn2i[rxn2]: 10 for rxn2 in self.coupling_dependencies(self.i2rxn[r_i], dir)})
            # iactive.update(iactive_coupled)

            # iactive_dual = set(r_i for r_i, d in enumerate(m.solution.get_dual_values()) if d != 0)
            # iactive_dual = iactive_dual.intersection(tested_reactions)
            # iactive = iactive.union(iactive_dual)

            iactive = {r_i: val for r_i, val in iactive.iteritems() if r_i not in active_reactions}
            active_reactions.update(iactive)
            tested_reactions = tested_reactions_init - set(active_reactions)
            added_reactions.append(len(iactive))
            # print "{}[{}]: {} {}".format(i, sample_n, dir, len(active_reactions))

            if i > 5:
                if chunked == 0 and sample_n < 2 and sum(added_reactions[-2:]) == 0:
                    force_dir = "fwd" if dir == "rev" else "rev"
                    sample_n = 5
                elif sample_n > 1 :
                    sample_n = 1
                #     force_dir = "rev"
                # elif sample_n > 1 and chunked == 2:
                #     sample_n = 1
                elif sum(added_reactions[-2:]) == 0:
                    break

        # active_reactions2 = self.active_reactions2(tested_reactions_init)
        # print "VALIDATE {} (expected {})".format(len(active_reactions), len(active_reactions2))

        # print "CORRECT {}:".format(len(set(active_reactions2).union(set(active_reactions))))
        # for r_i in set(active_reactions2).union(set(active_reactions)):
        #     print "{}: {}[{}]".format(r_i, active_reactions.get(r_i, 0), active_reactions2.get(r_i, 0))

        # print "MISSING {}:".format(len(set(active_reactions2) - set(active_reactions)))
        # for r_i in set(active_reactions2) - set(active_reactions):
        #     print "{}: {} {}".format(self.i2rxn[r_i], active_reactions2.get(r_i, 0), "reversible" if r_i in reversible_reactions else "")

        # print "UNEXPECTED {}:".format(len(set(active_reactions) - set(active_reactions2)))
        # for r_i in set(active_reactions) - set(active_reactions2):
        #     print "{}: {}".format(r_i, active_reactions[r_i])

        # print "{}: {}".format(i, len(active_reactions))

        m.variables.set_lower_bounds(backup_lb.items())
        m.variables.set_upper_bounds(backup_ub.items())
        m.objective.set_linear([(r_i, val) for r_i, val in enumerate(backup_obj)])

        return active_reactions

def sbml2cplex(model, objective=None):
    obj_reaction = objective
    all_reactions, columns, lb, ub, obj = [], [], [], [], []

    all_compounds_dict = {m.id: m for m in model.species}
    all_compounds_set = set(m.id for m in model.species)
    all_compounds = list(all_compounds_set)
    all_compounds_ind = {cpd: i for i, cpd in enumerate(all_compounds)}
    all_reactions_ind, rxn2external, rxn2compound, rxn2bounds = {}, {}, {}, {}

    rxn2dir = {}
    r_i = 0
    for r in model.reactions:
        rxn2dir[r.id] = {dir_fwd: r.id}
        rxn2external[r.id] = False
        all_reactions.append(r.id)
        all_reactions_ind[r.id] = r_i

        cpds = []
        coefs = []
        for m in r.reactants:
            cpds.append(all_compounds_ind[m.species])
            coefs.append(-m.stoichiometry)
        for m in r.products:
            cpds.append(all_compounds_ind[m.species])
            coefs.append(m.stoichiometry)

        columns.append(cplex.SparsePair(cpds, coefs))
        r_bounds = FbaBounds(0 if r.reversible else -cplex.infinity, cplex.infinity)
        rxn2bounds[r.id] = r_bounds
        lb.append(r_bounds.lb)
        ub.append(r_bounds.ub)
        obj.append(float(r.id == obj_reaction))

        r_i += 1

    for cpd in all_compounds:
        if not all_compounds_dict[cpd].boundary_condition:
            continue

        r_id = "EX_{}".format(cpd)
        all_reactions.append(r_id)
        rxn2external[r_id] = True
        all_reactions_ind[r_id] = len(all_reactions_ind)
        rxn2bounds[r_id] = FbaBounds(-cplex.infinity, cplex.infinity)
        rxn2compound[r_id] = cpd
        columns.append(cplex.SparsePair([all_compounds_ind[cpd]], [1]))
        lb.append(-cplex.infinity)
        ub.append(cplex.infinity)
        obj.append(0)

    prob = cplex.Cplex()
    prob.objective.set_sense(prob.objective.sense.maximize)
    prob.linear_constraints.add(rhs=[0]*len(all_compounds), senses='E'*len(all_compounds), names=all_compounds)
    prob.variables.add(ub=ub, lb=lb, names=all_reactions, columns=columns, obj=obj)
    #prob.parameters.lpmethod.set(prob.parameters.lpmethod.values.auto)
    #prob.parameters.preprocessing.reduce.set(3)
    prob.parameters.advance.set(0)
    prob.set_log_stream(None)
    prob.set_error_stream(None)
    prob.set_warning_stream(None)
    prob.set_results_stream(None)

    return FbaProblem(model=prob, rxn2i=all_reactions_ind, cpd2i=all_compounds_ind, rxn2ext=rxn2external, rxn2cpd=rxn2compound,
               rxn2bounds=rxn2bounds, rxn2dir=rxn2dir, obj=obj_reaction)


def bioopt2cplex(bioopt, split_reversible=False, objective=None):
    obj_reaction = bioopt.objective.operands[0].name if objective is None and bioopt.objective is not None else objective
    all_reactions, columns, lb, ub, obj = [], [], [], [], []

    all_compounds_dict = {m.name: m for m in bioopt.find_metabolites()}
    all_compounds_set = set(m.name for m in bioopt.find_metabolites())
    all_compounds = list(all_compounds_set)
    all_compounds_ind = {cpd: i for i, cpd in enumerate(all_compounds)}
    all_reactions_ind, rxn2external, rxn2compound, rxn2bounds = {}, {}, {}, {}

    rxn2dir = {}
    r_i = 0
    for r in bioopt.reactions:
        if split_reversible and r.direction == dir_rev and obj_reaction != r.name:
            rxn2dir[r.name] = {}
            for dir in [dir_fwd, dir_rev]:
                r_id_dir = "{}_{}".format(r.name, dir)
                rxn2dir[r.name][dir] = r.name_dir
                rxn2external[r.name_dir] = False

                all_reactions.append(r_id_dir)
                all_reactions_ind[r_id_dir] = r_i

                cpds = []
                coefs = []
                for m in r.participants:
                    cpds.append(all_compounds_ind[m.metabolite.name])
                    coefs.append((1-int(dir == dir_rev)*2)*m.coefficient)

                r_lb = 0.0
                r_ub = r.bounds.ub if dir == dir_fwd else r.bounds.ub
                r_ub = r_ub if r_ub != r.bounds.inf() else cplex.infinity

                columns.append(cplex.SparsePair(cpds, coefs))
                lb.append(r_lb)
                ub.append(r_ub)
                rxn2bounds["{}_{}".format(r.name, dir)] = FbaBounds(r_lb, r_ub)
                obj.append(float(r.name == obj_reaction))

                r_i += 1
        else:
            rxn2dir[r.name] = {dir_fwd: r.name}
            rxn2external[r.name] = False
            all_reactions.append(r.name)
            all_reactions_ind[r.name] = r_i

            cpds = []
            coefs = []
            for m in r.participants:
                cpds.append(all_compounds_ind[m.metabolite.name])
                coefs.append(m.coefficient)

            columns.append(cplex.SparsePair(cpds, coefs))
            r_lb = r.bounds.lb if r.bounds.lb_is_finite else -cplex.infinity
            r_ub = r.bounds.ub if r.bounds.ub_is_finite else cplex.infinity
            lb.append(r_lb)
            ub.append(r_ub)
            rxn2bounds[r.name] = FbaBounds(r.bounds.lb, r.bounds.ub)
            obj.append(float(r.name == obj_reaction))

            r_i += 1

    for cpd in all_compounds:
        if not all_compounds_dict[cpd].boundary:
            continue

        r_id = "EX_{}".format(cpd)
        if split_reversible:
            rxn2dir[r_id] = {}
            for dir in [model.Direction.forward(), dir_rev]:
                r_id_dir = "{}_{}".format(r_id, dir)
                rxn2dir[r_id][dir] = r_id_dir
                all_reactions.append(r_id_dir)
                rxn2external[r_id_dir] = True
                all_reactions_ind[r_id_dir] = len(all_reactions_ind)
                rxn2bounds["{}_{}".format(r_id, dir)] = model.Bounds(0, cplex.infinity)
                rxn2compound[r_id_dir] = cpd
                columns.append(cplex.SparsePair([all_compounds_ind[cpd]], [1]))
                lb.append(0)
                ub.append(cplex.infinity)
                obj.append(0)
        else:
            r_id = "EX_{}".format(cpd)
            all_reactions.append(r_id)
            rxn2external[r_id] = True
            all_reactions_ind[r_id] = len(all_reactions_ind)
            rxn2bounds[r_id] = model.Bounds(-cplex.infinity, cplex.infinity)
            rxn2compound[r_id] = cpd
            columns.append(cplex.SparsePair([all_compounds_ind[cpd]], [1]))
            lb.append(-cplex.infinity)
            ub.append(cplex.infinity)
            obj.append(0)


    prob = cplex.Cplex()
    prob.objective.set_sense(prob.objective.sense.maximize)
    prob.linear_constraints.add(rhs=[0]*len(all_compounds), senses='E'*len(all_compounds), names=all_compounds)
    prob.variables.add(ub=ub, lb=lb, names=all_reactions, columns=columns, obj=obj)
    #prob.parameters.lpmethod.set(prob.parameters.lpmethod.values.auto)
    #prob.parameters.preprocessing.reduce.set(3)
    prob.parameters.advance.set(0)
    prob.set_log_stream(None)
    prob.set_error_stream(None)
    prob.set_warning_stream(None)
    prob.set_results_stream(None)

    return FbaProblem(model=prob, rxn2i=all_reactions_ind, cpd2i=all_compounds_ind, rxn2ext=rxn2external, rxn2cpd=rxn2compound,
               rxn2bounds=rxn2bounds, rxn2dir=rxn2dir, obj=obj_reaction)


def summary_dual(model, map={}):
    str = ""
    res = zip([map.get(n, n) for n in model.linear_constraints.get_names()], model.solution.get_dual_values())
    res = sorted([(n, dual) for n, dual in res if dual != 0], key=lambda x: x[0])
    res_str = "{{:<{}}} --> {{}}\n".format(max(len(n) for n, dual in res))

    for name, dual in res:
        str += res_str.format(name, dual)

    return str


def summary_primal(model, map={}):
    str = ""
    for rxn, val in zip(model.variables.get_names(), model.solution.get_values()):
        if val != 0:
            str += "{:<20}{}\n".format(map.get(rxn, rxn), val)

    return str


def summary(model, primal=False, dual=False, map={}):
    status = model.solution.get_status_string()

    if status == "optimal":
        obj = model.solution.get_objective_value()
        ret = "{} ({})".format(status, obj)
        if primal:
            ret += "\nPrimal\n=================\n"
            ret += summary_primal(model, map=map)
        if dual:
            ret += "\nDual\n=================\n"
            ret += summary_dual(model, map=map)
    else:
        ret = status

    return ret


def reaction_precursors(prob, reaction, hide_inf=True, map={}):
    obj_bck = list(enumerate(prob.model.objective.get_linear()))
    obj = [(i, 0.0) for i in xrange(prob.rxnnum)]
    prob.model.objective.set_linear(obj)

    rcol = prob.model.variables.get_cols(prob.rxn2i[reaction])
    res = []
    for cpd_i, coef in zip(rcol.ind, rcol.val):
        if coef == 0:
            continue

        cpd = prob.i2cpd[cpd_i]

        prob.model.variables.add(lb=[0], ub=[1000], names=["MET_EXPORT_TEST"], columns=[cplex.SparsePair([cpd_i], [-1])], obj=[1])
        prob.model.solve()

        obj = prob.model.solution.get_objective_value()
        if not is_optimal(prob.model) or (is_optimal(prob.model) and (not hide_inf or obj < 100)):
            res.append((cpd, coef, summary(prob.model)))
        prob.model.variables.delete(prob.rxnnum)

    prob.model.objective.set_linear(obj_bck)

    res = [(cpd + ":" + map.get(cpd, cpd), coef, sol) for cpd, coef, sol in res]
    cmd_max = max([len(c[0]) for c in res])
    res_str = "{{:>{}}} {{:>8}}: {{}}".format(cmd_max)

    res = sorted(res, key=lambda x: x[0])
    for cpd, coef, sol in res:
        print res_str.format(cpd, coef, sol)

