import argparse
from bioopt_parser import BiooptParser
from converter import Bioopt2CobraPyConverter
import cobra.io, cobra.flux_analysis, cobra.manipulation
import moma # when run on cluster
from deap import base, creator, tools
import random
from cPickle import load, dump

class OptGene(object):
    def __init__(self, objective_function='Yield', max_mutation=3, target_type='reaction',
                 target_list=None, generations=5000, population_size=125, mutation_rate=None, cx_fraction=0.8,
                 cx_method='Two', record_file=None, flux_calculation='FBA', wt_flux=None, reduced=False, record_interval=0):
        self.objective_function = objective_function
        self.max_mutation = max_mutation
        self.target_type = target_type
        self.target_list = target_list
        self.generations = generations
        self.population_size = population_size
        self.mutation_rate = mutation_rate
        self.cx_fraction = cx_fraction # crossover fraction
        self.cx_method = cx_method # crossover method, 'One', 'Two' or 'Uni'
        self.record_file = record_file # when start from previous result
        self.flux_calculation = flux_calculation
        self.wt_flux = wt_flux
        self.reduced = reduced
        self.record_interval = record_interval


    def optimize(self, cobra_model, objective_reaction):
        '''
        Main function which uses genetic algorithm
        :param cobra_model: cobra.Model object
        :param objective_reaction: reaction ID of design objective
        :return: best individual
        '''
        global model, target
        # model: reduced model
        if self.reduced:
            model = cobra_model
        else:
            model = self.reduceModel(cobra_model) # remove blocked reactions and reactions by isozymes

        if self.flux_calculation.lower() == 'moma':
            if self.wt_flux:
                if self.wt_flux.endswith('.pickle'):
                    with open(self.wt_flux, 'rb') as f:
                        self.wt_flux = load(f)
                elif self.wt_flux.endswith('.txt'):
                    with open(self.wt_flux) as f:
                        self.wt_flux = {}
                        for line in f:
                            rid, flux = line.split()
                            self.wt_flux[rid] = float(flux)
                else:
                    raise ValueError("wild type flux data should be .txt or .pickle file")

            else:
                # ToDo: figure out why cplex returns 'infeasible' when this is run on cluster
                self.wt_flux = moma.optimize_minimum_flux(model)

        # define target of deletion
        if self.target_list:
            if self.target_list.endswith('.pickle'):
                with open(self.target_list, 'rb') as f:
                    target = load(f)
            elif self.target_list.endswith('.txt'):
                with open(self.target_list) as f:
                    target = [line.strip() for line in f]
            else:
                raise ValueError("target list should be .txt or .pickle file")

        else:
            target = self.preprocessing(model) # define deletion target list, exchange reactions and lethal reactions are removed from target
        print 'number of target: %s' % len(target)

        if not self.mutation_rate:
            self.mutation_rate = 1.0/len(target) # see paper
        print 'mutation rate: %s' % self.mutation_rate

        # for details, see the documentation of DEAP (https://code.google.com/p/deap/)
        # create classes
        creator.create("FitnessMax", base.Fitness, weights=(1.0,))
        creator.create("Individual", list, fitness=creator.FitnessMax)
        toolbox = base.Toolbox()

        # Attribute generator
        toolbox.register("attr_bool", int, 1)
        # Structure initializers
        toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_bool, len(target))
        toolbox.register("population", tools.initRepeat, list, toolbox.individual)

        # Operator registering
        toolbox.register("evaluate", self.__evaluation, objective_reaction=objective_reaction)

        if self.cx_method == 'Two':
            toolbox.register("mate", tools.cxTwoPoint)
        elif self.cx_method == 'One':
            toolbox.register("mate", tools.cxOnePoint)
        elif self.cx_method == 'Uni':
            toolbox.register("mate", tools.cxUniform, indpb=0.5)
        else:
            print 'cx_method of your choice is not registered'
        toolbox.register("mutate", tools.mutFlipBit, indpb=self.mutation_rate)
        toolbox.register("repair", self.__repair)
        toolbox.register("select", tools.selRoulette)

        # random.seed() # set seed when needed
        if self.record_file is None: # start from WT
            pop = toolbox.population(n=self.population_size)
            hof = tools.HallOfFame(5)  # the list of five best individuals found so far
            progress = [] # record the best fitness values in each generation here

        else: # start from previous result
            with open(self.record_file) as f:
                Rec = load(f)
            pop = Rec['population']
            hof = Rec['hof']
            progress = Rec['progress']


        print("Start of evolution")
        # Evaluate the entire population
        Evaluated = [] # store already evaluated individuals
        for ind in pop:
            if ind in Evaluated:
                 ind.fitness.values = Evaluated[Evaluated.index(ind)].fitness.values
            else:
                ind.fitness.values = toolbox.evaluate(ind)
                Evaluated.append(ind)
        hof.update(pop)

        print("  Evaluated %i individuals" % len(pop))

        # Begin the evolution
        for g in range(self.generations):
            print("-- Generation %i --" % g)

            # Select the next generation individuals
            offspring = toolbox.select(pop, len(pop)-2)
            if len(hof) >= 2:
                offspring.extend(hof[:2]) # two best individuals are automatically included
            else:
                offspring.extend([hof[0], hof[0]])

            # Clone the selected individuals
            offspring = list(map(toolbox.clone, offspring))

            # Apply crossover and mutation on the offspring
            random.shuffle(offspring) # to crossover random combination
            for child1, child2 in zip(offspring[::2], offspring[1::2]):
                if random.random() < self.cx_fraction:
                    toolbox.mate(child1, child2)
                    del child1.fitness.values
                    del child2.fitness.values

            for mutant in offspring:
                toolbox.mutate(mutant)
                toolbox.repair(mutant) # repair too many mutations
                del mutant.fitness.values


            # Evaluate the individuals with an invalid fitness
            invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
            for ind in invalid_ind:
                if ind in Evaluated:
                    ind.fitness.values = Evaluated[Evaluated.index(ind)].fitness.values
                else:
                    ind.fitness.values = toolbox.evaluate(ind)
                    Evaluated.append(ind)

            print("  Evaluated %i individuals" % len(invalid_ind))

            # The population is entirely replaced by the offspring
            pop[:] = offspring

            hof.update(pop)
            currentBest = tools.selBest(pop, 1)[0]
            progress.append(currentBest.fitness.values[0]) # fitness value of the best individual in this generation is recorded

            # Gather all the fitnesses in one list and print the stats
            fits = [ind.fitness.values[0] for ind in pop]

            length = len(pop)
            mean = sum(fits) / length
            sum2 = sum(x*x for x in fits)
            std = abs(sum2 / length - mean**2)**0.5

            print("  Min %s" % min(fits))
            print("  Max %s" % max(fits))
            print("  Avg %s" % mean)
            print("  Std %s" % std)

            # record results every specified generation
            if self.record_interval:
                try:
                    if g != 0 and (g+1) % self.record_interval == 0:
                        with open('Rec_%s_%s_m%s_%s_%s.pickle' % (objective_reaction, self.objective_function,
                                                                      self.max_mutation, self.target_type,
                                                                      self.flux_calculation), 'wb') as o1:
                            dump({'population': pop, 'hof': hof, 'progress': progress}, o1)
                except:
                    print('improper record interval')

        print("-- End of (successful) evolution --")
        print("%s different mutants were evaluated" % len(Evaluated))

        best_ind = hof[0]
        delList = [target[i] for i, v in enumerate(best_ind) if not v]  # list of deletion for best individual
        print("Best mutation combination is %s, %s" % (delList, best_ind.fitness.values))

        # show the progress until the end of the program
        try:
            import matplotlib.pyplot as plt
            plt.plot(progress)
            plt.xlabel('Generation')
            plt.ylabel('Objective value')
            # plt.axis([0, Generations, 0, Hof[0].fitness.values[0]])
            # plt.show()
            file_name = 'progress_%s_%s_m%s_%s_%s.png' % (objective_reaction, self.objective_function,
                                                        self.max_mutation, self.target_type, self.flux_calculation)
            plt.savefig(file_name)
            print("Progress curve was saved as %s" % file_name)
        except:
            pass


    def reduceModel(self, model, zero_cutoff=1e-12):
        '''
        :param model: cobra.Model object
        :param zero_cutoff: reactions which can't have flux larger than this value are removed
        :return: reduced model
        '''
        Reactions = model.reactions

        # find same reaction with different names
        ReactionFormulaList = [r.reaction for r in Reactions]
        unique, result = [], [] # reaction formula of overlapping reactions
        for r in ReactionFormulaList:
            if r not in unique:
                unique.append(r)
            else:
                if r not in result:
                    result.append(r)

        # get IDs of overlapping reactions
        def reactionInList(x):
            return x in result

        overLappingReactions = Reactions.query(reactionInList, 'reaction')

        uniqueReactions, uniqueID, toRemove = [], [], []
        for r in overLappingReactions:
            if r.reaction not in uniqueReactions:
                uniqueReactions.append(r.reaction)
                uniqueID.append(r.id)
            else:
                toRemove.append(r.id)

        # remove unnecessarily reactions
        model.remove_reactions(toRemove)


        # find reactions that can't carry a flux in current conditions (blocked reactions)
        flux_span_dict = cobra.flux_analysis.flux_variability_analysis(model, fraction_of_optimum=0.)
        blocked_reactions = [k for k, v in flux_span_dict.items() if max(map(abs, v.values())) < zero_cutoff]
        model.remove_reactions(blocked_reactions)
        print 'model reduced successfully: %s reactions removed' % (len(toRemove) + len(blocked_reactions))
        return model


    def preprocessing(self, model, target_type=None, method=None, tolerant_growth=None):
        '''
        Remove exchange reactions and computationally lethal reactions / genes from the target
        :param model: cobra.Model object
        :param  tolerant_growth: To find lethal reactions/genes (default: 1% of maximum growth)
        :return: target list
        '''
        global target
        if not target_type:
            target_type = self.target_type.lower()

        if not method:
            method = self.flux_calculation.lower()

        if target_type.lower() == 'reaction':
            # ToDo: properly deal with reactions controlled by same genes
            # geneList = [r.genes for r in Reactions if r.genes]
            # uniqueGeneList = list(set(geneList))

            # list of reactions, exchange reactions are not included
            target = [r.id for r in model.reactions if not len(r.metabolites) == 1]

        elif target_type.lower() == 'gene':
            target = [g.id for g in model.genes if g.reactions]

        # remove lethal reactions from the target
        if not tolerant_growth:
            sol = model.optimize()
            tolerant_growth = sol.f * 0.01

        if method.lower() == 'fba':
            growth_rates, statuses = cobra.flux_analysis.single_deletion(model, target, method='fba',
                                                                         element_type=target_type)

        elif method.lower() == 'moma':
            growth_rates, statuses = {}, {}
            for t in target:
                mutant = model.copy()
                if target_type.lower() == 'gene':
                    cobra.manipulation.delete_model_genes(mutant, [t], cumulative_deletions=False)
                elif target_type.lower() == 'reaction':
                    mutant.reactions.get_by_id(t).knock_out()

                sol = moma.moma(model, mutant, norm_flux_dict=self.wt_flux)
                try:
                    growth_rates[t] = sol['objective_value']
                    statuses[t] = sol['status']
                except:
                    growth_rates[t] = None
                    statuses[t] = 'failed'

        for k in growth_rates.iterkeys():
            if growth_rates[k] < tolerant_growth:
                target.remove(k)

        # ToDo: further reduce the number of target (non gene-associated reactions, transport reactions, reactions in unassociated subsystem )

        print "pre-process finished successfully"
        return target


    def __evaluation(self, individual, objective_reaction):
        global model, target
        # get the index of deleted elements
        deletion_list = [target[i] for i, v in enumerate(individual) if not v]

        # deletion of genes (set upper and lower bound of reaction to zero according to gene-reaction rules)
        delModel = model.copy()
        if self.target_type.lower() == 'gene':
            if not deletion_list == []: # avoid error in delete_model_genes
                cobra.manipulation.delete_model_genes(delModel, deletion_list, cumulative_deletions=False)
        elif self.target_type.lower() == 'reaction':
            for r in deletion_list:
                delModel.reactions.get_by_id(r).knock_out()
        else:
            print "target_type should be 'Gene' or 'Reaction'"

        # caluculate the flux distribution by FBA
        if self.flux_calculation.lower() == 'fba':
            sol = delModel.optimize()
            growth = sol.f
            status = sol.status
            try:
                objective_flux = sol.x_dict[objective_reaction]
            except TypeError:
                return 1e-16,

        elif self.flux_calculation.lower() == 'moma':
            sol_dict = moma.moma(model, delModel, minimize_norm=True, norm_flux_dict=self.wt_flux)
            try:
                growth = sol_dict['objective_value']
                status = sol_dict['status']
                objective_flux = abs(sol_dict['the_problem'].x_dict[objective_reaction])
            except:
                return 1e-16,

        if status == "infeasible" or objective_flux <= 1e-16 or growth < 1e-16:
                return 1e-16, # fitness shouldn't be zero to use selRoulette for the selection

        # return objective values
        if self.objective_function.lower() == "yield":
            return objective_flux,
        elif self.objective_function.lower() == "bpcy":
            return objective_flux * growth,
        else:
            print "object_function should be 'Yield' or 'BPCY'"

    def __repair(self, individual):
        # get the index of mutated elements
        mutList = [i for i, v in enumerate(individual) if not v]
        # remove mutation randomly while the number of mutation is larger than maximum
        while len(mutList) > self.max_mutation:
            toRepair = random.choice(mutList)
            individual[toRepair] = type(individual[toRepair])(not individual[toRepair])
            mutList.remove(toRepair)
        return individual,


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run OptGene script')
    parser.add_argument('bioopt', action='store', help='Path to bioopt file')
    parser.add_argument('--objective_reaction', '-objr', dest="objective_reaction", action='store',
                        help='Reaction ID of reaction to optimize', type=str)
    parser.add_argument('--objective_function', '-objf', dest="objective_function", action='store',
                        help='"Yield" or "BPCY"(only valid when growth is objective)', type=str)
    parser.add_argument('--max_mutation', '-m', dest="max_mutation", default=3, action='store',
                        help='Maximum mutation number (default: 3)', type=int)
    parser.add_argument('--target_type', '-t', dest="target_type", default='Reaction', action='store',
                        help='"Reaction" or "Gene" (default: "Reaction")', type=str)
    parser.add_argument('--target_list', '-tl', dest="target_list", default=None, action='store',
                        help='file name of list of target id for deletion (.txt (one reaction ID in each line) or .pickle)', type=str)
    parser.add_argument('--generations', '-g', dest="generations", default=5000, action='store',
                        help='Number of generations (default: 5000)', type=int)
    parser.add_argument('--population-size', '-p', dest="population_size", default=125, action='store',
                        help='Population size (default: 125)', type=int)
    parser.add_argument('--mutation_rate', '-mRate', dest="mutation_rate", default=None, action='store',
                        help='Mutation rate of each element (default: 1/(length of target))', type=float)
    parser.add_argument('--cx_fraction', '-cxf', dest="cx_fraction", default=0.8, action='store',
                        help='Fraction of children generated by crossover (default: 0.8)', type=float)
    parser.add_argument('--cx_method', '-cxm', dest="cx_method", default='Two', action='store',
                        help='Crossover method ("One", "Two"(default) or "Uni")', type=str)
    parser.add_argument('--record_file', '-rec', dest="record_file", default=None, action='store',
                        help='File name of record file if run from past record', type=str)
    parser.add_argument('--flux_calculation', '-fc', dest="flux_calculation", default='FBA', action='store',
                        help='"Flux calculation method ("FBA"(default) or "MOMA")', type=str)
    parser.add_argument('--wt_flux', '-wt', dest="wt_flux", default=None, action='store',
                        help='filename of wt-model flux distribution used as the reference in MOMA', type=str)
    parser.add_argument('--sbml', dest="is_sbml", action='store_true', help='Is SBML file')
    parser.add_argument('--cobra', dest="is_cobra", action='store_true', help='Is COBRA model file')
    parser.add_argument('--reduced', dest="is_reduced", action='store_true', help='Model is already reduced')
    parser.add_argument('--record_interval', '-rec_step', dest="record_interval", default=0, action='store',
                        help='interval of recording results (default: no recording)', type=int)

    args = parser.parse_args()

    if args.is_sbml:
        # Convert to bioopt
        pass
    elif args.is_cobra:
        with open(args.bioopt, 'rb') as f:
            cobra_model = load(f)
    else:
        bioopt_model = BiooptParser().parse_file(args.bioopt)
        cobra_model = Bioopt2CobraPyConverter().convert(bioopt_model)
        print 'bioopt model successfully converted to cobra model'


    ## Use pickled model
    # from cPickle import load
    # with open(args.bioopt, 'rb') as f:
    #     model = load(f)

    optgene = OptGene(objective_function=args.objective_function,
                      max_mutation=args.max_mutation, target_type=args.target_type, target_list=args.target_list,
                      generations=args.generations, population_size=args.population_size, mutation_rate=args.mutation_rate,
                      cx_fraction=args.cx_fraction, cx_method=args.cx_method, flux_calculation=args.flux_calculation,
                      record_file=args.record_file, wt_flux=args.wt_flux, reduced=args.is_reduced, record_interval=args.record_interval)
    optgene.optimize(cobra_model, objective_reaction=args.objective_reaction)