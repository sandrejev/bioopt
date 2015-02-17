import argparse
from bioopt_parser import BiooptParser
from converter import Bioopt2CobraPyConverter
import cobra.io, cobra.flux_analysis, cobra.manipulation
import moma
from deap import base, creator, tools
import random

class OptGene(object):
    def __init__(self, objective_reaction, objective_function, max_mutation,
                 target_type, target_list, generations, population_size, mutation_rate,
                 cx_fraction, cx_method, population, flux_calculation, wt_flux):
        self.objective_reaction = objective_reaction
        self.objective_function = objective_function
        self.max_mutation = max_mutation
        self.target_type = target_type
        self.target_list = target_list
        self.generations = generations
        self.population_size = population_size
        self.mutation_rate = mutation_rate
        self.cx_fraction = cx_fraction # crossover fraction
        self.cx_method = cx_method # crossover method, 'One', 'Two' or 'Uni'
        self.population = population # when start from previous result
        self.flux_calculation = flux_calculation
        self.wt_flux = wt_flux


    def optimize(self, bioopt_model):
        global model, target
        #ToDo: modify to allow users to provide reduced model and target list
        model = self.reduceModel(bioopt_model) # remove blocked reactions and reactions by isozymes
        target = self.preprocessing(model) # define deletion target list, exchange reactions and lethal reactions are removed from target
        if not self.wt_flux:
            self.wt_flux = moma.optimize_minimum_flux(model)

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
        toolbox.register("evaluate", self.__evaluation)

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
        if self.population == None: # start from WT
            pop = toolbox.population(n=self.population_size)
        else: # start from previous result
            pop = self.population

        Rec = [] # record the best fitness values in each generation here
        hof = tools.HallOfFame(5)  # the list of five best individuals found so far

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
            offspring = toolbox.select(pop, len(pop))
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
            Rec.append(currentBest.fitness.values[0]) # fitness value of the best individual in this generation is recorded

            # Gather all the fitnesses in one list and print the stats
            fits = [ind.fitness.values[0] for ind in pop]

            length = len(pop)
            mean = sum(fits) / length
            sum2 = sum(x*x for x in fits)
            std = abs(sum2 / length - mean**2)**0.5

            # print("  Min %s" % min(fits))
            print("  Max %s" % max(fits))
            print("  Avg %s" % mean)
            # print("  Std %s" % std)


        print("-- End of (successful) evolution --")
        print("%s different mutants were evaluated" % len(Evaluated))

        best_ind = hof[0]
        delList = [target[i] for i, v in enumerate(best_ind) if not v]  # list of deletion for best individual
        print("Best individual is %s, %s" % (delList, best_ind.fitness.values))

        # show the progress until the end of the program
        # with open('Rec_SUCCxtO_Yield_m4_gene_uni.pickle', 'rb') as infile:
        #     Rec = load(infile)
        import matplotlib.pyplot as plt
        plt.plot(Rec)
        plt.xlabel('Generation')
        plt.ylabel('Objective value')
        # plt.axis([0, Generations, 0, Hof[0].fitness.values[0]])
        plt.show()


    def reduceModel(self, bioopt_model, zero_cutoff=1e-9):
        global model
        model = Bioopt2CobraPyConverter().convert(bioopt_model)
        print 'bioopt model successfully converted to cobra model'
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
        print 'model reduced successfully'
        return model


    def preprocessing(self, model):
        global target
        if self.target_type == 'Reaction':
            target = model.reactions.list_attr('id')

            # remove exchange reactions from deletion target
            # ToDo: make it possible to detect exchange reaction in model with other format
            target = [r for r in target if not ('xtO' in r or 'xtI' in r)]

            # remove lethal reactions from deletion target
            sol = model.optimize()
            tolerant_growth = sol.f * 0.01 # reaction is considered lethal when it's remocved the growth rate decrases to less than 1% of WT
            growth_rates, statuses = cobra.flux_analysis.single_deletion(model, target,
                                                                         method="fba", element_type='reaction')
            for k in growth_rates.iterkeys():
                if growth_rates[k] < tolerant_growth:
                    target.remove(k)

        # ToDo further reduce the number of target (non gene-associated reactions, transport reactions, reactions in unassociated subsystem )

        # ToDo preprocessing of genes
        elif self.target_type == 'Gene':
            genesToRemove = [g.id for g in model.genes if not g.reactions]
            target = model.genes

        else:
            print "target_type should be 'Reaction' or 'Gene'"

        print "pre-process finished successfully"
        return target



    def __evaluation(self, individual):
        global model, target
        # get the index of deleted elements
        deletion_list = [target[i] for i, v in enumerate(individual) if not v]

        # deletion of genes (set upper and lower bound of reaction to zero according to gene-reaction rules)
        delModel = model.copy()
        if self.target_type == 'Gene':
            if not deletion_list == []: # avoid error in delete_model_genes
                cobra.manipulation.delete_model_genes(delModel, deletion_list, cumulative_deletions=False)
        elif self.target_type == 'Reaction':
            for r in deletion_list:
                delModel.reactions.get_by_id(r).knock_out()
        else:
            print "target_type should be 'Gene' or 'Reaction'"

        # caluculate the flux distribution by FBA
        if self.flux_calculation == 'FBA':
            sol = delModel.optimize()
            growth = sol.f
            status = sol.status
            objective_flux = sol.x_dict[self.objective_reaction]

        if self.flux_calculation == 'MOMA':
            # ToDO: import moma module
            sol_dict = moma.moma(model, delModel, minimize_norm=True, norm_flux_dict=self.wt_flux)
            growth = sol_dict['objective_value']
            status = sol_dict['status']
            objective_flux = abs(sol_dict['the_problem'].x_dict[self.objective_reaction])

        if status == "infeasible" or objective_flux <= 1e-8 or growth < 1e-8:
                return 1e-8, # fitness shouldn't be zero to use selRoulette for the selection

        # return objective values
        if self.objective_function == "Yield":
            return objective_flux,
        elif self.objective_function == "BPCY":
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
                        help='"Yield" or "BPCY"', type=str)
    parser.add_argument('--max_mutation', '-m', dest="max_mutation", default=3, action='store',
                        help='Maximum mutation number (default: 3)', type=int)
    parser.add_argument('--target_type', '-t', dest="target_type", default='Reaction', action='store',
                        help='"Reaction" or "Gene" (default: "Reaction")', type=str)
    parser.add_argument('--target_list', '-tl', dest="target_list", default=None, action='store',
                        help='List of target id for deletion', type=list)
    parser.add_argument('--generations', '-g', dest="generations", default=10000, action='store',
                        help='Number of generations (default: 10000)', type=int)
    parser.add_argument('--population-size', '-p', dest="population_size", default=125, action='store',
                        help='Population size (default: 125)', type=int)
    parser.add_argument('--mutation_rate', '-mRate', dest="mutation_rate", default=0.002, action='store',
                        help='Mutation rate of each element (default: 0.002)', type=float)
    parser.add_argument('--cx_fraction', '-cxf', dest="cx_fraction", default=0.8, action='store',
                        help='Fraction of children generated by crossover (default: 0.8)', type=float)
    parser.add_argument('--cx_method', '-cxm', dest="cx_method", default='Two', action='store',
                        help='Crossover method ("One", "Two"(default) or "Uni")', type=str)
    parser.add_argument('--population', '-pop', dest="population", default=None, action='store',
                        help='Population from which evolution start', type=list)
    parser.add_argument('--flux_calculation', '-fc', dest="flux_calculation", default='FBA', action='store',
                        help='"Flux calculation method ("FBA"(default) or "MOMA")', type=str)
    parser.add_argument('--wt_flux', '-wt', dest="wt_flux", default=None, action='store',
                        help='x_dict of wt-model used as the reference in MOMA', type=dict)
    parser.add_argument('--sbml', dest="is_sbml", action='store_true', help='Is SBML file')

    args = parser.parse_args()

    if args.is_sbml:
        # Convert to bioopt
        pass
    else:
        bioopt_model = BiooptParser().parse_file(args.bioopt)

    ## Use pickled model
    # from cPickle import load
    # with open(args.bioopt, 'rb') as f:
    #     model = load(f)

    optgene = OptGene(objective_reaction=args.objective_reaction, objective_function=args.objective_function,
                      max_mutation=args.max_mutation, target_type=args.target_type, target_list=args.target_list,
                      generations=args.generations, population_size=args.population_size, mutation_rate=args.mutation_rate,
                      cx_fraction=args.cx_fraction, cx_method=args.cx_method, flux_calculation=args.flux_calculation,
                      population=args.population, wt_flux=args.wt_flux)
    optgene.optimize(bioopt_model)