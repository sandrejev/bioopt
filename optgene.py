import argparse
from bioopt_parser import BiooptParser
from converter import Bioopt2CobraPyConverter
import cobra.io, cobra.flux_analysis
from deap import base, creator, tools
import random

class OptGene(object):
    def __init__(self, objective_reaction, objective_function='Yield', max_mutation=3,
                 target_type='Reaction', generations=10000, population_size=100, mutation_rate=1.0/500,
                 cx_fraction=0.8, cxMethod='Two', population=None, flux_calculation='FBA'):
        self.objective_reaction = objective_reaction
        self.objective_function = objective_function
        self.max_mutation = max_mutation
        self.target_type = target_type
        self.generations = generations
        self.population_size = population_size
        self.mutation_rate = mutation_rate
        self.cx_fraction = cx_fraction # crossover fracton
        self.cxMethod = cxMethod # crossover method, 'One', 'Two' or 'Uni'
        self.population = population # when start from previous result
        self.flux_calculation = flux_calculation # ToDo: implement MOMA

    def optimize(self, bioopt_model):
        global model, target
        model = self.reduceModel(bioopt_model)
        target = self.preprocessing(model)

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

        if self.cxMethod == 'Two':
            toolbox.register("mate", tools.cxTwoPoint)
        elif self.cxMethod == 'One':
            toolbox.register("mate", tools.cxOnePoint)
        elif self.cxMethod == 'Uni':
            toolbox.register("mate", tools.cxUniform, indpb = 0.5)
        else:
            print 'cxMethod of your choice is not registered'
        toolbox.register("mutate", tools.mutFlipBit, indpb=self.mutation_rate)
        toolbox.register("repair", self.__repair)
        toolbox.register("select", tools.selRoulette)

        # ToDo: Convert from bioopt to cobrapy

        # random.seed()
        if self.population == None:
            pop = toolbox.population(n=self.population_size)
        else:
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
        # print len(Evaluated)
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
                toolbox.repair(mutant)
                del mutant.fitness.values


            # Evaluate the individuals with an invalid fitness
            invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
            for ind in invalid_ind:
                if ind in Evaluated:
                    ind.fitness.values = Evaluated[Evaluated.index(ind)].fitness.values
                else:
                    ind.fitness.values = toolbox.evaluate(ind)
                    Evaluated.append(ind)
            # print len(Evaluated)

            print("  Evaluated %i individuals" % len(invalid_ind))

            # The population is entirely replaced by the offspring
            pop[:] = offspring

            hof.update(pop)
            currentBest = tools.selBest(pop, 1)[0]
            Rec.append(currentBest.fitness.values[0])

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

        ## show the progress until the end of the program
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
        from bioopt.converter import Bioopt2SbmlConverter
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

        # find reactions that can't carry a flux in current conditions
        flux_span_dict = cobra.flux_analysis.flux_variability_analysis(model, fraction_of_optimum=0.)
        blocked_reactions = [k for k, v in flux_span_dict.items() if max(map(abs, v.values())) < zero_cutoff]
        model.remove_reactions(blocked_reactions)
        print 'model reduced successfully'
        return model


    def preprocessing(self, model):
        global target
        if self.target_type == 'Reaction':
            target = model.reactions.list_attr('id')

            # remove exchange reactions
            target = [r for r in target if not ('xtO' in r or 'xtI' in r)]

            # remove lethal reactions from the target
            sol = model.optimize()
            tolerant_growth = sol.f * 0.01
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
            print "target_type should be 'Gene' or 'Teaction'"

        # caluculate the flux distribution by FBA
        sol = delModel.optimize()
        if sol.status == "infeasible" or sol.x_dict[self.objective_reaction] <= 1e-8 or sol.f < 1e-8:
            return 1e-8, # fitness shouldn't be zero to use selRoulette for the selection

        # return objective values
        if self.objective_function == "Yield":
            return sol.x_dict[self.objective_reaction],
        elif self.objective_function == "BPCY":
            return sol.x_dict[self.objective_reaction] * sol.f,
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
    parser.add_argument('--objective_reaction', '-objr', dest="objective_reaction", action='store', help='Reaction ID of reaction to optimize', type=str)
    parser.add_argument('--objective_function', '-objf', dest="objective_function", action='store', help='"Yield" or "BPCY"', type=str)
    parser.add_argument('--max_mutation', '-m', dest="max_mutation", default=3, action='store', help='Maximum mutation number', type=int)
    parser.add_argument('--target_type', '-t', dest="target_type", default='Reaction', action='store', help='"Reaction" or "Gene"', type=str)
    parser.add_argument('--generations', '-g', dest="generations", default=10000, action='store', help='Number of generations', type=int)
    parser.add_argument('--population-size', '-p', dest="population_size", default=125, action='store', help='Population size', type=int)
    parser.add_argument('--mutation_rate', '-mRate', dest="mutation_rate", default=0.002, action='store', help='Mutation rate of each element', type=float)
    parser.add_argument('--cx_fraction', '-cxf', dest="cx_fraction", default=0.8, action='store', help='Fraction of children generated by crossover', type=float)
    parser.add_argument('--cxMethod', '-cxm', dest="cxMethod ", action='store', default='Two', help='"Crossover method ("One", "Two" or "Uni")', type=str)
    parser.add_argument('--flux_calculation', '-fc', dest="flux_calculation", default='FBA', action='store', help='"Flux calculation method ("FBA" or "MOMA")', type=str)
    parser.add_argument('--sbml', dest="is_sbml", action='store_true', help='Is SBML file')

    # optional
    # parser.add_argument('--population', '-pop', dest="population", action='store', help='Fraction of children generated by crossover', type=toolbox.population)

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
                      max_mutation=args.max_mutation, target_type=args.target_type, generations=args.generations,
                      population_size=args.population_size, mutation_rate=args.mutation_rate, cx_fraction=args.cx_fraction)
                      # cxMethod=args.cxMethod)
    optgene.optimize(bioopt_model)