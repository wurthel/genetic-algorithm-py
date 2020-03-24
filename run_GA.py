import configparser
from evolution import *
from constraints import Constraints, constraint_included, constraint_distances, constaint_charge
from functools import partial

# PARSING CONFIG
config = configparser.ConfigParser()
config.read('config.ini')

pdb_file = "6GUX_t.pdb"
cros_prob = float(config['PARAMS']['CrosProb'])
mut_prob = float(config['PARAMS']['MutProb'])
eval_param = float(config['PARAMS']['EvalParam'])
pop_size = int(config['PARAMS']['PopSize'])
compute_lmb_inf = config['COMPUTING']['ComputeLambdaInf']
compute_lmb_ouf = config['COMPUTING']['ComputeLambdaOuf']
computed_proteins_path = config['COMPUTING']['ComputedProteinsFileName']
result_file_name = config['COMPUTING']['ResultFileName']

# GENERATING CONSTRAINTS
constraints = Constraints()

coordinates = read_coordinates(pdb_file)

f1 = partial(constraint_included, aminoacids_set="DE", positions_set=PositionsSet1)
f2 = partial(constraint_distances, min_distance=5.0, coords=coordinates, positions_set=PositionsSetUnion)
f3 = partial(constaint_charge, max_charge=7)

constraints.add(f1)
constraints.add(f2)
constraints.add(f3)

# COMPUTING
computed_proteins = dict()
init_population = generate_population(pdb_file=pdb_file, pop_size=pop_size, computed_path=computed_proteins_path)
population = ProteinEvolution(init_population, mut_prob=mut_prob, cros_prob=cros_prob,
                              input_file=compute_lmb_inf, output_file=compute_lmb_ouf, save_file=computed_proteins_path,
                              checker=constraints)

iteration, step, stop_step = 1, 0, 5

the_best_value = 0
while step < stop_step:
    population.mutation(attempts_por_protein=500)
    population.crossover(attempts=500)
    population.compute()

    population.selection(eval_param=0.05, save_n_best=3)

    cur_best_value = population.get_best_protein().value
    if the_best_value < cur_best_value:
        the_best_value = cur_best_value
        step = 0
    else:
        step += 1

    population.save_to_file()

    print(f"Iteration: {iteration}\n"
          f"The best result: {the_best_value}\n"
          f"Step/Stop {step}/{stop_step}\n")

    iteration += 1

# WRITING RESULTS
