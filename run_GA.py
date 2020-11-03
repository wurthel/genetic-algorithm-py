import configparser
from functools import partial
from constraints import Constraints, constraint_included, constraint_distances, constraint_max_charge, constraint_max_num_changes
from evolution import *
from logger import FileLogger
from utils import *
import threading

# PARSING CONFIG
config = configparser.ConfigParser()
config.read('config.ini')

pdb_file = config['PDB']['File']
value = float(config['PDB']['Value'])
cros_prob = float(config['PARAMS']['CrosProb'])
mut_prob = float(config['PARAMS']['MutProb'])
mut_num = int(config['PARAMS']['MutNum'])
eval_param = float(config['PARAMS']['EvalParam'])
pop_num = int(config['PARAMS']['PopNum'])
pop_size = int(config['PARAMS']['PopSize'])
compute_lmb_dir = config['COMPUTING']['ComputeLambdaDir']
computed_proteins_file = config['COMPUTING']['ComputedProteinsFileName']
result_file_name = config['COMPUTING']['ResultFileName']

# GENERATING CONSTRAINTS
constraints = Constraints()

coordinates = read_coordinates(pdb_file)
sequence = read_sequence(pdb_file)

f1 = partial(constraint_included, aminoacids_set="DE", positions_set=PositionsSet1)
f2 = partial(constraint_distances, min_distance=5.0, coords=coordinates, positions_set=PositionsSetUnion)
f3 = partial(constraint_max_charge, max_charge=7)
f4 = partial(constraint_max_num_changes, max_num_changes=10)

constraints.add(f1)
constraints.add(f2)
constraints.add(f3)
constraints.add(f4)

# COMPUTING
evolution = []
computed_protein_saver = ProteinEvolutionSaver(computed_proteins_file)
for i in range(pop_num):
    working_dir = os.path.join(compute_lmb_dir, f'{i}')
    cur_evolution = ProteinEvolution(population=None, mut_prob=mut_prob, mut_num=mut_num, cros_prob=cros_prob,
                                     working_dir=working_dir, logger=FileLogger, save_function=computed_protein_saver, checker=constraints)
    cur_evolution.generate_population(default_sequence=sequence, default_value=value, pop_size=pop_size, from_computed=True)
    evolution.append(cur_evolution)


def evolution_step(e):
    e.mutation(attempts=4000)
    e.crossover(attempts=4000)
    e.compute()
    e.selection(eval_param=0.05, save_n_best=3)


logger = FileLogger('logout')
iteration, step, stop_step = 1, 0, 5
the_best_value = 0
while step < stop_step:
    logger(f"Iteration: {iteration}\n")

    tasks = [threading.Thread(target=evolution_step, args=(e,)) for e in evolution]
    for task in tasks:
        task.start()

    for task, e in zip(tasks, evolution):
        task.join()
        e.print_info(iteration)

    cur_best_value = max([e.get_best_protein().value for e in evolution])
    if the_best_value < cur_best_value:
        the_best_value = cur_best_value
        step = 0
    else:
        step += 1

    logger(f"The best value: {the_best_value}\n"
           f"Step/Stop {step}/{stop_step}\n\n")

    iteration += 1
