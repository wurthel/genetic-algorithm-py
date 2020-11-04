import configparser
from functools import partial

from constraints import Constraints, constraint_included, constraint_distances, constraint_max_charge, constraint_max_num_changes
from evolution import *
from logger import FileLogger
from utils import *

# PARSING CONFIG
config = configparser.ConfigParser()
config.read('config.ini')

pdb_file = config['PDB']['File']
value = float(config['PDB']['Value'])
cros_prob = float(config['PARAMS']['CrosProb'])
mut_prob = float(config['PARAMS']['MutProb'])
eval_param = float(config['PARAMS']['EvalParam'])
pop_count = int(config['PARAMS']['PopCount'])
pop_size = int(config['PARAMS']['PopSize'])
stop_step = int(config['PARAMS']['StopStep'])
compute_lmb_dir = config['COMPUTING']['ComputeLambdaDir']
computed_proteins_file = config['COMPUTING']['ComputedProteinsFileName']
result_file_name = config['COMPUTING']['ResultFileName']
population_from_computed = config['COMPUTING']['PopulationFromComputed']

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
from_computed = None
if population_from_computed.lower() == 'true':
    from_computed = True
elif population_from_computed.lower() == 'false':
    from_computed = False
computed_protein_saver = ProteinEvolutionSaver(computed_proteins_file)
evolution = []
for i in range(pop_count):
    working_dir = os.path.join(compute_lmb_dir, f'{i}')
    e = ProteinEvolution(population=None, mut_prob=mut_prob, cros_prob=cros_prob,
                         working_dir=working_dir, logger=FileLogger, save_function=computed_protein_saver, checker=constraints)
    e.generate_population(default_sequence=sequence, default_value=value, pop_size=pop_size, from_computed=from_computed)
    evolution.append(e)


async def main():
    async def evolution_step(e):
        e.mutation(attempts=4000)
        e.crossover(attempts=4000)
        await e.compute()
        e.selection(eval_param=0.05, save_n_best=3)

    logger = FileLogger('logout')
    iteration, step = 1, 0
    the_best_value = 0.0
    while step < stop_step:
        logger(f"Iteration: {iteration}\n")

        await asyncio.gather(*(evolution_step(e) for e in evolution))
        for e in evolution:
            e.print_info(iter=iteration)

        cur_best_value = max([e.get_best_protein().value for e in evolution])
        if the_best_value < cur_best_value:
            the_best_value = cur_best_value
            step = 0
        else:
            step += 1

        logger(f"The best value: {the_best_value}\n"
               f"Step/Stop {step}/{stop_step}\n\n")

        iteration += 1


asyncio.run(main())
