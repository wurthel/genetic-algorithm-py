import configparser
from evolution import *

# PARSING CONFIG
config = configparser.ConfigParser()
config.read('config.ini')

pattern_sequence = config['TEMPLATE']['Sequence']
cros_prob = float(config['PARAMS']['CrosProb'])
mut_prob = float(config['PARAMS']['MutProb'])
eval_param = float(config['PARAMS']['EvalParam'])
pop_size = int(config['PARAMS']['PopSize'])
compute_lmb_inf = config['COMPUTING']['ComputeLambdaInf']
compute_lmb_ouf = config['COMPUTING']['ComputeLambdaOuf']
computed_proteins_path = config['COMPUTING']['ComputedProteinsFileName']
result_file_name = config['COMPUTING']['ResultFileName']

# COMPUTING
computed_proteins = dict()
population = generate_population(sequence=pattern_sequence, pop_size=pop_size,
                                 compute_lmb_ouf=compute_lmb_ouf, compute_lmb_inf=compute_lmb_inf)
# save_computing(population, computed_proteins, computed_proteins_path)
# best_variance = None
# if os.path.exists(computed_proteins_path):
#     computed_proteins = read_computed_proteins(computed_proteins_path)
#     if computed_proteins:
#         best_variance = max(computed_proteins, key=computed_proteins.get)

with open("cavity2.txt", "r") as file:
    for line in file.readlines():
        pos = int(line)
        for protein in population:
            protein[pos - 1].cros_prob = 0.1
            protein[pos - 1].mut_prob = 0.2

iteration, step, stop_step = 1, 0, 5
while step < stop_step:
    best_proteins, population = selection(population, eval_param)
    crossover(population, cros_prob)
    mutation(population, mut_prob)
    compute_lambda(population, pattern_sequence, computed_proteins, compute_lmb_ouf, compute_lmb_inf)
    
    population = sorted(population, key=lambda p: p.value, reverse=True)
    population = best_proteins + population[0:pop_size-len(best_proteins)]

    cur_best_protein = get_best_protein(population)
    if best_variance is None:
        best_variance = ''.join(cur_best_protein.variance.values())
    elif cur_best_protein.value <= computed_proteins[best_variance]:
        step += 1
    else:
        best_variance = ''.join(cur_best_protein.variance.values())
        step = 0
    save_computing(population, computed_proteins, computed_proteins_path)
    line =  f'Iter: {iteration}. The best result: {best_variance}, {computed_proteins[best_variance]}'
    line += f' | {step}/{stop_step}'
    print(line)
    iteration += 1

# WRITING RESULTS
results = sorted(computed_proteins.items(), key=lambda x: x[1], reverse=True)
protein = Protein(pattern_sequence, genes)
with open(result_file_name, 'w') as ouf:
    for variance, value in results:
        protein.set_variance(variance)
        line = f'{value}, {variance}\n{protein.get_sequence(pattern_sequence)}\n\n'
        ouf.write(line)