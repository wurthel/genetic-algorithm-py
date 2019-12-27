import configparser
from evolution import *

# PARSING CONFIG
config = configparser.ConfigParser()
config.read('config.ini')

genes = read_genes(config['GENES'])
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

population = generate_population(pop_size, genes, compute_lmb_ouf, compute_lmb_inf)
save_computing(population, computed_proteins, computed_proteins_path)

best_variance = None
if os.path.exists(computed_proteins_path):
    computed_proteins = read_computed_proteins(computed_proteins_path)
    if computed_proteins:
        best_variance = max(computed_proteins, key=computed_proteins.get)

iteration, step, stop_step = 1, 0, 5
while step < stop_step:
    population = selection(population, eval_param)
    crossover(population, genes, cros_prob)
    mutation(population, genes, mut_prob)
    compute_lambda(population, pattern_sequence, computed_proteins, compute_lmb_ouf, compute_lmb_inf)
    cur_best_protein = get_best_protein(population)
    if best_variance == None:
        best_variance = cur_best_protein.get_variance()
    elif cur_best_protein.mlambda <= computed_proteins[best_variance]:
        step += 1
    else:
        best_variance = cur_best_protein.get_variance()
        step = 0
    save_computing(population, computed_proteins, computed_proteins_path)
    line =  f'Iter: {iteration}. The best result: {best_variance}, {computed_proteins[best_variance]}'
    line += f' | {step}/{stop_step}'
    print(line)
    iteration += 1

# WRITING RESULTS
results = sorted(computed_proteins.items(), key=lambda x: x[1], reverse=True)
protein = Protein(genes)
with open(result_file_name, 'w') as ouf:
    for variance, lmb in results:
        protein.set_variance(variance)
        line = f'{lmb}, {variance}\n{protein.get_sequence(pattern_sequence)}\n\n'
        ouf.write(line)