import pandas as pd

target_name = '1_Amino_2_propanol'

pathway_sequence = ['pks','bio']
num_bio_steps = 2 # enter None if pathway sequence is ['pks'] only

input_filepath = None
reaction_rules_filepath = '../data/coreactants_and_rules/JN3604IMT_rules.tsv'
reaction_rules_df = pd.read_csv(reaction_rules_filepath,delimiter='\t')

if pathway_sequence == ['pks']:
    input_filepath = f'../data/results/{target_name}_PKS_only.txt'

if pathway_sequence == ['pks','bio']:
    input_filepath = f'../data/results_logs/{target_name}_PKS_BIO{num_bio_steps}.txt'

output_filepath = f"{input_filepath.rstrip('.txt')}_with_enzymes.txt"

with open(input_filepath, 'r') as file:
    lines = file.readlines()

rxn_rules_count = 0

for i,line in enumerate(lines):
    if 'reaction rules:' in line:
        rxn_rules_count += 1
        rxn_rules = line.lstrip('reaction rules: ')[1:-2].split(',')
        enzymes_dict = {}
        for rule in rxn_rules:
            rule_name = rule.strip()[1:-1]
            enzymes = list(reaction_rules_df[reaction_rules_df['Name']==rule_name]['Comments'])[0]
            enzymes_dict.update({f'{rule_name} enzyme UNIPROT IDs':enzymes})
        lines[i] = f'    {str(enzymes_dict)}\n'

if rxn_rules_count != 0:
    with open(output_filepath,'w') as new_file:
        new_file.writelines(lines)

else:
    print(f'\nNo pathways found. As such, no enzymes were queried.')


