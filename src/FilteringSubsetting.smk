import pandas as pd
import re
### functions
def remove_pattern(string, pattern):
    match = re.search(pattern, string).group()
    if match:
        return match[0]
    else:
        return None

rule get_bait_prey_subset:
    input:
        intact_human="data/human/human.txt"
    output:
        bait_prey_filename=""
    run:
        intact_df = pd.read_csv(input.intact_human, sep = "\t")
