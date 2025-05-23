from collections import defaultdict

def nested_dict():
    return defaultdict(nested_dict)

def set_or_update_nested(c_dict, keys, value):
    for key in keys[:-1]:
        c_dict = c_dict[key]
    if keys[-1] in c_dict:
        c_dict[keys[-1]] += value
    else:
        c_dict[keys[-1]] = value

count_dict = nested_dict()

with open("work_folder/inferred_search_space/aggregated/methods/ms_y2h_experimental_wise.csv", "r") as f:
    header=True
    for line in f:
        if header:
            header=False
        else:
            _, _, n_test, n_obs = line.strip().split("\t")

            set_or_update_nested(count_dict, [n_test, n_obs], 1)


with open("negatome_plot/n_tested_obs.csv", "w") as w:
    for n_test, obs_dict in count_dict.items():
        for n_obs, n_occ in obs_dict.items():
            w.write(
                "\t".join(
                    map(str,
                        [
                            n_test,
                            n_obs,
                            n_occ
                        ]
                        )
                ) + "\n"
            )