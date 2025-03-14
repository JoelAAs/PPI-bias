
rule placeholder_POD_calc:
    output:
        prey_pod = "work_folder/cell_type_pod/{cell_line}_pod.csv"
    run:
        pob_file = config["cell_lines"][wildcards.cell_line]["pod"]
        if pob_file:
            shell(f"cp {pob_file} {output.prey_pod}")
        else:
            shell(f"touch {output.prey_pod}")