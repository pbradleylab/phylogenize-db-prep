def get_targets(wildcards):
    """Get the path to the target database"""
    target = config["target_db"]["paths"].get(wildcards.target_db)
    if target:
        return target
    else:
        raise ValueError(f"Target not found for database {wildcards.target_db}")


# Create databases from targets
rule make_targets:
    input: get_targets
    output:
        index="results/{database}/make_targets/{target_db}/{target_db}.index",
        target_path=directory("results/{database}/make_targets/{target_db}/")
    conda: "../envs/target.yml"
    log: "logs/{database}/make_targets/{target_db}.log"
    threads: config["mmseqs2"]["createdb"]["threads"]
    shell:
        """
        mmseqs createdb {input} {output.target_path}/{wildcards.target_db} --dbtype 1 2> {log}
        mmseqs createindex {output.target_path}/{wildcards.target_db} /tmp 2> {log}
        """
