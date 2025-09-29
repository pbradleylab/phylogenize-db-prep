def get_targets(wildcards):
    """Get the path to the target database"""
    target = config["target_db"]["paths"].get(wildcards.target_db)
    if target:
        return(target)
    else:
        raise ValueError(f"Target not found for database {wildcards.target_db}")


# Create databases from targets
rule make_targets:
    input: get_targets
    output:
        index="results/{database}/target/make_targets/{target_db}/{target_db}.index",
        target_path=directory("results/{database}/target/make_targets/{target_db}/")
    params:
        split=config["mmseqs2"]["map"]["split"],
        split_memory_limit=config["mmseqs2"]["map"]["split_memory_limit"]
    conda: "../envs/target.yml"
    log: "logs/{database}/target/make_targets/{target_db}.log"
    threads: config["mmseqs2"]["createdb"]["threads"]
    shell:
        """
        mmseqs createdb {input} {output.target_path}/{wildcards.target_db} --dbtype 1 2> {log}
        mmseqs createindex {output.target_path}/{wildcards.target_db} {params.split} {params.split_memory_limit} /tmp 2> {log}
        """
