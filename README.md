# Generation of Protein-Clustered Pangenome Databases

This Snakemake workflow prepares species-level gene-family profiles and phylogenetic
trees for [Phylogenize](https://bitbucket.org/pbradz/phylogenize). It starts with an
existing clustered protein catalogue, its protein membership table, genome metadata,
and a source tree. Catalogue proteins are matched to reference protein families;
unmatched proteins are clustered with MMseqs2 Linclust. The resulting family
assignments are used to build binary and continuous profiles for each phylum.

The workflow also contains KO/GO annotation and optional 16S branches. KO and GO
annotations retain the family and species identifiers used in the profiles, including
de novo families. The optional 16S branch still has execution issues described in
the [run guide](workflow/documentation.md#known-limitations). A successful dry run
checks dependencies and configuration; it does not validate the external tools or
the biological contents of your inputs.
Raw-genome pangenome construction is not implemented in the active workflow.

## Installation

Clone the branch containing this workflow and run all commands from the repository
root:

```bash
git clone --branch dev https://github.com/pbradleylab/phylogenize-db-prep.git
cd phylogenize-db-prep
```

Install Conda first, then create the workflow environment from the repository's
environment file:

```bash
CONDA_SOLVER=libmamba conda env create --yes --file env.yml
conda activate phylogenize-pipeline
snakemake --version
```

The launcher includes Snakemake, Conda, `peppy`, and `pandas`. The `--use-conda` option
automatically creates separate environments from `workflow/envs/` for the other
tools. Successful installations are reused on later runs; changing an environment
definition creates a new environment. The first run needs network access and disk
space for packages and reference data. See
[environment versions and reuse](workflow/documentation.md#environment-versions-and-reuse)
for the compatibility constraints and installation checks.

## Prepare the inputs

For each catalogue, provide the following files. Use uncompressed FASTA files and
paths without whitespace, since several shell rules do not quote their arguments.

| Input | Required format |
| --- | --- |
| Clustered catalogue | Amino-acid FASTA of protein centroids, with unique sequence IDs. |
| Protein membership | Headerless, two-column TSV: catalogue centroid ID, then member protein ID. |
| Genome metadata | TSV with a header containing `Genome`, `Species_rep`, and `Lineage`. |
| Source tree | One Newick tree with branch lengths and tips matching the species representatives. An additional tree can cover a separate set of species. |
| Reference families | Amino-acid FASTA with the target family IDs used in the output matrices. Supply at least one reference database. |

Member protein IDs must encode the genome ID before a configured delimiter. For
example, with `split_char` set to `_`, `G001.1_00001` belongs to genome `G001.1`.
Genome IDs containing `_` need a different delimiter in the member IDs and a matching
`create_species_matrix.split_char` setting in `config/tools.json`.

See [input formats and matching](workflow/documentation.md#input-formats-and-matching)
for small examples and the rules for species representatives, taxonomy, and trees.
The PEP files and `.ffn` translation rules are not required for this catalogue route.

## Configure a run

Copy [config/config.example.json](config/config.example.json) and replace every
`/path/to/` path with your own:

```bash
cp config/config.example.json config/config.local.json
```

The example uses `database: "example"` for the result directory and `catalogue` as the
catalogue key. It selects full mode, disables 16S, and uses one reference database.
Keep the catalogue key identical across `files.fasta`, `files.mapping`,
`files.taxonomy`, and `files.tree`. Change the result paths in the commands below if
you rename either `example` or `catalogue`.

After editing, choose an unused backup filename, back up the existing default
configuration, and install your run configuration:

```bash
cp config/config.json config/config.backup.json
cp config/config.local.json config/config.json
python -m json.tool config/config.json > /dev/null
```

The workflow loads `config/config.json` and `config/tools.json` automatically. Using
`--configfile` with a different catalogue does not replace the defaults: Snakemake
merges nested dictionaries, which can leave unwanted catalogues or reference
databases active. Replacing the default configuration as shown above avoids this.
See the [configuration reference](workflow/documentation.md#configuration-reference) for
optional settings and tool parameters.

## Build the profiles and trees

Use explicit output targets to build the profiles without scheduling annotation or
16S. Start with a dry run:

```bash
snakemake \
    results/example/binary/get_binary/catalogue-binary.rds \
    results/example/binary/get_binary/get_continuous/catalogue-continuous.rds \
    results/example/binary/get_tree/catalogue-tree.rds \
    --snakefile workflow/Snakefile --cores 8 --resources mem_mb=64000 --dry-run
```

After inspecting the planned jobs, run them with Conda environments enabled:

```bash
snakemake \
    results/example/binary/get_binary/catalogue-binary.rds \
    results/example/binary/get_binary/get_continuous/catalogue-continuous.rds \
    results/example/binary/get_tree/catalogue-tree.rds \
    --snakefile workflow/Snakefile --cores 8 --resources mem_mb=64000 \
    --use-conda --rerun-incomplete
```

The continuous-matrix rule requests 8 CPUs and 64,000 MB of memory. Its DuckDB
aggregation is capped at 16 GB, but constructing the sparse matrices uses additional
memory. On a cluster, run Snakemake inside a job with enough allocated memory and
CPUs; `--resources` limits Snakemake's scheduling and does not allocate cluster
resources. Large catalogues and reference databases also need substantial disk space.

Always provide explicit file targets or `all`. There is no root `Snakefile`, and
running the entrypoint without a target does not select the intended aggregate rule.

## Outputs

Paths below are relative to `results/{database}/`; `{mapping_db}` is the catalogue
key, such as `catalogue` in the example.

| Output | Path |
| --- | --- |
| Binary profiles | `binary/get_binary/{mapping_db}-binary.rds` |
| Continuous profiles | `binary/get_binary/get_continuous/{mapping_db}-continuous.rds` |
| Continuous counts and fractions | `binary/get_binary/get_continuous/{mapping_db}-continuous-full.csv` |
| Per-phylum source trees | `binary/get_tree/{mapping_db}-tree.rds` |
| Canonical species taxonomy | `binary/get_taxonomy/{mapping_db}-taxonomy.csv` |
| Species/family pairs retained for annotation reuse | `binary/get_taxonomy/{mapping_db}-binary-temp.csv` |
| Combined centroid/family hits | `binary/combined_species_hits/{mapping_db}_{database}.tsv` |

Each profile RDS contains a named list of sparse matrices, one per phylum. Rows are
gene families and columns are canonical `Species_rep` IDs. Binary values indicate
presence in any member genome. Continuous values are the fraction of source genomes
in the species containing the family; multiple copies in one genome count once.
Tree RDS files contain a named list of RED-normalized trees for the same taxonomy
phyla. Genome versions are normalized for matching, while canonical species IDs and
taxonomy underscores are preserved in the outputs.

There is no automatic `results/{database}/final/` packaging step. Keep the source
configuration with the output files. See the [output reference](workflow/documentation.md#output-reference)
for annotation and 16S paths and an R inspection example.

## Annotation and optional 16S

Full mode (`mode.annotate_only: false`) with the explicit `all` target requests
profiles, trees, KO/GO annotations, and any enabled 16S outputs. Annotation-only mode
(`mode.annotate_only: true`) requests annotations from existing mapping intermediates
and stops if those intermediates are missing. It does not rebuild profiles or run
sequence searches. JSON booleans and the legacy `"True"`/`"False"` strings are accepted;
the legacy `annotation.annotate_only` key is also supported.

Before running these branches, read the
[mode instructions](workflow/documentation.md#workflow-modes),
[annotation output formats](workflow/documentation.md#annotation-handoffs), and
[known limitations](workflow/documentation.md#known-limitations). The annotation
completion marker requires both species-linked KO and GO files. GO evidence is
exported separately; downstream GO enrichment requires support in the analysis tool.

## Troubleshooting and development

Rerun the same command with `--rerun-incomplete` after fixing an interrupted job.
Check the failing rule, its printed shell command (`--printshellcmds`), and any
rule-specific files under `logs/{database}/`. Changing the output `database` name
keeps a new run separate from existing results; retain the same configuration when
resuming a run. See the [run guide](workflow/documentation.md) for matching errors,
annotation prerequisites, scratch-directory limitations, and regression checks.

## Submission Process

To discuss contributing a database or transferring large result files, contact
[Kathryn Kananen](mailto:kananen.13@osu.edu) or
[Patrick Bradley](mailto:bradley.720@osu.edu). Include the configuration, catalogue and
reference versions, and the generated profiles, taxonomy, and trees.
