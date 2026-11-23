# Workflow run guide

The [README](../README.md) gives the installation steps and profile-build commands.
This guide describes the input contracts, configuration, optional branches, and
failure messages. Commands run from the repository root with the workflow Conda
environment activated.

## Input formats and matching

All files for a catalogue must describe the same release of the source data. The
catalogue key, called `{mapping_db}` in output patterns, identifies this set of files.
It is distinct from `database`, which names the run's output directory, and
`{target_db}`, which identifies a reference protein-family database.

### Protein catalogue and membership

The catalogue is an amino-acid FASTA of centroid sequences. Sequence IDs are the
first whitespace-delimited part of each FASTA header and must be unique. For example:

```text
>C001
MAVKQSTAAAG
>C002
MELKQSTVVAG
```

The membership file has exactly two tab-separated columns and **no header**. The
first column is a centroid ID from the catalogue; the second is an original member
protein ID. Include every member protein, including centroid proteins when they are
members of their own clusters:

```text
C001	G001.1_00001
C001	G002.1_00003
C002	G001.1_00002
```

These are real tab characters in the example. A centroid can have members from many
genomes. The genome that supplied the centroid does not determine the species
membership of the entire family.

`create_species_matrix.split_char` in `config/tools.json` separates the genome and
protein portions of a member ID. The default `_` extracts `G001.1` from
`G001.1_00001`. The scripts take the portion before the **first** occurrence of the
delimiter. If genome IDs contain underscores, use member IDs such as
`GCA_000001.1|00001` and set `split_char` to `|`. Changing only the configuration
cannot repair IDs that do not contain an unambiguous delimiter.

### Genome metadata

Supply a TSV with a header and at least these three columns:

```text
Genome	Species_rep	Lineage
G001.1	G001.1	d__Bacteria;p__Firmicutes_A;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Bacillus;s__Bacillus_example
G002.1	G001.1	d__Bacteria;p__Firmicutes_A;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Bacillus;s__Bacillus_example
```

- `Genome` identifies each source genome. IDs must be unique after a trailing
  numeric version, such as `.1`, is removed for matching.
- `Species_rep` is the canonical species representative. Every representative must
  have its own row with `Genome == Species_rep`.
- `Lineage` contains seven semicolon-separated ranks in domain, phylum, class,
  order, family, genus, and species order. Rank prefixes such as `p__` are removed
  from the prepared taxonomy. Underscores within names are preserved. An empty
  phylum becomes `unclassified`.

Member genome IDs are joined to `Genome` using version-normalized matching keys.
Output matrix columns and tree tips use the original canonical `Species_rep` values.
Missing or ambiguous IDs cause an error. Every source species must occur in the
member-derived profiles. For continuous profiles, every genome in the metadata
contributes to its species denominator, even when it lacks a particular family.

### Source trees

`files.tree.<mapping_db>` supplies one Newick tree with finite, nonnegative branch
lengths. Its tips must cover the canonical species representatives, allowing the
same trailing-version normalization for matching. Extra source-tree tips are
allowed; the output is restricted to the configured taxonomy.

An optional `files.alttree.<mapping_db>.tree` can cover species from another source
tree, for example archaeal representatives absent from a bacterial tree. Add this
object inside `files` for the example catalogue:

```json
"alttree": {
    "catalogue": {
        "tree": "/path/to/archaeal_species_tree.nwk"
    }
}
```

Each species must be covered by exactly one tree, and every phylum must be contained
within one tree. The workflow roots and RED-normalizes each complete source tree
separately before extracting per-phylum trees. The legacy alternate-tree `name`
setting does not rename the phylum outputs.

The direct `make_tree.R` CLI also supports `--allow_missing_species` with a CSV
containing a `cluster` column listing exactly the audited species absent from the
source trees. This option is not exposed by the Snakemake rule and leaves those
species absent from the tree RDS; it is not a general workaround for mismatched data.

### Reference databases

Provide an amino-acid FASTA for each `target_db.paths` entry. Its sequence IDs become
the target family IDs. The example uses one reference database; at least one entry
is required. Additional entries are searched in the order described below.

Reference databases are processed in JSON insertion order. Queries passing the
hard-coded identity filter for an earlier reference are removed before the next
search. The search minimum identity is 0.50, and the converted-hit filter retains
`pident > 50`. Within a reference, top-hit selection keeps the highest percentage
identity per query; equal values retain the first hit. Unmatched queries from the
last reference are clustered with Linclust. Its TSV is interpreted as
`representative,member`, so the member is the query and the representative is the
new family ID. This is a protein-family mapping threshold, not a species boundary.

Family IDs are stored as literal strings, including leading zeros. No reference-name
prefix is added to family IDs. Keep IDs unique across references and de novo
representatives if they denote different families. The workflow expects reference
FASTAs to exist; `target_db.urls` and `skip_db_download` do not implement reference
FASTA downloading in the active rules.

## Configuration reference

Use [config/config.example.json](../config/config.example.json) as a starting point.
The file is a template with placeholder paths, not bundled test data.

| Setting | Meaning |
| --- | --- |
| `database` | Run label used under `results/` and `logs/`. Use a new label for a different source/reference release. |
| `mode.annotate_only` | `false` builds profiles when requested; `true` makes the annotation branch reuse existing mapping intermediates. |
| `files.fasta.<mapping_db>` | Catalogue centroid protein FASTA. |
| `files.mapping.<mapping_db>` | Headerless centroid/member TSV. |
| `files.taxonomy.<mapping_db>` | Source genome metadata TSV. |
| `files.tree.<mapping_db>` | Primary Newick tree. |
| `files.alttree.<mapping_db>.tree` | Optional additional Newick tree. Omit the entry when unused. |
| `files.16S.<mapping_db>.generate_for_16S` | Whether full mode with `all` requests the 16S completion target. Defaults to `false`. |
| `files.16S.<mapping_db>.fasta_dir` | Directory traversed for genome GFF files when 16S is enabled. |
| `target_db.paths.<target_db>` | Reference family protein FASTA. At least one is required. |

Use simple run, catalogue, and reference labels containing letters, numbers, or
hyphens. All relative paths are resolved from the repository working directory;
absolute paths make input locations explicit. Avoid whitespace in paths because
some shell rules do not quote file arguments.

The workflow always reads `config/config.json` and `config/tools.json`. The README
recommends backing up the first file and replacing it with your edited example.
There is no need to add `indir`, `pepfile`, or PEP sample tables for the protein
catalogue route. `indir` is required only when explicitly requesting
`reduce_complexity`; the legacy PEP-dependent translation route is not wired into
the default catalogue build.

If you use `--configfile` for overrides, Snakemake recursively merges dictionaries
with the default files. For example, adding a `catalogue` entry does not delete the
shipped `marine` entry. Review and replace the default `files` and `target_db.paths`
dictionaries when changing datasets. Mode overrides can use JSON booleans or
`"True"`/`"False"` strings. An explicit command-line configuration mode overrides the
default mode. Within the same configuration, `mode.annotate_only` takes precedence
over the legacy `annotation.annotate_only` key.

Tool settings remain in [config/tools.json](../config/tools.json):

| Setting | Effect |
| --- | --- |
| `create_species_matrix.split_char` | Member protein-to-genome delimiter; must not be empty. |
| `mmseqs2.createdb.threads` | CPU request for query/target database preparation. |
| `mmseqs2.map.threads`, `mmseqs2.map.sensitivity` | Search CPU request and MMseqs2 sensitivity arguments. |
| `mmseqs2.convertalis.threads` | CPU request for conversion rules. |
| `mmseqs2.linclust.threads`, `mmseqs2.linclust.tmp_dir` | Linclust CPU request and scratch directory. |
| `anvio.contigs_db.threads`, `anvio.run_kegg_annotations.threads` | CPU requests for the KO branch. |
| `emapper.threads`, `emapper.go_evidence` | CPU request and GO evidence selection for eggNOG-mapper. |

Snakemake caps declared threads at the available `--cores`. Not every tool uses its
rule's thread count explicitly. The legacy `mmseqs2.linclust.seq_id_precent` setting
is not passed to Linclust; editing it does not change the clustering threshold.
The continuous rule currently fixes its request at 8 CPUs and 64,000 MB and passes
`--memory 16 --speedy TRUE` to its R script.

## Environment versions and reuse

The launcher selects the libmamba solver and strict channel priority when activated.
The environment files use conda-forge followed by Bioconda and exclude implicit
default channels with `nodefaults`. The tracked MIDAS 2 test environment also uses
the package author's `zhaoc1` channel. Versions were checked against Linux package
metadata and solved with strict channel priority on 2026-11-23.

| Environment | Selected versions and compatibility constraints |
| --- | --- |
| Root `env.yml` | Python 3.14.10, Snakemake 9.27.2, interface-common 1.23.0, Conda 26.7.3, eido 0.2.6, peppy 0.40.8, and pandas below 3. Newer Conda and interface-common releases require `packaging >=26`, while this Snakemake release requires `packaging <26`. |
| Query, target, mapping, clustering, and target-download environments | MMseqs2 20.ae60e throughout, so database creation, search, and conversion use the same version. Mapping also supplies samtools 1.25; clustering supplies `faSomeRecords` 502. |
| Matrix | Polars 2.1.0, R 4.5.3, Matrix 1.7-6, castor 1.8.8, and duckplyr 1.2.1. The current compatible R package builds select R 4.5. |
| Anvi'o annotation | The amino-acid-contig fork requires Python 3.10, NumPy 1.24.1, pandas 1.4.4, matplotlib 3.5.1, scikit-learn 1.2.2, and PuLP 2.7.0. Setuptools 81 retains its required `pkg_resources` module. Those constraints are retained; the surrounding packages use current compatible releases. |
| eggNOG/GO | eggNOG-mapper 2.1.15 and GOATOOLS 1.6.9. GOATOOLS is installed from PyPI because Bioconda's older release requires an Excel reader incompatible with current pandas. DIAMOND below 2.1 requires the older zlib ABI, so this separate environment uses compatible wget 1.21.4 and pandas 2.3.3. Its transitive MMseqs2 installation is separate from the main workflow's MMseqs2 databases. |
| 16S | PASTA 1.9.3, FastTree 2.2.0, vsearch 2.32.0, Biopython 1.88, and GFF3 1.0.1, plus the R packages imported by the filtering scripts. Tree-fixing, rescaling, and pruning rules activate this environment too. |
| Phylogenetics | Phylogenize 0.91 and optparse 1.7.3 preserve the existing `install.data.figshare()` setup API and its R 4.1 package builds. Phylogenize 2.0.4 replaces that downloader with `install_data(path)`, which requires a local archive; migrating the setup rule is a separate change. |
| Translation | EMBOSS 6.6.0, its current release, and BBMap 40.05. |
| Legacy tests | MIDAS 1.3.2 and MIDAS 2 1.1.3 remain the latest stable releases of their respective tools. MIDAS 2 release candidates are not selected. The read-processing environment uses the current Snakemake, fastp, Kraken2, and Bracken packages. |

Anvi'o and eggNOG are separate because their Python and compiled-library
requirements differ. The Anvi'o post-deployment hook installs the exact custom
revision `21aad5f33a2c3c2c1e48f428dda08a0ef86673c0` without prompts or dependency
replacement; its Python and build dependencies are supplied by `annotation.yml`.
The hook checks that Anvi'o imports successfully before finishing. It runs only
when Snakemake creates that environment. The workflow relies on the fork's commands
and `--allow-amino-acid-contig-db` option; moving to a public Anvi'o release requires
validating that interface separately.
Its Python 3.10 dependency stack includes Snakemake 7.32.4; the workflow launcher
uses Snakemake 9.27.2 from the root environment.
Creating `annotation.yml` with Conda alone prepares the dependencies; Snakemake
runs the post-deployment hook that installs the custom Anvi'o package.

Install the root environment once, then activate it before running the workflow.
For an existing launcher environment, update it from the base environment instead
of attempting to create the same name again:

```bash
CONDA_SOLVER=libmamba CONDA_ALWAYS_YES=true conda env update --name phylogenize-pipeline --file env.yml --prune
conda activate phylogenize-pipeline
```

Use the actual environment name if your launcher was created with a different one.
Snakemake caches rule environments under `.snakemake/conda/` and reuses them across
runs. Keep that directory to avoid repeating installation. Changes to YAML files
or post-deployment hooks trigger a new environment; a failed installation must be
retried. Package and reference-data installation are separate: the annotation
resource rules retain their downloads under `resources/`.

After configuring the input files, software environments can be prepared before
starting the analysis. In full mode, this creates the environments needed by the
`all` target and runs their post-deployment hooks without executing analysis jobs:

```bash
snakemake all --snakefile workflow/Snakefile --cores 8 \
    --resources mem_mb=64000 --use-conda --conda-create-envs-only
```

Normal runs then reuse those environments. This command prepares software;
annotation reference-data downloads are still handled by the resource rules.

After an environment upgrade, use a new `database` label for a fresh build or inspect
Snakemake's proposed reruns carefully. Existing results and MMseqs2 database files
can have been produced with older tools. Do not mix old database indexes with a
new MMseqs2 executable manually.

## Workflow modes

### Profiles and source trees

Set `mode.annotate_only` to `false`, then request the three file targets shown in the
README. This builds searches, de novo clusters, species/family mappings, binary and
continuous matrices, and source-tree subsets. It does not schedule the annotation
or 16S branches. Each catalogue needs its own three output targets.

### Full aggregate target

With `mode.annotate_only: false`, this command plans the full run:

```bash
snakemake all --snakefile workflow/Snakefile --cores 8 \
    --resources mem_mb=64000 --dry-run
```

The `all` target requests binary profiles, continuous profiles, source-tree RDS,
species-linked KO annotations, and GO annotations for every catalogue/reference
combination. If 16S is enabled for a catalogue, it also requests the final 16S
alignment and tree through the completion marker. Use the explicit `all` target:
the first included rule is not this aggregate rule.

Read the known limitations below before replacing `--dry-run` with
`--use-conda --rerun-incomplete`. The full aggregate target is not currently a
validated end-to-end annotation run.

### Annotation only

Set `mode.annotate_only` to `true` in the same run configuration, then use:

```bash
snakemake all --snakefile workflow/Snakefile --cores 8 --dry-run
```

Annotation-only mode includes no sequence-search, clustering, profile, or 16S
producer rules. It requires these files from the original mapping run:

| Intermediate | Path relative to `results/{database}/` |
| --- | --- |
| Combined centroid/family TSV | `binary/combined_species_hits/{mapping_db}_{database}.tsv` |
| Species/family CSV | `binary/get_taxonomy/{mapping_db}-binary-temp.csv` |
| De novo representatives | `clustering/mmseqs2_linclust/{mapping_db}/unaligned_linclust_rep_seq.fasta` |
| Full alignment for each reference | `clustering/mmseqs2/convertalis_full_blast/{target_db}_{mapping_db}_convertlis.8` |

Keep the original `database`, catalogue/reference keys, reference FASTAs, and mapping
intermediates together. Missing files cause an `annotate_only requires an existing
mapping intermediate` error. Earlier builds that deleted the species/family CSV
must regenerate it in full mode. Existing annotation outputs are reused by
Snakemake when up to date; selecting annotation-only mode does not force a rerun.
This mode may still download annotation resources and run annotation tools.

### Optional 16S

Set `files.16S.<mapping_db>.generate_for_16S` to `true` and supply `fasta_dir` to
request this branch with full-mode `all`. The extractor recursively searches for
`.gff` and `.gff.gz` files with embedded FASTA and SSU annotations from barrnap or
INFERNAL. Genome file basenames are matched to the metadata `Genome` values.
The setting is ignored in annotation-only mode. The example keeps it disabled
because of the tool, path, and filtering problems below.

## Output reference

The README lists the profile, taxonomy, and combined-hit outputs. There is no final
packaging directory or automatic installation into Phylogenize. Preserve the run
configuration and record catalogue/reference releases alongside the results.

Binary RDS files contain logical sparse matrices. Continuous RDS files contain
sparse numeric matrices of fractions in `[0, 1]`. In both, a zero means the family is
absent, rows name families, and columns name canonical species. The continuous CSV
includes `column0` (family ID), `cluster` (species ID), taxonomy ranks,
`genomes_with_family`, `total_genomes`, and `frac_observed`. The denominator includes
all source genomes assigned to that species.

Inspect the example's output in R with `Matrix` and `ape` installed:

```r
library(Matrix)
library(ape)
binary <- readRDS("results/example/binary/get_binary/catalogue-binary.rds")
continuous <- readRDS("results/example/binary/get_binary/get_continuous/catalogue-continuous.rds")
trees <- readRDS("results/example/binary/get_tree/catalogue-tree.rds")
names(binary)
lapply(binary, dim)
phylum <- names(binary)[1]
head(rownames(binary[[phylum]]))
head(colnames(binary[[phylum]]))
head(trees[[phylum]]$tip.label)
```

Match matrices and tree tips by their IDs before combining them. Family row sets
can differ between phyla; the workflow does not promise identical list or column
ordering across every output.

Additional outputs, when their branches succeed, are relative to
`results/{database}/`:

| Output | Path |
| --- | --- |
| Species-linked KO annotations | `annotation/link_nodes/{mapping_db}/node_annotations.tsv` |
| KO completion marker | `annotation/run_anntoation/{mapping_db}.done` |
| GO prediction triples before term annotation | `annotation/parse_blast/{target_db}_{mapping_db}.tsv` |
| GO term annotations | `annotation/reannotate_blast_results/{target_db}_{mapping_db}.tsv` |
| Final 16S alignment | `16S/{mapping_db}-final/{mapping_db}-alignment.aln` |
| Final 16S tree | `16S/{mapping_db}-final/{mapping_db}-tree.phy` |
| 16S completion marker | `16S/.{mapping_db}.complete` |

The spellings `run_anntoation` and `convertlis` match the current rules. Annotation
resources live under `resources/anvio/kofams/`, `resources/eggnog/`, and
`resources/go_ontology/`. They are shared across runs in the same checkout. The
annotation environment's post-deploy script installs a pinned custom Anvi'o revision that
supports amino-acid contig databases and requires network access.

## Known limitations

These are implementation issues in the current rules, not additional setup steps.
A dry run can succeed without exercising them.

| Area | Current issue and impact |
| --- | --- |
| Shared mapping scratch | `combine_species_hits` appends to `/tmp/{mapping_db}.1` before deduplication. An interrupted run or concurrent run using the same catalogue key can contaminate assignments. Confirm no active job is using the scratch file before dealing with leftovers; changing only `database` does not isolate this path. |
| Empty unmapped set | `faSomeRecords` uses `grep` under strict shell handling, and Linclust has no explicit empty-input branch. A run in which all queries map can stop instead of producing an empty de novo result. |
| GO resource download | `download_go_ontology` reverses the `wget -O` arguments and uses `mkdir` without `-p`; resource setup can fail. |
| GO term annotation | `reannotate_blast_results` calls `workflow/scripts/annotate_go_terms.py`, which is absent from the tracked repository. A clean checkout cannot complete that rule. |
| GO filtering | `make_go_db` filters every line containing `-`, which can also discard valid identifiers. |
| KO merging | `combine_link_info` uses `{input}[0]` instead of indexing the first input file, so header extraction can fail. |
| Optional 16S | Several rules invoke `scripts/...` from the repository root instead of `workflow/scripts/...`. The taxonomy filter also needs review for missing matches and duplicate filtering. The environment supplies the required tools, but those path and processing issues still need fixes. |
| Raw genomes | The pangenome include is disabled, PEP initialization is commented out, and there is no active raw-genome aggregate path. Use an existing protein catalogue. |

The example and README profile commands avoid scheduling GO, KO, and 16S. They still
use the mapping and Linclust rules, so the scratch and empty-input limitations apply.
The full mode and annotation-only graphs have regression tests; production-scale
execution of every external tool and database download has not been validated by
those tests.

## Troubleshooting

| Message or symptom | What to check |
| --- | --- |
| Unexpected catalogue/reference or paths from the original configuration | Nested `--configfile` merges kept default keys. Replace `config/config.json` with the complete intended run configuration after making a backup. |
| `Target rules may not contain wildcards` before a build starts | Provide `all` or explicit output paths and `--snakefile workflow/Snakefile`. |
| Missing member genomes or ambiguous genome IDs | Check the member-ID delimiter, metadata release, and version-normalized `Genome` uniqueness. |
| Missing species representatives | Every `Species_rep` needs its own metadata row and member-derived species/family pairs. |
| Centroids missing from combined hits | Check catalogue/member ID agreement and the mapping/de novo results. Do not substitute centroid genomes for member genomes. |
| Missing or overlapping tree species | Check the primary and alternate trees against canonical `Species_rep` IDs. A species or phylum must not span source trees. |
| Old four-column binary input | Regenerate the taxonomy step; the current species/family CSV uses `cluster,target` columns. |
| Annotation-only missing intermediate | Restore the original mapping outputs or switch to full mode and request the needed profile outputs. |
| Incomplete output after an interrupted run | Resolve the failing tool/input problem and rerun the same targets with `--rerun-incomplete`. |
| Out of memory or insufficient memory resource | Allocate enough physical memory and set Snakemake's memory limit to at least the continuous rule's 64,000 MB request. The 16 GB DuckDB cap does not cap total R memory. |

Use `--printshellcmds` to see commands. Rule logs are under `logs/{database}/` when
a rule declares a log; otherwise inspect Snakemake's captured console output and
`.snakemake/log/`. An interrupted combine step also needs the scratch-file check
above before resuming. Keep retained mapping intermediates if you plan to annotate
later; do not mark them as temporary merely to save space.

## Regression checks

Run the Python/top-hit and workflow-mode checks in the workflow environment:

```bash
python -B -m unittest discover -s test -p test_get_top_hits.py -v
python -B -m unittest discover -s test -p test_workflow_modes.py -v
```

The mode checks create temporary inputs and run Snakemake dry runs; they do not run
MMseqs2 or annotation tools. `test/test_clustering.py` also checks database creation,
Linclust representative/member output, and literal FASTA IDs when MMseqs2 and
`faSomeRecords` are available on `PATH` alongside Python. To run the binary,
continuous, and tree checks, create
the matrix environment and provide its executables to the test process:

```bash
conda env create --yes --name phylogenize-matrix --file workflow/envs/matrix.yml
conda run --no-capture-output --name phylogenize-matrix \
    python -B -m unittest discover -s test -p test_binary_matrices.py -v
conda run --no-capture-output --name phylogenize-matrix \
    python -B -m unittest discover -s test -p test_continuous_pangenomes.py -v
conda run --no-capture-output --name phylogenize-matrix \
    python -B -m unittest discover -s test -p test_make_tree.py -v
```

The matrix environment supplies Python/Polars and R with the matrix/tree packages.
Tests skip when required executables are unavailable; inspect the reported skips
before treating a run as validation. The older `test/taraOceans/Snakefile` is a
separate data-generation workflow, not the main entrypoint or a self-contained
example dataset.
