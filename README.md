# Generation of Protein-Clustered Pangenome Databases
For looking a species level pangenomic comparisons, it may be of use to have protein-clustered pangenomes in a taxa matrix. It is also easier to manage a matrix of clusters than individual pangenomes; hence why using this while developing pangenomic analysis tools is useful. [Phylogenize](https://bitbucket.org/pbradz/phylogenize) is a tools that allows, as descibed by it's developers, "links genes in microbial genomes to either microbial prevalence in, or specificity for, a given environment, while also taking into account an important potential confounder: the phylogenetic relationships between microbes". Protein level databases are used for the tool, which is why it is important as a compliment to have a way to easily generate new databses to increase its reach easily. We develope this workflow to efficiently work with nucleotides pangenomes to create new databases. We also encourage the community to contribute to this effort by submitting PR requests for databases to include or submitting the final databases generated via this workflow to the developers of Phylogenize or this repository.

## Running The Workflow
This workflow expects that you have conda installed prior to starting. Conda is very easy to install in general and will allow you to easily install the other dependencies needed in this workflow. Then you will need to download this repository via git clone such as `git clone git@github.com:Kekananen/phylogenize-db-prep.git`.

1. Edit the config/pepconfig.yml file's `raw_data:` string to be where your files are located at. Make sure to use the full path to avoid any errors. All your files will need to be in the same directory for this workflow. If you have a lot of files, you can symlink them into one directory to avoid taking up any more space. Another assumption if that the files will end with `.ffn`. The files won't be seen if they don't end with this; however you can rename then with the symlink when generated to avoid editing any actual file names prior to running this workflow. 
2. 
3. 
4.

### Dependencies and Installation
## Output Generated
## Submission Process
Please contact either [Kathryn Kananen](kananen.13@osu.edu) or [Patrick Bradley](bradley.720@osu.edu) if you wish to submit to the phylogenize databases being used.
