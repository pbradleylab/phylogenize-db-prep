include: "clustering.smk"
include: "mapping.smk"

rule mmseqs2_convertalis_unmapped_blast_uhgp50_db:
     input:rules.mmseqs2_createsubdb_uhgp50.output.outdir
     output:"results/{database}/uhgp50/mmseqs2/convert2fasta/{database}_convertlis.txt"
     params:
         prefix=rules.mmseqs2_createsubdb_uhgp50.params.prefix,
     conda: "../envs/blast.yml"
     shell:
         """
         mmseqs convert2fasta {input}/{params.prefix} /tmp/tmp.fa
         grep '>' /tmp/tmp.fa | sed 's/>//g' > /tmp/tmp.3
         paste /tmp/tmp.3 /tmp/tmp.3 > {output}
         """
