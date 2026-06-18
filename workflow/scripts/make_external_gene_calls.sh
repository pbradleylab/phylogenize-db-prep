#!/bin/bash

fasta_file=$1
output_file=$2
counter=1

# Write header
echo -e "gene_callers_id\tcontig\tstart\tstop\tdirection\tpartial\tcall_type\tsource\tversion\taa_sequence" > "$output_file"

# Extract sequence names, lengths, and sequences per contig
awk '/^>/ {
    if (seq) {
        seq_length = length(gensub(/\s/, "", "g", seq))
        print counter++ "\t" name "\t0\t" seq_length "\tf\t0\t1\tGeneCaller\tv1.0\t" seq
    }
    name=$0
    gsub(/^>/, "", name)
    seq=""
    next
}
{
    seq = seq $0
}
END {
    if (seq) {
        seq_length = length(gensub(/\s/, "", "g", seq))
        print counter++ "\t" name "\t0\t" seq_length "\tf\t0\t1\tGeneCaller\tv1.0\t" seq
    }
}' "$fasta_file" >> $2
