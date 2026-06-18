#awk -F',' 'NR > 1 { gsub(/^ +| +$/, "", $25); species[$25] = $2 }
#          END { for (s in species) print s "\t" species[s] }' \
#    /fs/project/bradley.720/db/phylogenize/v2.0/human-gut/human-gut-taxonomy.csv > species_cluster.tsv

awk -F',' '
  NR==1 {
    for (i=1; i<=NF; i++) {
      if ($i == "species") sp_col = i
      else if ($i == "cluster") cl_col = i
    }
    next
  }
  NR>1 {
    gsub(/^ +| +$/, "", $sp_col)
    gsub(/^ +| +$/, "", $cl_col)
    if ($sp_col != "" && $sp_col != "\"\"")   # skip empty species
      print $sp_col "\t" $cl_col
  }
' /fs/project/bradley.720/db/phylogenize/v2.0/human-gut/human-gut-taxonomy.csv > species_cluster.tsv
# Step 2: Loop through each species/cluster and process FASTA
while IFS=$'\t' read -r species cluster; do
  fasta_species="s__${species}"

  awk -v sp="$fasta_species" -v cl="$cluster" '
    BEGIN { keep = 0; found = 0 }
    /^>/ {
      keep = index($0, sp) > 0
      if (keep) {
        found = 1
        header = $0
        # Step 1: Extract accession (first word after ">")
        match(header, /^>\S+/)
        acc = substr(header, RSTART + 1, RLENGTH - 1)

        # Step 2: Extract species name after s__
        match(header, /s__[^[]+/)
	species = substr(header, RSTART + 3, RLENGTH - 3)
	gsub(/^[ \t]+|[ \t]+$/, "", species)
        # Step 3: Construct clean header
        print ">" acc ";;" species";;" cl
        next
      }
    }
    {
      if (keep) print
    }
    END { if (!found) exit 1 }
  ' /fs/project/bradley.720/db/gtdb/v226.0/genomic_files_all/ssu_all_r226.fna >> filtered.fasta

  if [ $? -ne 0 ]; then
    echo "No match found for: $species"
  fi

done < species_cluster.tsv
