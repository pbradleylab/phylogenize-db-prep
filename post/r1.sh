echo $(cut -f1 binary.header) > binary.all
cat protein_all.txt >> binary.all
# Generate the skeleton for the binary with the headers. Stopped at 100111 (46GB)
for ((i=2; i <= 111; i++))
do
	echo $(cut -f$i binary.header) > /tmp/col.tmp
	yes 0 | head -n215352134 >> /tmp/col.tmp
	paste binary.all /tmp/col.tmp > /tmp/binary.tmp
	mv /tmp/binary.tmp binary.all
	rm /tmp/col.tmp
done
cp binary.all binary.template

awk 'NF{NF--};1' binary.all > binary.template
mv binary.template rounds/r1/
cp rounds/r1/binary.template rounds/r1/binary.copy
for l in $(cat rounds/r1/species.txt); do 
  grep "^$l" phylogenize_sorted.protein >> rounds/r1/phylogenize_sorted.protein
done
cut -f2 rounds/r1/phylogenize_sorted.protein | sort -u > rounds/r1/proteins.txt


# Search for the protein with the largest number of shared hits first and then remove
# those entries from the file that is being iterated through. Start with the largest
# number of hits found across species.
for l in $(cat proteins.txt); do
# Index of where the protein is found in the file
declare -a rows=($(grep -n $l rounds/r1/binary.template | sed 's/:.*//g'))
# Name of the protein found
declare -a proteins=($(grep -n $l phylogenize_sorted.protein | sed 's/.*://g' | cut -f2))
# Name of the species found
declare -a species=($(grep -n $l phylogenize_sorted.protein | sed 's/.*://g' | cut -f1))

len=${#rows[@]}
  for ((i=1; i <= $len; i++)); do
    specie=$(echo "${species[$i-1]}")
    row=$(echo "${rows[$i-1]}")
    protein=$(echo "${proteins[$i-1]}")
    value=1

    # Find the index of the species (header) for the column 
    col=$(head -n1 rounds/r1/binary.template | tr "\t" "\n" | grep -nx "$specie" | cut -d":" -f1)

    echo $specie; echo $row; echo $protein; echo $col
    awk -v value=1 -v row=$row -v col=$col 'NR==row {$col=value}1' rounds/r1/binary.template > rounds/r1/binary.tmp
    mv rounds/r1/binary.tmp rounds/r1/binary.template
  done
done
