readonly tab=$'\t'
sort -u phylogenize.0 > phylogenize_uniq.0

cut -f1 phylogenize_uniq.0 > col1.txt
cut -f2 phylogenize_uniq.0 > col2.txt
sort -u col1.txt > col1_uniq.cnts
sort -u col2.txt > col2_uniq.cnts

# This is necessary to have items in the same order for the binary
sort -k2 -n phylogenize_uniq.0 > phylogenize_sorted.protein

# Format the output of the species so that the most 
grep '      [0-9]' col2_uniq.cnts | sed 's/      .//g' | sed 's/ //g' > species_1.hits
grep '     [0-9][0-9]' col2_uniq.cnts | sed 's/     ..//g' | sed 's/ //g' > species_2.hits
grep '    [0-9][0-9][0-9]' col2_uniq.cnts | sed 's/    ...//g' | sed 's/ //g' > species_3.hits
grep '   [0-9][0-9][0-9][0-9]' col2_uniq.cnts | sed 's/   ....//g' | sed 's/ //g' > species_4.hits
grep '  [0-9][0-9][0-9][0-9][0-9]' col2_uniq.cnts | sed 's/  .....//g' | sed 's/ //g' > species_5.hits
cat species_*.hits > protein_all.txt
grep ' [0-9][0-9][0-9][0-9][0-9][0-9]' col1_uniq.cnts | sed -r "s/.*( .*?)$/\1/" > species_all.txt

# Make the protein column for the binary
cut -f1 binary.header > protein_col.txt
cat species_*.hits >> protein_col.txt

# Make the first binary column with only 0's
yes 0 | head -n215352134 > blank.template
cp blank.template blank.0

# Generate the header for the file
header=""
for h in $(cat col1_uniq.cnts)
do
	header="${header}_split_${h}"
	echo $header
	echo $header > binary.header
done
sed -i 's/_split_......_split_/\t/g' binary.header
sed -i 's/_split_....._split_/\t/g' binary.header
sed -i 's/_split_...._split_/\t/g' binary.header
sed -i 's/_split_..._split_/\t/g' binary.header
sed -i 's/_split_.._split_/\t/g' binary.header
sed -i 's/_split_._split_/\t/g' binary.header

bash r1.sh
