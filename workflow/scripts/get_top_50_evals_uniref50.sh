#touch $1
awk '$3>50 {{print}}' $3 > $2
awk -F\\t '{print>$1}' $1
for f in (ls $4*_1)
do
     prefix=$(basename $f)
     python workflow/scripts/get_top_hits.py -i $PWD/$f -o $PWD/$4$prefix"_top50.tsv"
done
echo query     target  pident > $2
cat *_top50.tsv >> $2
