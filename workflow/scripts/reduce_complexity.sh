indir=$1
out=$2
database=$3

mkdir -p results/$database/renamed/
mkdir -p results/$database/transeq/
mkdir -p results/$database/combined_fasta/

#for file in $(ls $indir)
#do 
#name=$(echo $file | sed 's/.ffn//g')
#bbrename.sh in=$indir/$name.ffn out=results/$database/renamed/$name.ffn prefix=$name addprefix=t
#transeq -clean true results/$database/renamed/$name.ffn results/$database/transeq/$name.ffn
#done
cat results/$database/transeq/*.ffn > results/$database/combined_fasta/$database.fa
