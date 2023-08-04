files=$(find $1 -name $2)
cat $files > /tmp/tmp.1 
uniq /tmp/tmp.1 > $3
