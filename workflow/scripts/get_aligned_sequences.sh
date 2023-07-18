grep ">" $1 > /tmp/headers.0
awk '(NR==FNR) { remove[$1]; next } /^>/ { p=1; for(h in remove) if ( h ~ $0) p=0 }p' /tmp/headers.0 $2 > $3
