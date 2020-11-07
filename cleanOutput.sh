

 cat $1 | sed 's/\"//g' > tt1
 cat tt1 | sed 's/,/COMMA/g' > tt2
 cat tt2 | sed 's/-/HYPHEN/g' > tt3
 cat tt3 | sed "s/\'//g" > tt4
 mv tt4 $1
rm tt1
rm tt2
rm tt3

