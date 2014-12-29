file=$1



while read line; do
	#echo $line
	echo $line > line
	
	chr=`cat line | awk '{print $1}'` 
	germ=`cat line | awk '{print $2}' |sed "s/,/\n/g" | sort | uniq | grep -v void | wc -l` 
	tumor=`cat line | awk '{print $3}' |sed "s/,/\n/g" | sort | uniq | grep -v void | wc -l` 
	rm line
	echo $chr" "$germ" "$tumor
done < $file
