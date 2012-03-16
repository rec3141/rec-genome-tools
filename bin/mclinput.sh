for file in `find ./NC* | grep 'hits$'`; 
    do cat $file; 
done
