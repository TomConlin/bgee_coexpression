#! /usr/bin/gawk -f

# uberon_subclass_class.tab  tissue_index.txt
BEGIN{FS="\t"}

FNR==NR {parent_of[$1] = $2}
FNR!=NR {
	tissue[FNR] = leaf = step = $1
	while(step in parent_of){
		path[leaf] = path[leaf] "," step
		step = parent_of[step]
	}
}
END{
	for(i=1; i<FNR; i++){
		split(path[tissue[i]], left_path, ",")
		left_len = length(left_path)
		for(j=i+1; j<=FNR; j++){
			split(path[tissue[j]], right_path, ",")
			right_len = length(right_path)
			for(ii=2; ii<=left_len; ii++){
				for(jj=2; jj<=right_len; jj++){
					if( left_path[ii] == right_path[jj] && \
					    left_path[ii] != tissue[i] && \
					    left_path[ii] != tissue[j] ){
						print tissue[i] "\t" tissue[j] "\t" left_path[ii]
						# ii += left_len; jj += right_len  # bail
					}
				}
			}
		}
	}
}


