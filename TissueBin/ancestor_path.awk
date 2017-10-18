#! /usr/bin/gawk -f

# uberon_subclass_class.tab

function ancestor(up, here, path){
	leaf = here
	there = here
	while(there in up){
		here = there
		path[there] = count+=1
		there = up[here]
	}
}

BEGIN{FS="\t"}

FNR==NR {parent_of[$1] = $2}

FNR!=NR {name[$1] = $2}

END{
	for(node in name){
	  	# printf("%s",node )
		there = node
		while(there in parent_of){
			here = there
			there = parent_of[here]
			if(there in name){
				printf("%s\t%s\t#  %s\t%s\n",node, there,name[node],name[there])
			}
		}
	}
}
