#! /usr/bin/gawk -f

# called as
# ./ancestor_path.awk uberon_subclass_class.tab expression_profile_lables.txt
# 


BEGIN{FS="\t"}

FNR==NR {parent_of[$1] = $2}

FNR!=NR {name[$1] = $2}

END{
	for(node in name){
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
