using StatsBase
using SharedArrays
include("macros_TEC_minimal.jl")

"""
    read data files 
"""
function init(tissue::String, gene::String, genetissue::String,)
	 tissue_array = readdlm(tissue, '\n') 
	 geneid_array = readdlm(gene, '\n')
	 data = readdlm(genetissue ,'\t', String)
	 tissue_index = Dict(zip(tissue_array, 1:length(tissue_array)))
	 geneid_index = Dict(zip(geneid_array, 1:length(geneid_array)))
	#  rotating the file make attributes contigious in memory
	return tissue_index, geneid_index, geneid_array, rotr90(data)
end 

"""
    Collapse expresion into a profile vector
"""
function build_profile(tissue_index, geneid_index, data) 
	profile = BitArray(zeros(length(tissue_index) ,length(geneid_index)))
	for i in 2:length(data[1,:])
		if data[3,i] == "present"
			profile[
				tissue_index[data[2,i]], 
				geneid_index[data[1,i]]
			] = 1
		end
	end
	return profile
end


"""
    wondering if there are any empty profiles
    and what the typical expression tissue profile looks like
  
"""
function check_profile(profile::BitArray)
    sumstat = Array{Int32}(size(profile, 2))
    tisdist = zeros(Int32, size(profile, 1))
    empty = 0
    for g in 1:size(profile,2)
        tisdist .+= profile[:,g]
        sumstat[g] = sum(profile[:,g])
        if 0 >= sumstat[g]
           empty += 1
        end
    end
    println(STDERR, "Gene expression profile curve statistics\n")
    println(STDERR, "How many tissues per gene")
    describe(sumstat)
    println(STDERR, "----------------------------------------------------") 
    println(STDERR, "How many genes per tissue")
    describe(tisdist)
    println(STDERR, "Flat profiles:\t", empty)
end

"""
   where the time goes 
"""
function pairwise_comparison(geneid_array, geneid_index, profile, threshold::Float64) 
	# for parallel writable array only simple numeric arrays are allowed
     elbow_room = 5000000
	 assoc_len = elbow_room * nprocs()
 	 g_len = length(geneid_index)

    association = SharedArray{UInt16}(3, assoc_len)	

	println(STDERR, now(), "\t Comparison Starting")
	println(STDERR, "Using:\t", nprocs(), " processes")

	# An explicit balanced work load option would be nice 
	# one approach could be @parallel_LT @parallel_UT anotations
    # for Lower|Upper Triangular matrix
	# another more general approach	would allow passing a partition function

    mypart = 0
	@sync @parallelLT for a in 1:g_len-1 
		mystart = elbow_room*(myid()-2)
		if 0 == a % 10000 
			println( STDERR, 
                a, "\t", now(), "\t Process: ", myid(), "  Produced: ", mypart)
		end
		for b in a+1:g_len
			intersect = sum(profile[:,a] .& profile[:,b])
			together = sum(profile[:,a] .| profile[:,b])
			jaccard = intersect / together
			if threshold <= jaccard
				mypart+=1
				association[1,mystart + mypart] = a 
				association[2,mystart + mypart] = b
				association[3,mystart + mypart] = round(jaccard*100)
			end
		end
        if mypart > elbow_room
            println(STDERR, 
                "Woopsie! processor ", myid()," is over running its space at ", 
                 mypart)
                exit(-1)
        end
 	end
    println(STDERR, now(), "\t Comparison done")
	# try to eliminate the [0,0,0]s
	# [if(association[:,i]!=[0,0,0]) association[:,i] end for i=1:assoc_len]
 	# turns them into sweet nothings which is no better
    sm = SparseMatrixCSC(association)
    println(STDERR, "size of assoc as soarce matrix", size(sm))
    nz = nonzeros(sm)
    println(STDERR, "size of assoc as non zeros", size(nz))
    rs = reshape(nz,3,:)
    println(STDERR, "size of assoc as reshaped", size(rs))
	return rs
end



function report(association, geneid_array)
    # println(STDERR, "Associations: ", size(association,2))
	writedlm("bgee-association.raw", rotl90(association), '\t')
    reitute = [
        [geneid_array[association[1,i]] * "\t" * 
         geneid_array[association[2,i]] * "\t" * 
         string(association[3,i])] 
            for i = 1:size(association,2)
    ]
    writedlm("bgee-association.tab", reitute, '\n')

end


######################################################################

"""
    Wrapper function
"""
function main()
    println(STDERR,"Begin  \t", now())
    (tissue, gene, genetissue) = 
        ("tissue_index.txt", "geneid_index.txt", "geneid_bintissueid_exists.tab")

    (tissue_index, geneid_index, geneid_array, data) = 
        init(tissue, gene, genetissue)

    println(STDERR, "Loaded \t", now())

    profile = build_profile(tissue_index, geneid_index, data) # rotr90(data))
    println(STDERR, "Profile\t", now())
    println(STDERR, "####################################################")
    check_profile(profile)
    println(STDERR, "####################################################")
    # 40-ish minutes with 4 procs
    association = pairwise_comparison(geneid_array, geneid_index, profile, 0.97)
    println(STDERR, "AllPairs\t", now())
    println(STDERR, "####################################################")
    report(association, geneid_array)
    println(STDERR, "Dumped \t", now())
end 


# Call the wrpper
main()

###############################################################################
#

function pairwise_comparisonII(
    geneid_array, geneid_index, profile, threshold::Float64) 

    @everywhere  elbow_room = 4000000
	@everywhere  assoc_len = elbow_room * nprocs()
 	@everywhere  g_len = length(geneid_index)
    @everywhere association = SharedArray{UInt16}(3, assoc_len)	
    println(STDERR, now(), "\t Comparison Starting")
	println(STDERR, "Using:\t", nprocs(), " processes")
    mypart = 0
   
    for k in workers() 
        @spawn k worker()
     end
    println(STDERR, now(), "\t Comparison done")
	# try to eliminate the [0,0,0]s
	# [if(association[:,i]!=[0,0,0]) association[:,i] end for i=1:assoc_len]
 	# turns them into sweet nothings which is no better
    sm = SparseMatrixCSC(association)
    println(STDERR, "size of assoc as soarce matrix", size(sm))
    nz = nonzeros(sm)
    println(STDERR, "size of assoc as non zeros", size(nz))
    rs = reshape(nz,3,:)
    println(STDERR, "size of assoc as reshaped", size(rs))
	return rs
end


@everywhere function worker(myrange,  )
	 for a in myrange
		mystart = elbow_room*(myid()-2)
		if 0 == a % 10000 
			println( STDERR, 
                a, "\t", now(), "\t Process: ", myid(), "  Produced: ", mypart)
		end
		for b in a+1:g_len
			intersect = sum(profile[:,a] .& profile[:,b])
			together = sum(profile[:,a] .| profile[:,b])
			jaccard = intersect / together
			if threshold <= jaccard
				mypart+=1
				association[1,mystart + mypart] = a 
				association[2,mystart + mypart] = b
				association[3,mystart + mypart] = round(jaccard*100)
			end
		end
        if mypart > elbow_room
            println(STDERR, 
                "Woopsie! processor %i is over running its space at %i", 
                myid(), mypart)
        end
 	end
end
