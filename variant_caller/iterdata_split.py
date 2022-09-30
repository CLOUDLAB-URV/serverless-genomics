import ast
concur_fun=1000
fasta_set=147

iterdata_file="iterdata_test.txt"
iterdata_sets=[]
iterdata_set=[]
with open(iterdata_file, 'r') as f:
    s = f.read()
    s = ast.literal_eval(s)
    iterdata_n=len(s)
    print("iterdata_n "+str(iterdata_n))
    
    count=0
    count_tot=0
    for el in s:
        count+=1
        count_tot+=1
        fastq_count=el["fastq_chunk"][1]["number"]
        if ((count+fasta_set)>concur_fun and fastq_count==1):
            print("starting new iterdata set")
            print("count total: \t"+str(count_tot)+"\tcount: "+str(count)+"\tfastq_count "+str(fastq_count))
            print("length of finished set: "+ str(len(iterdata_set)))
            iterdata_sets.append(iterdata_set)
            iterdata_set=[]
            count=1
        iterdata_set.append(el)
        if count_tot==iterdata_n:
            iterdata_sets.append(iterdata_set)

    print("number of iterdata sets: "+str(len(iterdata_sets)))
    print("size of each iterdata set: ")
    for el in iterdata_sets:
        print(str(len(el)))

        #print("el: "+str(count)+"\t"+str(el))




