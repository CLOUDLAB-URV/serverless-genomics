import Metadata
import subprocess as sp
import lithops

cmd = "sudo python3.8 FastaPartitioner.py --mybucket ayman-lithops-meta-cloudbutton-hutton --chunk_size 134217728 --overlap 300 --data_location s3://ayman-lithops-meta-cloudbutton-hutton/fasta/hg19.fa --fasta_folder fasta/"

total_objects = sp.check_output(cmd.split(" ")).decode('utf-8')
print(total_objects)

i = 0

storage = lithops.Storage()
data = ""
chr_title = ""
ini = 1
while i <= int(total_objects):
    
    ret = storage.get_object('ayman-lithops-meta-cloudbutton-hutton', f'cache/obj{i}.data', stream = True).read().decode('utf-8')
    for line in ret.splitlines():
        if line == '\n' or line == '':
            pass
        elif '-' in line and ',' in line:
            d = line.split(',')
            range = d[1].split('-')
            abs_pos = int(range[1]) -  int(range[0])+1
            range = '\t'.join(range)
            data =   data + d[0].replace('>','')+ '\t' +chr_title[1].replace('>','') + '\t' + str(ini) + '\t' +  str(abs_pos)  + '\n'
            ini = abs_pos 
        elif '-' in line:
            range = line.split('-')
            abs_pos = int(range[1]) -  int(range[0])+1
            data =  data + chr_title[0].replace('>','') + '\t' + chr_title[1].replace('>','') + '\t' + str(ini) + '\t' +  str(abs_pos) +'\n'
            ini = abs_pos
        else:
            chr_title = line.split(',')
            ini = 1


            
    
    i= i+1

print(data)



#Test of efetch api
#Api key needs to be added as an environment variable: export ENA_API_KEY=<your_api_key>
accessions = ['SRR6052133']
meta = Metadata.SraMetadata()
accesion = meta.efetch_sra_from_accessions(accessions)
print(accesion['pairing'])
print(accesion['spots'])
