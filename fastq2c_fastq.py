import sys
import os
import gzip
input_neme=sys.argv[1]
if input_neme.strip().split(".")[-1]=="gz": 
   input_fastqc=gzip.GzipFile(input_neme,"r")
else:
   input_fastqc=open(input_neme,"r")

output_fastq=open(sys.argv[2]+".fastq","w")
huan=[]
#for line in input_fastqc.xreadlines():
#with input_fastqc as pf:
#  for line in pf:
for line in input_fastqc:
    id_t=line[0]
    line=line.strip().split("\n")[0]
    if id_t=="@":
       if len(huan)==4:
          output_fastq.write(str(huan[0])+"\n"+str(huan[1])+"\n"+str(huan[2])+"\n"+str(huan[3])+"\n")
       elif len(huan)>0 and len(huan)<4:
          print huan[0] 
       huan=[]
       huan.append(line)
    elif id_t=="A" or id_t=="T" or id_t=="C" or id_t=="G":
       line=line.strip().split("@")
       if len(line)>=2:
          huan=[]
          huan.append("@"+line[-1])
       else:
          if line[0][-1]=="A" or line[0][-1]=="T"  or  line[0][-1]=="C" or line[0][-1]=="G": 
             huan.append(line[0])
    else:
       huan.append(line)
if len(huan)==4:
    output_fastq.write(str(huan[0])+"\n"+str(huan[1])+"\n"+str(huan[2])+"\n"+str(huan[3])+"\n")
elif len(huan)>0 and len(huan)<4:
    print huan[0]
input_fastqc.close()
output_fastq.close()
