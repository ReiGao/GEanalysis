import os
import sys
import datetime
a=sys.argv[1]+"/"
listname="list.txt"
listfile=open(a+listname).readlines()
Results=open(a+"DesiredandByproductsResults.txt","w")
Results.write("Group"+"\t"+"Filename"+"\t"+"Mutated"+"\t"+"ReferenceSequence"+"\t"+"All"+"\t"+"Desired"+"\t"+"Wildtype"+"\t"+"Byproducts"+"\n")


for desired in listfile:
    #print(desired)
    filelist=os.listdir(a)
    for i in sorted(filelist):
        if i.endswith(".detil.csv") and i.find(desired.split()[0])!=-1:
            snpfile=open(a+i).readlines()
            indelfile=open(a+i.split(".")[0]+".merge.Indeldetil.csv").readlines()
            print("\n",i,"...........")
            print("desired edit is"+" "+str(desired.split()[1].strip()))
            print("snpfile is"+" "+str(i))
            print("indelfile is"+" "+str(a+i.split(".")[0]+".merge.Indeldetil.csv"))
            #print(type(indelfile))
            #print(indelfile[1])
            #print(a+i.split(".")[0]+".merge.Indeldetil.csv")          
            allreadsnumber=0
            desirednumber=0
            Wildtypenumber=0
            for j in snpfile:
                #print(j.split(",")[0])
                if j.split(",")[0] !="Key" and (j.split(",")[0].find("N")==-1):            
                    allreadsnumber=allreadsnumber+int(j.split(",")[1])
                    #print(j.split(",")[1])
                    if str(j.split(",")[0].strip()).find(str(desired.split()[1].strip()))!=-1:
                        desirednumber=desirednumber+int(j.split(",")[1])
                    if str(j.split(",")[0].strip()).find(str(desired.split()[2].strip()))!=-1:
                        Wildtypenumber=Wildtypenumber+int(j.split(",")[1])
            for k in indelfile:
                if k.split(",")[0] !="Key" and (k.split(",")[0].find("N")==-1):            
                    allreadsnumber=allreadsnumber+int(k.split(",")[2])
                    print("founding desired edits of "+ str(desired.split()[1].strip())+" in indels files.....")
                    if str(k.split(",")[0].strip()).find(str(desired.split()[1].strip()))!=-1:
                        desirednumber=desirednumber+int(k.split(",")[2])
                        #print("Desirededits founded")
                    if str(k.split(",")[0].strip()).find(str(desired.split()[2].strip()))!=-1:
                        Wildtypenumber=Wildtypenumber+int(k.split(",")[2])
            print(allreadsnumber)
            print(desirednumber)
            print(Wildtypenumber)
            Byproductsnumber=str((int(allreadsnumber)-int(desirednumber)-int(Wildtypenumber)))
            Results.write(str(desired.split()[0])+"\t"+str(i)+"\t"+str(desired.split()[1].strip())+"\t"+str(desired.split()[2].strip())+"\t"+str(allreadsnumber)+"\t"+str(desirednumber)+"\t"+str(Wildtypenumber)+"\t"+str(Byproductsnumber)+"\n")
Results.close()     
print("Congratulation! Desired mutation and byproducts searching succeed~\nShuai Jin")
print(datetime.date.today())

    