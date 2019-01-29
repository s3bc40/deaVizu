import urllib
import re


Locus = open("Liste_Locus.txt","r").read()
Locus = Locus.split("\n")


Listt = []
for i in range (len(Locus)):
	Listt.append(Locus[i])
	

geneProd = urllib.urlopen("http://regulondb.ccg.unam.mx/menu/download/datasets/files/GeneProductSet.txt")
growth = urllib.urlopen("http://regulondb.ccg.unam.mx/menu/download/datasets/files/GCSet.txt")

bdd = {}
classes = []
fonctions = ["Growth on non-optimal carbon source", "growth on no-optimal carbon source", "Glucose starvation", "growth on glucose", "growth with glucose", "supplemented with glucose", "Stationary phase", "Culture in late exponential growth phase", "exponential growth"]
for locus in range(len(Listt)):
	loci = {
	
		"gene" : None,
		"product" : None,
		"growth" : []
	    }
	bdd[Listt[locus]]= loci
for line in geneProd:
    if(re.match("^#",line)):
	continue
    line = line.strip()
    lineSplit = line.split("\t")
    for loci in bdd:
	if(loci == lineSplit[0] and len(lineSplit) >= 7):
	    bdd[loci]["gene"] = lineSplit[1]
	    bdd[loci]["product"] = lineSplit[6]
for line in growth:
    if(re.match("^#",line)):
	continue
    line = line.strip()
    lineSplit = line.split("\t")
    for loci in bdd:
	if(bdd[loci]["gene"] == lineSplit[5]):
		conditionGrowth = {
	        "condition" : None,
	        "effect" : None
	        }
		if lineSplit[1] in fonctions:        
			conditionGrowth["condition"] = lineSplit[1]
			conditionGrowth["effect"] = lineSplit[6]
			bdd[loci]["growth"].append(conditionGrowth)
		#else:
			#del bdd[loci]
print(bdd)
print (len(bdd))





	
