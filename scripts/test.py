import urllib
import re

geneProd = urllib.urlopen("http://regulondb.ccg.unam.mx/menu/download/datasets/files/GeneProductSet.txt")
growth = urllib.urlopen("http://regulondb.ccg.unam.mx/menu/download/datasets/files/GCSet.txt")

bdd = {
    "ECK120002020" : {
        "gene" : None,
        "product" : None,
        "growth" : []
    }
}
for line in geneProd:
    if(re.match("^#",line)):
        continue
    line = line.strip()
    lineSplit = line.split("\t")
    for loci in bdd:
        if(loci == lineSplit[0]):
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
            conditionGrowth["condition"] = lineSplit[1]
            conditionGrowth["effect"] = lineSplit[6]
            bdd[loci]["growth"].append(conditionGrowth)
print(bdd)