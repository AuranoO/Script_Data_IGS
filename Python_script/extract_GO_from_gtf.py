with open("../ref_kamtchatka/A_castellanii_HiC-Seq2.gtf", "r") as gtf:
    gtf_read = gtf.readlines()

# go_to_append = open("../final_data/gene2GO_amibe_hic2.map", "w")
for line in gtf_read:
    if "Ontology_term" in line:
        gene_id = line.split("gene_id ")[1].split('"')
        # go_to_append.write(gene_id[1]+"\t")
        # print(goterm)
        goterm = line.split("Ontology_term")[1].split(";")
        # go_to_append.write(goterm[0]+"\t")
        goterm = line.split("\t")[8].split(";")[len(goterm)]
        start = ": "
        end = " -"
        counter = goterm.count("-")
        if counter == 1:
            go = goterm[goterm.find(start) + len(start):goterm.rfind(end)]
            # print(go)
            # go_to_append.write(go+"\n")
        else:
            goterm = goterm.split("GO_")
            del goterm[0]
            cat_go = ""
            for i in range(0, len(goterm)):
                go = goterm[i].split(": ")[1].split(" -")[0]
                cat_go += go + " "
            # go_to_append.write(cat_go+"\n")
    else:
        pass

# go_to_append.close()
