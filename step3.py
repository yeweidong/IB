#!/usr/bin/env python3
import sys
value = 0.1
with open("analysis.out") as f,open("final.out","w") as outfile,open("uniqgene.txt","w") as outfile2,open("similarity_rate_sort.txt","w") as outfile3:
    yesno = {}
    diff = {}
    uniq_gene = set()
    outfile.write("bone\tGene_ID\tAccession\tQuery_cover\tSubject_cover\tIdentity\tscore\tsimilarity_rate\n")
    outfile3.write("gene\tsimilarity_rate_diff\n")
    for line in f:
        data = line.strip().split("\t")
        if data[0] != "bone":
            if data[1].split("_")[0] not in yesno:
                yesno[data[1].split("_")[0]] = {}
                if data[0] not in yesno[data[1].split("_")[0]]:
                    yesno[data[1].split("_")[0]][data[0]] = {}
                    yesno[data[1].split("_")[0]][data[0]][line] = float(data[3])*float(data[4])*float(data[5])*0.000001
                else:
                    yesno[data[1].split("_")[0]][data[0]][line] = float(data[3])*float(data[4])*float(data[5])*0.000001
            else:
                if data[0] not in yesno[data[1].split("_")[0]]:
                    yesno[data[1].split("_")[0]][data[0]] = {}
                    yesno[data[1].split("_")[0]][data[0]][line] = float(data[3])*float(data[4])*float(data[5])*0.000001
                else:
                    yesno[data[1].split("_")[0]][data[0]][line] = float(data[3])*float(data[4])*float(data[5])*0.000001
    for gene in yesno.keys():
        if len(yesno[gene]) == 2:
            similarity_rate_diff = min(yesno[gene]["yes"].values()) - max(yesno[gene]["no"].values())
            if similarity_rate_diff > value:
                diff[gene] = similarity_rate_diff
                uniq_gene.add(gene)
                outfile2.write(gene + "\n")
                for yn in yesno[gene].keys():
                    for fish in yesno[gene][yn].keys():
                        outfile.write(fish.strip() + "\t" + str(yesno[gene][yn][fish]) + "\n")
    all_gene = set()
    print(name,len(uniq_gene))
    for value in sorted(list(diff.values()), reverse=True):
        for gene in diff.keys():
            if diff[gene] == value and gene not in all_gene:
                outfile3.write(gene + "\t" + str(value) + "\n")
                all_gene.add(gene)
