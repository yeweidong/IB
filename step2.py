#!/usr/bin/env python3
with open("score_rank.out") as f,open("analysis.out","w") as outfile:
    data = {}
    for line in f:
        if line[0] == "b":
            outfile.write(line)
        else:
            lists = line.strip().split("\t")
            bone = lists[0]
            gene = lists[1].split("_")[0]
            ident = lists[5]
            score = float(lists[6])
            fish = lists[2].split("-")[1]
            if gene not in data:
                data[gene] = {}
                if bone == 'yes':
                    data[gene]['yes'] = {}
                    data[gene]['yes'][fish] = {}
                    data[gene]['yes'][fish][score] = ident + "#" + line
                else:
                    data[gene]['no'] = {}
                    data[gene]['no'][fish] = {}
                    data[gene]['no'][fish][score] = ident + "#" + line
            else:
                if bone == 'yes':
                    if 'yes' not in data[gene]:
                        data[gene]['yes'] = {}
                        data[gene]['yes'][fish] = {}
                        data[gene]['yes'][fish][score] = ident + "#" + line
                    else:
                        if fish not in data[gene]['yes']:
                            data[gene]['yes'][fish] = {}
                            data[gene]['yes'][fish][score] = ident + "#" + line
                        else:
                            data[gene]['yes'][fish][score] = ident + "#" + line
                else:
                    if 'no' not in data[gene]:
                        data[gene]['no'] = {}
                        data[gene]['no'][fish] = {}
                        data[gene]['no'][fish][score] = ident + "#" + line
                    else:
                        if fish not in data[gene]['no']:
                            data[gene]['no'][fish] = {}
                            data[gene]['no'][fish][score] = ident + "#" + line
                        else:
                            data[gene]['no'][fish][score] = ident + "#" + line
    test = open("gene2gene.txt","w")#
    for gene in data.keys():
        for yn in data[gene].keys():
            for fish in list(data[gene][yn].keys()):
                max_fish_score = max(data[gene][yn][fish])
                #Obtain the sequence with the highest mapping score to a gene in zebrafish
                data[gene][yn][max_fish_score] = data[gene][yn][fish][max_fish_score]
                #test.write(gene + "\t" + data[gene][yn][fish][max_fish_score])
                test.write(gene + "\t" + data[gene][yn][fish][max_fish_score].split("\t")[2].replace("-","\t").replace("[","").replace("]","") + "\n")
                del data[gene][yn][fish]
    #remove duplicate
    sort_data = {}
    for gene in data.keys():
        if 'no' in data[gene] and 'yes' in data[gene] and len(data[gene]['yes']) == 3:
            if min(data[gene]['yes']) > max(data[gene]['no']):
                diff = float(data[gene]['yes'][min(data[gene]['yes'])].split("#")[0]) - float(data[gene]['no'][max(data[gene]['no'])].split("#")[0])
                ident = float(data[gene]['yes'][min(data[gene]['yes'])].split("#")[0])
                if diff > 10 and ident > 90:
                    #Here are the screening requirements
                    if diff not in sort_data:
                        sort_data[diff] = data[gene]
                    else:
                        sort_data[diff + 0.01] = data[gene]
        elif 'no' not in data[gene] and len(data[gene]['yes']) == 3 and float(data[gene]['yes'][min(data[gene]['yes'])].split("#")[0]) > 50:
            for score in data[gene]['yes'].keys():
                outfile.write(data[gene]['yes'][score].split("#")[1])
        elif 'yes' not in data[gene]:
            pass
            #print(gene)
    #sort the results
    for diff in sorted(sort_data,reverse = True):
        for yn in sort_data[diff].keys():
            for score in sort_data[diff][yn].keys():
                outfile.write(sort_data[diff][yn][score].split("#")[1])


                
