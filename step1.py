#!/usr/bin/env python3
import os,sys

def arg(argv):
    import argparse
    parser = argparse.ArgumentParser(description = "This script is used to compere different expression gene protein to yes_bone and no_bone fish protein",epilog = "")
    parser.add_argument("-f",help="The file of protein.fa")
    args = parser.parse_args()
    fa = args.f
    blastp_bone(fa,"no_bone_gene.protein.fa","protein_compere_no_bone.out")#Compared with fish species without intermuscular bones
    blastp_bone(fa,"yes_bone_gene.protein.fa","protein_compere_yes_bone.out")#Compared with fish species with intermuscular bones
    score_rank("clean_protein_compere_yes_bone.out","clean_protein_compere_no_bone.out")

def blastp_bone(fa,inputf,outputf):

    os.system("blastp -query %s -db %s -out %s -evalue 0.0001 -matrix BLOSUM62 -outfmt '6 qseqid sseqid qlen slen pident length mismatch nident gapopen qstart qend sstart send evalue bitscore' -num_threads 16 " % (fa,inputf,outputf))
    os.system("blastp_modify.py -f %s -o clean_%s" % (outputf,outputf))

def score_rank(F1,F2):

    def get_score(f,outfile,ys):
        for line in f:
            if line[0:4] != "Gene":
                gene = line.split("\t")[0]
                fish_name = line.split("\t")[1].split("-")[-1]
                Query_cover = float(line.split("\t")[2])
                Subject_cover = float(line.split("\t")[3])
                Identity = float(line.strip().split("\t")[4])
                score = float(line.strip().split("\t")[5])
                if gene not in gene_score:
                    gene_score[gene] = {}
                    gene_score[gene][score] = ys + "\t" + line
                    gene_score[gene]['fish'] = set()
                    gene_score[gene]['fish'].add(fish_name)
                else: 
                    if score not in gene_score[gene] and fish_name not in gene_score[gene]['fish']:
                        gene_score[gene][score] = ys + "\t" + line
                        gene_score[gene]['fish'].add(fish_name)
                    elif score in gene_score[gene] and fish_name not in gene_score[gene]['fish']:
                        gene_score[gene][score - 0.01] = ys + "\t" + line
                        gene_score[gene]['fish'].add(fish_name)
    def rank():
        for gene in list(gene_score.keys()):
            del gene_score[gene]['fish']
            for score in sorted(gene_score[gene],reverse=True):
                outfile.write(gene_score[gene][score])


    with open(F1) as f1,open(F2) as f2,open("score_rank.out","w") as outfile:
        gene_score = {}
        outfile.write("bone\tGene_ID\tAccession\tQuery_cover\tSubject_cover\tIdentity\tscore\n")
        get_score(f1,outfile,'yes')
        get_score(f2,outfile,'no')
        rank()

if __name__ == "__main__":
    arg(sys.argv)

