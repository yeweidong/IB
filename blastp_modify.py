#!/usr/bin/env python3
import sys

def arg(argv):
    import argparse
    parser = argparse.ArgumentParser(description = "Recalculated the blastp results file, including Cover Percent and Identity",epilog = "")
    parser.add_argument("-f",help="The file of protein_compere_bone.out")
    parser.add_argument("-o",help="outfile name")
    args = parser.parse_args()
    blastp = args.f
    out = args.o
    with open(blastp) as f1,open(out,"w") as outfile:
        outfile.write("Gene_ID\tAccession\tQuery_cover\tSubject_cover\tIdentity\tscore\n")
        all_data = {}         
        clean_blastp_out(f1,all_data)
        statistics_cover_identity(all_data,outfile)

def clean_blastp_out(f1,all_data):
    for line in f1:
        data = line.strip().split("\t")
        align = int(data[5])
        mac = int(data[7])
        qlen = int(data[2])
        qstart = int(data[9])
        qend = int(data[10])
        slen = int(data[3])
        sstart = int(data[11])
        send = int(data[12])
        score = float(data[14])
        if data[0] not in all_data:
            all_data[data[0]] = {}
            all_data[data[0]][data[1]] = []
            all_data[data[0]][data[1]].append((qlen,slen,align,mac,abs(qend-qstart)+1,abs(send-sstart)+1,score,qstart,qend))
        else:
            if data[1] not in all_data[data[0]]:
                all_data[data[0]][data[1]] = []
                all_data[data[0]][data[1]].append((qlen,slen,align,mac,abs(qend-qstart)+1,abs(send-sstart)+1,score,qstart,qend))
            else:
                n = 0
                for i in range(len(all_data[data[0]][data[1]])):
                    start = all_data[data[0]][data[1]][i][-2]
                    end = all_data[data[0]][data[1]][i][-1]
                    if qstart < start and qend < start or qstart > end and qend > end:
                        n += 1
                if n == len(all_data[data[0]][data[1]]):
                    all_data[data[0]][data[1]].append((qlen,slen,align,mac,abs(qend-qstart)+1,abs(send-sstart)+1,score,qstart,qend))
                        

def statistics_cover_identity(all_data,outfile):
    for key1 in all_data.keys():
        for key2 in all_data[key1].keys():
            if len(all_data[key1][key2]) == 1:
                d = all_data[key1][key2][0]
                outfile.write(key1 + "\t" + key2 + "\t" + str(round(d[4]*100/d[0],2)) + "\t" + str(round(d[5]*100/d[1],2)) + "\t" + str(round(d[3]*100/d[2],2)) + "\t" + str(d[6]) + "\n")
            else:
                d = all_data[key1][key2]
                qlength = d[0][0]
                slength = d[0][1]
                alignlen = 0
                match = 0
                qover = 0
                sover = 0
                score = 0
                n = 0
                for i in all_data[key1][key2]:
                    alignlen += d[n][2]
                    match += d[n][3]
                    qover += d[n][4]
                    sover += d[n][5]
                    score += d[n][6]
                    n += 1
                if sover >= slength:
                    sover = slength
                outfile.write(key1 + "\t" + key2 + "\t" + str(round(qover*100/qlength,2)) + "\t" + str(round(sover*100/slength,2)) + "\t" + str(round(match*100/alignlen,2)) + "\t" + str(score) + "\n")

if __name__ == "__main__":
    arg(sys.argv)

