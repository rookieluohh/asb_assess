#!/data/luohaohui/Biosoft/miniconda3/envs/myenv/bin/python
# -*- coding: UTF-8 -*-

import argparse

def ArgParse():
    group = argparse.ArgumentParser(description='A python script for genome assessment.')
    group.add_argument('-i','--input',help='assemble result with fasta format.',required=True)
    group.add_argument('-r','--reference',help='reference sequence with fasta format.',required=True)
    group.add_argument('-k','--kmer-length',type=int,help='the kmer length used in assessment, default=21',default=21)
    group.add_argument('-o','--out-prefix',help='prefix of transforming assemble result and reference sequence into fasta format',required=True)

    return group.parse_args()

def integrateReadLine(fa,prefix):
    fi = open(fa,"r")
    lines = fi.readlines()
    fi.close()
    dic1 = {}
    base_line = ""
    read_name = ""
    for line in lines:
        data = line.strip()
        if data != "":
            if data[0] == ">":
                a = 1
                read_name = data
                base_line = ""
            else:
                a = 0
            if a == 0:
                base_line += data
                dic1[read_name] = base_line
    fa_list = list(dic1.items())
    with open(prefix + "_" + fa,"w") as fo:
        for i in range(len(fa_list)):
            if i != len(fa_list)-1:
                fo.write("{}\n{}\n".format(fa_list[i][0],fa_list[i][1]))
            else:
                fo.write("{}\n{}".format(fa_list[i][0],fa_list[i][1]))
    fo.close()

def refUniqueKmerSearch(ref,kmer_length):
    fi = open(ref,"r")
    lines = fi.readlines()
    fi.close()
    external_dic = {}
    internal_dic = {}
    chr_name = ""
    for line in lines:
        line = line.strip().split("\t")
        if line != "":
            if line[0][0] == ">":
                chr_name = line[0][1:]
                internal_dic = {}
            else:
                for i in range(len(line[0])-kmer_length+1):
                    internal_dic[line[0][i:i+kmer_length]] = internal_dic.get(line[0][i:i+kmer_length],0) + 1
                external_dic[chr_name] = internal_dic
    sumKmer_dic = {}
    for line in lines:
        line = line.strip()
        if line[0] != ">":
            for i in range(len(line)-kmer_length+1):
                sumKmer_dic[line[i:i+kmer_length]] = sumKmer_dic.get(line[i:i+kmer_length],0) + 1
    unique_kmer = []
    for key,value in sumKmer_dic.items():
        if value == 1:
            unique_kmer.append(key)
    for K,V in external_dic.items():
        dic = {}
        for key,value in V.items():
            if key in unique_kmer:
                dic[key] = value
        external_dic[K] = dic
    return [external_dic,unique_kmer]  #external_dic={"chr":{"ATCG":5}}   unique_kmer=["ATCG",...,"ATCG"]

def refUniqueKmer_in_contig(ref,contig,kmer_length):
    ref_unikmer_list = refUniqueKmerSearch(ref,kmer_length)
    ref_unikmer = ref_unikmer_list[1]
    contig_data = open(contig,"r")
    contig_lines = contig_data.readlines()
    contig_data.close()
    external_dic = {}
    internal_dic = {}
    contig_unikmer = []
    sumSameKmerCount_dic = {}
    contig_name = ""
    for line in contig_lines:
        line = line.strip().split("\t")[0]
        if line[0] == ">":
            contig_name = line[1:]
            internal_dic = {}
        else:
            for i in range(len(line) - kmer_length + 1):
                for kmer in ref_unikmer:
                    if line[i:i+kmer_length] == kmer:
                        internal_dic[kmer] = internal_dic.get(kmer,0) + 1
                        sumSameKmerCount_dic[kmer] = sumSameKmerCount_dic.get(kmer,0) + 1
            external_dic[contig_name] = internal_dic
    for line in contig_lines:
        line = line.strip().split("\t")[0]
        if line[0] != ">":
            for i in range(len(line) - kmer_length + 1):
                if line[i:i+kmer_length] in ref_unikmer:
                    if line[i:i+kmer_length] not in contig_unikmer:
                        contig_unikmer.append(line[i:i+kmer_length])
    return [external_dic,sumSameKmerCount_dic,contig_unikmer]   #external_dic={"ctg0001":{"ATCG":5}}  sumSameKmerCount_dic={"ATCG":2,...,"ATCG":5}   contig_unikmer=["ATCG",...,"ATCG"]

def getRefUniKmerPos(ref,kmer_length):
    eachchr_unikmer = refUniqueKmerSearch(ref,kmer_length)[0]
    data = open(ref,"r")
    lines = data.readlines()
    data.close()
    kmerpos_dic = {}
    chr_name = ""
    for line in lines:
        line = line.strip().split("\t")[0]
        if line[0] == ">":
            chr_name = line[1:]
            kmerpos_dic[chr_name] = []
    for line in lines:
        line = line.strip().split("\t")[0]
        if line[0] == ">":
            chr_name = line[1:]
        else:
            for i in range(len(line)-kmer_length+1):
                if line[i:i+kmer_length] in eachchr_unikmer[chr_name]:
                    kmerpos_dic[chr_name].append((i+1,i+kmer_length))
    return kmerpos_dic          #kmerpos_dic={"chr1":[(1,21),...,(41,51)]}

def getCtgUniKmerPos(ref,contig,kmer_length):
    ctgUniKmer_list = refUniqueKmer_in_contig(ref,contig,kmer_length)[2]
    data = open(contig,"r")
    lines = data.readlines()
    data.close()
    kmerpos_dic = {}
    contig_name = ""
    for line in lines:
        line = line.strip().split("\t")[0]
        if line[0] == ">":
            contig_name = line[1:]
            kmerpos_dic[contig_name] = []
    for line in lines:
        line = line.strip().split("\t")[0]
        if line[0] == ">":
            contig_name = line[1:]
        else:
            for i in range(len(line)-kmer_length+1):
                if line[i:i+kmer_length] in ctgUniKmer_list:
                    kmerpos_dic[contig_name].append((i+1,i+kmer_length))
    return kmerpos_dic       #kmerpos_dic={"ctg0001":[(1,21),...,(41,61)]}

def removeOverlapKmer(classed_kmer_dic):
    removeOL_kmerpos_dic = {}
    for key,value in classed_kmer_dic.items():
        my_list1 = value[:]
        while True:
            my_list2 = my_list1[:]
            my_list1 = []
            if len(my_list2) != 1:
                for i in range(len(my_list2)-1):
                    if my_list2[i][1] >= my_list2[i+1][0] and my_list2[i][1] != my_list2[i+1][1]:
                        my_list1.append((my_list2[i][0],my_list2[i+1][1]))
                    else:
                        my_list1.append(my_list2[i])
                        if i == len(my_list2)-2:
                            my_list1.append(my_list2[i+1])
            else:
                my_list1 = my_list2[:]
            if my_list1 == my_list2:
                removeOL_kmerpos_dic[key] = my_list1
                break
    removeOL_kmerpos = {}
    for k in removeOL_kmerpos_dic.keys():
        removeOL_kmerpos[k] = []
    for key,value in removeOL_kmerpos_dic.items():
        list = range(value[0][0]+1,value[0][1]+1)
        for i in range(len(value)):
            if (value[i][0] >= list[0] and value[i][1] <= list[-1]):
                continue
            elif (value[i][0] <= list[0] and value[i][1] >= list[-1]):
                list = range(value[i][0],value[i][1]+1)
            else:
                removeOL_kmerpos[key].append((list[0],list[-1]))
                list = range(value[i][0],value[i][1]+1)
        removeOL_kmerpos[key].append((list[0], list[-1]))
    return removeOL_kmerpos     #removeOL_kmerpos_dic={"chr/ctg":[(1,21),...,(31,51)]}

def getHeaderKmer(fasta,headerkmerpos_dic,kmer_length):
    data = open(fasta,"r")
    lines = data.readlines()
    data.close()
    kmer_pos_dic = {}
    seqname_kmer_dic = {}
    sequence_name = ""
    for line in lines:
        line = line.strip().split("\t")[0]
        if line[0] == ">":
            sequence_name = line[1:]
            kmer_pos_dic = {}
        else:
            for key,value in headerkmerpos_dic.items():
                if sequence_name == key:
                    for pos in value:
                        kmer_pos_dic[line[pos[0] - 1:pos[0] + kmer_length + 1]] = []
            seqname_kmer_dic[sequence_name] = kmer_pos_dic
    for line in lines:
        line = line.strip().split("\t")[0]
        if line[0] == ">":
            sequence_name = line[1:]
        else:
            for key,value in headerkmerpos_dic.items():
                if sequence_name == key:
                    for pos in value:
                        seqname_kmer_dic[sequence_name][line[pos[0]-1 : pos[0]+kmer_length+1]].append(pos)
    return seqname_kmer_dic    #seqname_kmer_dic={"chr/ctg":{"ATCG":[(1,21),...,(31,51)]}}

def assembleAssessment(ref_fa,asb_fa,kmer_length,prefix):
    integrateReadLine(ref_fa,prefix)
    integrateReadLine(asb_fa,prefix)
    refUniKmerPos = getRefUniKmerPos(prefix + "_" + ref_fa,kmer_length)
    asbUniKmerPos = getCtgUniKmerPos(prefix + "_" + asb_fa,prefix + "_" + asb_fa,kmer_length)
    refHeaderKmerPos = removeOverlapKmer(refUniKmerPos)
    refUniKmerPos = None
    asbHeaderKmerPos = removeOverlapKmer(asbUniKmerPos)
    asbUniKmerPos = None
    refHeaderKmer = getHeaderKmer(prefix + "_" + ref_fa,refHeaderKmerPos,kmer_length)
    refHeaderKmerPos = None
    asbHeaderKmer = getHeaderKmer(prefix + "_" + asb_fa,asbHeaderKmerPos,kmer_length)
    asbHeaderKmerPos = None
    asbKmerCount = {}
    for Key,Value in asbHeaderKmer.items():
        for k,v in Value.items():
            asbKmerCount[k] = asbKmerCount.get(k,0) + len(v)
    singleCopyKmerNum = 0
    duplicateKmerNum = 0
    colKmerNum = 0
    for key,value in refHeaderKmer.items():
        for k,v in value.items():
            colKmerNum += 1
    for key,value in asbKmerCount.items():
        if value == 1:
            singleCopyKmerNum += 1
        else:
            duplicateKmerNum += 1
    singleCopy = singleCopyKmerNum / colKmerNum
    duplicateRate = duplicateKmerNum / colKmerNum
    ctgkmer_dic = {}
    for key,value in asbHeaderKmer.items():
        ctgkmer_dic[key] = []
    for key,value in asbHeaderKmer.items():
        for k,v in value.items():
            if k not in ctgkmer_dic[key]:
                ctgkmer_dic[key].append(k)
    chrkmer_dic = {}
    for key,value in refHeaderKmer.items():
        chrkmer_dic[key] = []
    for key,value in refHeaderKmer.items():
        for k,v in value.items():
            if k not in chrkmer_dic[key]:
                chrkmer_dic[key].append(k)
    chr_to_kmerlist = {}
    contig_to_chr = {}
    for ctg, ctgkmer in asbHeaderKmer.items():
        chr_to_kmerlist = {}
        for chr, chrkmer in refHeaderKmer.items():
            chr_to_kmerlist[chr] = []
        contig_to_chr[ctg] = chr_to_kmerlist
    for asb_key, asb_value in ctgkmer_dic.items():
        for i in asb_value:
            for ref_key, ref_value in chrkmer_dic.items():
                if i in ref_value and i not in contig_to_chr[asb_key][ref_key]:
                    contig_to_chr[asb_key][ref_key].append(i)        #contig_to_chr={"ctg0001":{"chr01":["ATCG",...,"ATCG"]}}chr和ctg共有kmer
    ctgkmer_dic = None
    chrkmer_dic = None
    proportion_of_the_largest_categories = 0
    largest_categories_ex = {}
    largest_categories_in = {}
    for key,value in contig_to_chr.items():               #找出最大类
        largest_categories_num = 0
        largest_categories_in = {}
        for k,v in value.items():
            if len(v) > largest_categories_num:
                largest_categories_num = len(v)
        for k,v in value.items():
            if len(v) == largest_categories_num:
                largest_categories_in[k] = v
                largest_categories_ex[key] = largest_categories_in          #largest_categories_ex={"ctg0001":{"chr01":["ATCG",...,"ATCG"]}}
        proportion_of_the_largest_categories += float(largest_categories_num) / float(colKmerNum)
    aveEachCtgDistance_sum = 0
    contig_to_chr = None
    for key,value in largest_categories_ex.items():
        ref_base_pos_dic = {}
        asb_base_pos_dic = {}
        t = 0
        eachCtgDistance_sum = 0
        for k,v in value.items():
            for i in range(len(v)):
                ref_base_pos_dic[v[i]] = refHeaderKmer[k][v[i]]        #ref_base_pos_dic={"ATCG":[(1,21),...,(31,51)]}
                asb_base_pos_dic[v[i]] = asbHeaderKmer[key][v[i]]
            ref_base_list = list(ref_base_pos_dic.items())
            ref_base_list.sort(key = lambda x:x[1][0][0], reverse=False)
            for i in range(len(ref_base_list)-1):
                refKmerDistance = abs(ref_base_list[i+1][1][0][0] - ref_base_list[i][1][0][0])
                asbKmerDistance = []
                for j in range(len(asb_base_pos_dic[ref_base_list[i][0]])):
                    for m in range(len(asb_base_pos_dic[ref_base_list[i+1][0]])):
                        distance = abs(asb_base_pos_dic[ref_base_list[i][0]][j][0]-asb_base_pos_dic[ref_base_list[i+1][0]][m][0])
                        asbKmerDistance.append(distance)
                for j in asbKmerDistance:
                    eachCtgDistance_sum += abs(j-refKmerDistance)
                    t += 1
        aveEachCtgDistance_sum += eachCtgDistance_sum / t
    contig_nums = len(largest_categories_ex)
    ave_distance_diff = aveEachCtgDistance_sum / contig_nums
    print("{:<12.3f}Complete\n"
          "{:<12.3f}Complete and single-copy\n"
          "{:<12.3f}Complete and duplicated\n"
          "{:<12.3f}Proportion of the largest categories\n"
          "{:<12.1f}ave distance diff\n".format(singleCopy+duplicateRate,singleCopy,duplicateRate,proportion_of_the_largest_categories,ave_distance_diff))
    with open(prefix + ".result","w") as fi:
        fi.write("{:<12.3f}Complete\n"
          "{:<12.3f}Complete and single-copy\n"
          "{:<12.3f}Complete and duplicated\n"
          "{:<12.3f}Proportion of the largest categories\n"
          "{:<12.1f}ave distance diff\n".format(singleCopy+duplicateRate,singleCopy,duplicateRate,proportion_of_the_largest_categories,ave_distance_diff))
        fi.close()

'''if __name__ == "__main__":
    opt = ArgParse()
    ref_seq = opt.input
    asb_seq = opt.reference
    kmer_len = opt.kmer_length
    prefix = opt.out_prefix
    assembleAssessment(ref_seq,asb_seq,kmer_len,prefix)
'''
assembleAssessment("ref_test.txt","asb_test.txt",21,"test")






