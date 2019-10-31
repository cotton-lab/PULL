#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Pull.py 1.0.0 -- Plant full length lncRNA identity

Usage: Pull.py [options]  

Options: 
    -h                               Show this screen.
    -v                               Show version.
    -g GENO                          Genome fasta file    
    -r REF                           Genome annotation file in gtf format.
    -c CAGE                          A directory contain cage sample.ctss files
    -p PAS                           A directory contain pas.ctss files
    --ngs NGS_BED                    A bed12 format file assembled by NGS sequencing data, you should assembel all sample to a final merged.bed 
    --ngs_bam NGS_BAM                A directory contain all bam files for different samples.
    --pb PACBIO                      A bed12 format file assembled by Pacbio sequencing data.
"""


__author__ = 'Yanjun Chen (chenyj@whu.edu.cn)'
__version__ = '1.0.0'

import sys
import os
from docopt import docopt
import glob
from math import sqrt
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from collections import OrderedDict

def multipl(a,b):
    sumofab=0.0
    for i in range(len(a)):
        temp=a[i]*b[i]
        sumofab+=temp
    return sumofab

def corrcoef(x,y):
    n=len(x)
    sum1=sum(x)
    sum2=sum(y)
    sumofxy=multipl(x,y)
    sumofx2 = sum([pow(i,2) for i in x])
    sum2=sum(y)
    sumofxy=multipl(x,y)
    sumofx2 = sum([pow(i,2) for i in x])
    sumofy2 = sum([pow(j,2) for j in y])
    num=sumofxy-(float(sum1)*float(sum2)/n)
    den=sqrt((sumofx2-float(sum1**2)/n)*(sumofy2-float(sum2**2)/n))
    return num/den

def get_BSgenome(geno_seq):
    return_code = os.system("mkdir -p tmp/01-build_X_BSgenome/seqs tmp/01-build_X_BSgenome/seed_file") >> 8
    return_code = os.system("faToTwoBit " + geno_seq + " tmp/01-build_X_BSgenome/seqs/X.2bit") >> 8
    x = []
    for r in SeqIO.parse(geno_seq,"fasta"):
        id = str(r.id)
        x.append("\""+id+"\"")
        y ="c(" +  ','.join(x) + ")"
        pwd_dir = str(os.path.abspath('.'))
        w = open('tmp/01-build_X_BSgenome/seed_file/BSgenome.X-seed','w')
        w.write("Package: BSgenome.X.ucsc\nTitle: X full genome\nDescription: X full genome as provided by Pacbio and stored in Biostrings objects.\nVersion: 1.1.1\norganism: X\ncommon_name: XX\nprovider: UCSC\nprovider_version: V1\nrelease_date: May 2018\nrelease_name: Gs May 2018 sequence\norganism_biocview: X\nBSgenomeObjname: XX\nseqnames: " +y+ "\nSrcDataFiles: X.2bit\nPkgExamples: genome$"+x[0][1:-1]+"\nseqs_srcdir: " + str(os.path.abspath('.')) + "/tmp/01-build_X_BSgenome/seqs \nseqfile_name: X.2bit" )
        next_dir = str(os.popen("find ~/ -name Py3_pull ").read().strip()) + "/lib/R/library/BSgenome/extdata/GentlemanLab"
        os.chdir(next_dir)
        w = open('tmp.r','w')
        w.write("library(BSgenome)\nforgeBSgenomeDataPkg(\"" + pwd_dir + "/tmp/01-build_X_BSgenome/seed_file/BSgenome.X-seed\") \n")
        w.close()
        return_code = os.system("Rscript tmp.r")  >> 8
        return_code = os.system("R CMD build BSgenome.X.ucsc/")
        return_code = os.system("R CMD check BSgenome.X.ucsc_1.1.1.tar.gz")
        return_code = os.system("R CMD INSTALL BSgenome.X.ucsc_1.1.1.tar.gz")
        return_code = os.system("rm tmp.r")  >> 8
        os.chdir(pwd_dir)

def extract_dominant_ctss(target_dir,x):
    return_code = os.system("mkdir -p tmp/" + x ) >> 8
    return_code = os.system("cat " + target_dir + "*ctss  > tmp/" + x +"/MIX.bed4") >> 8
    pwd_dir = str(os.path.abspath('.'))
    next_dir = pwd_dir + "/tmp/" + x +"/"
    os.chdir(next_dir)
    w= open('ALL.ctss','w')
    ID = OrderedDict()
    with open('MIX.bed4','r') as f:
        for line in f:
            id = '\t'.join(line.strip().split('\t')[0:3])
            count = int(line.strip().split('\t')[3])
            if id not in ID:
                ID[id] = count
            else :
                ID[id] = ID[id] + count

        for k in ID.keys():
            w.write(k + '\t' + str(ID[k]) +'\n')
    w.close()
    w = open('tmp1.r','w')
    w.write("library(CAGEr)\nlibrary(BSgenome.X.ucsc)\ninputFile <- \"ALL.ctss\"\nmyCAGEset <- new(\"CAGEset\", genomeName = \"BSgenome.X.ucsc\", inputFiles = inputFile, inputFilesType = \"ctss\", sampleLabels = \"ALL\")\ngetCTSS(myCAGEset)\nnormalizeTagCount(myCAGEset, method = \"none\")\nclusterCTSS(object = myCAGEset, threshold = 5, thresholdIsTpm = TRUE, nrPassThreshold = 1, method = \"distclu\", maxDist = 20,removeSingletons = TRUE, keepSingletonsAbove = 5)\ntc <- tagClusters(myCAGEset, sample = \"ALL\")\nwrite.table(tc,\"ALL.tc\", sep = \",\", row.names = TRUE)")
    w.close()
    return_code = os.system("Rscript tmp1.r")  >> 8
    w=open('dominant_ctss-summary','w')  
    with open('ALL.tc','r') as f:
        for line in f:
            if 'cluster' not in line:
                s=line.strip().split(',')
                id = s[2][1:6] + '\t' + s[7] +'\t' + s[5][1]+'\t'+s[9] #chr,dominant_ctss_pos,strand,dominant_ctss_counts
                w.write(id+'\n')
    w.close()
    os.chdir(pwd_dir)

def filter_fot_txs(f1):
    filter_file=str(f1)+"-filtered"
    ID_pos = OrderedDict()
    ID_count = OrderedDict()
    w=open(filter_file,'w')
    with open(f1,'r') as f:
        for line in f:
            s = line.strip().split('\t')
            if s[0] not in ID_pos:
                ID_pos[id] =  int(s[1])
                ID_count[id] =  int(s[2])
                w.write(line)
            else :
                pos = float(s[1])
                count = float(s[2])
                if (abs(pos-ID_pos[id]) >= 10) and ((count/ID_count[id]) >= 0.8):
                    w.write(line)
    w.close()

def get_txs_tpm(f2):
    Chr_site = OrderedDict()
    TPM = OrderedDict()
    ID = OrderedDict()
    with open('consensusClusters.bed','r') as f:
        for line in f:
            s=line.strip().split(',')
            cluster_id = s[1]
            chr=s[2][1:6]
            Chr_site[cluster_id] = chr
    w=open(str(f2.split('.')[0])+".tpm",'w')
    with open('consensusClusters.tpm','r') as f:
        for line in f:
            s=line.strip().split(',')
            if 'Anther' in line:
                for i in range(0,len(s)) :
                    s[i] = s[i][1:-1]
                w.write('      '+'\t'.join(s) +'\n')
            else :
                cluster_id = s[0][1:-1]
                count ='\t'.join(s[1:len(s)])
                TPM[cluster_id] = count
    with open(f2,'r') as f:
        for line in f:
            s=line.strip().split('\t')
            cluster_id = s[0]
            site = s[1]
            count = s[2]
            if cluster_id not in ID:
                ID[cluster_id] = site+'\t'+count
            else :
                if (abs(int(site) - int(ID[cluster_id].split('\t')[0])) >= 5) and (int(count) >= 10) :
                    ID[cluster_id] =  ID[cluster_id] + ',' + site+'\t'+count 
    for k in ID.keys():
        x = ID[k].split(',')
        for n in range(0,len(x)) :
            chr = Chr_site[k]
            tss_site = x[n].split('\t')[0]
            count = TPM[k]
            w.write(k+'_'+chr+'_'+tss_site+'\t' +count+'\n')

def get_txs_cluster(target_dir,x):
    pwd_dir = str(os.path.abspath('.'))
    next_dir = pwd_dir + "/tmp/" + x +"/"
    ctss_file_path = str(os.path.abspath((target_dir))) 
    os.chdir(next_dir)
    w = open('tmp2.r','w')
    w.write("library(CAGEr)\nlibrary(BSgenome.X.ucsc)\ninputFiles = list.files(path=\""+ ctss_file_path +"\",full.names = TRUE)\nce<-new(\"CAGEset\",genomeName=\"BSgenome.X.ucsc\",inputFiles=inputFiles,inputFilesType =\"ctss\",sampleLabels= sub(\".ctss\",\"\",basename(inputFiles)))\ngetCTSS(ce)\nnormalizeTagCount(ce,method=\"simpleTpm\")\nclusterCTSS(object=ce,threshold = 1,thresholdIsTpm = TRUE, nrPassThreshold = 1, method = \"distclu\", maxDist = 20, removeSingletons = TRUE, keepSingletonsAbove = 5)\ncumulativeCTSSdistribution(ce, clusters = \"tagClusters\", useMulticore = T,nrCores=20)\nquantilePositions(ce, clusters = \"tagClusters\", qLow = 0.1, qUp = 0.9)\naggregateTagClusters(ce, tpmThreshold = 5, qLow = 0.1, qUp = 0.9, maxDist = 20)\ncl <- consensusClusters(ce)\nwrite.table(cl, file=\"consensusClusters.bed\", sep = \",\", row.names = TRUE)\ntpm <- consensusClustersTpm(ce)\nwrite.table(tpm, file=\"consensusClusters.tpm\", sep = \",\", row.names = TRUE)\n")
    w.close()
    #return_code = os.system("Rscript tmp2.r")  >>8 
    #define vars for different chr dictionary
    with open('dominant_ctss-summary','r') as f:
        for line in f:
            exec ("%s={}"%line.strip().split('\t')[0])
    #assignment for different chr dictionary
    with open('dominant_ctss-summary','r') as f:
        for line in f:
            s =line.strip().split('\t')
            if s[2] == '+' :
                exec ("%s[s[1]]=int(s[3])"%s[0])
    w=open('Plus','w')
    with open('consensusClusters.bed','r') as f:
        for line in f:
            if 'consensus' not in line:
                s =line.strip().split(',')
                if s[5][1] == '+' : #only keep plus strand line
                    cluster_num = s[1]
                    chr = s[2][1:6]
                    st = int(s[3])
                    et = int(s[4])
                    for k in locals()[chr]:
                        if (int(k) >= st) and (int(k) <= et) :
                            w.write(cluster_num +'\t' + k +'\t' +str(locals()[chr][k]) +'\n')
    w.close()
    with open('dominant_ctss-summary','r') as f:
        for line in f:
            exec ("%s={}"%line.strip().split('\t')[0])
    with open('dominant_ctss-summary','r') as f:
        for line in f:
            s =line.strip().split('\t')
            if s[2] == '-' :
                exec ("%s[s[1]]=int(s[3])"%s[0])
    w=open('Minus','w')
    with open('consensusClusters.bed','r') as f:
        for line in f:
            if 'consensus' not in line:
                s =line.strip().split(',')
                if s[5][1] == '-' : #only keep plus strand line
                    cluster_num = s[1]
                    chr = s[2][1:6]
                    st = int(s[3])
                    et = int(s[4])
                    for k in locals()[chr]:
                        if (int(k) >= st) and (int(k) <= et) :
                            w.write(cluster_num +'\t' + k +'\t' +str(locals()[chr][k]) +'\n')
    w.close()
    return_code = os.system("sort -k1,1n -k3,3nr Minus > Minus.sorted")  >>8 
    return_code = os.system("sort -k1,1n -k3,3nr Plus > Plus.sorted")  >>8 
    filter_fot_txs("Minus.sorted")
    filter_fot_txs("Plus.sorted")
    get_txs_tpm("Minus.sorted-filtered")
    get_txs_tpm("Plus.sorted-filtered")
    os.chdir(pwd_dir)

def combin_ngs_and_pacbio(f_ngs,f_pb,ref):
    pwd_dir = str(os.path.abspath('.'))
    f1 = str(os.path.abspath((f_ngs)))
    f2 = str(os.path.abspath((f_pb)))
    return_code = os.system("mkdir -p tmp/Transcript") >> 8
    next_dir = pwd_dir + "/tmp/Transcript"
    os.chdir(next_dir)
    script_dir = str(os.popen("find / -name Pull.py ").read().strip()).replace("Pull.py","script") + "/"
    return_code = os.system(str(script_dir) + "bedToGtf.sh " + f1 + " > NGS-Raw.gtf")
    return_code = os.system(str(script_dir) + "bedToGtf.sh " + f2 + " > Pacbio-Raw.gtf")
    return_code = os.system("echo \"NGS-Raw.gtf\" > mergelist")
    return_code = os.system("echo \"Pacbio-Raw.gtf\" >> mergelist")
    return_code = os.system("stringtie --merge -G " + str(ref) + " -o Raw.gtf mergelist")
    return_code = os.system(str(script_dir)+"gtf2bed Raw.gtf > Raw.bed")    
    os.chdir(pwd_dir)

def only_ngs(f_ngs):
    pwd_dir = str(os.path.abspath('.'))
    return_code = os.system("mkdir -p tmp/Transcript") >> 8
    next_dir = pwd_dir + "/tmp/Transcript"
    return_code = os.system("cp " + f_ngs + " tmp/Transcript/Raw.bed")    

def only_pacbio(f_pb):
    pwd_dir = str(os.path.abspath('.'))
    return_code = os.system("mkdir -p tmp/Transcript") >> 8
    next_dir = pwd_dir + "/tmp/Transcript"
    return_code = os.system("cp " + f_pb + " tmp/Transcript/Raw.bed")

def calculate_fpkm_for_iso(ngs_bam_dir):
    pwd_dir = str(os.path.abspath('.'))
    next_dir = pwd_dir + "/tmp/Transcript"
    os.chdir(next_dir)
    return_code = os.system("echo \"for i in \`ls " + ngs_bam_dir + "/*bam\`; do name=\${i/*\/} ;name=\${name/.bam/}; stringtie -i \${i} -G Raw.gtf -e -p 4 --rf -o ./Gtfiles/\${name}.gtf; done\" > tmp.sh")
    #return_code = os.system("sh tmp.sh")
    list_gene = OrderedDict()
    sample='sample'
    list_gene[sample] = ''
    text_filenames = glob.glob('Gtfiles/*gtf')
    for filename in text_filenames:
        with open(filename,'r') as fp:
            for line in fp:
                if ('#' not in line) and (line.strip().split('\t')[2] == 'transcript'):
                    id =line.strip().split('\t')[8].split('\"')[3]
                    if id not in list_gene:
                        list_gene[id] = ''
    for filename in text_filenames:
        s=filename.split('.')[0].split('/')[1]
        list_gene[sample] = list_gene[sample]+s+'\t'
        with open(filename,'r') as fp:
            for line in fp:
                if ('#' not in line) and (line.strip().split('\t')[2] == 'transcript') and (line.strip().split('\t')[8].split('\"')[3] in list_gene):
                    fpkm = str(round(float(line.strip().split('\t')[8].split('\"')[7]),2))
                    list_gene[line.strip().split('\t')[8].split('\"')[3]] = list_gene[line.strip().split('\t')[8].split('\"')[3]]+ fpkm + '\t'
    with open('FPKM.txt','w') as fp_out:
        for key,value in list_gene.items():
            fp_out.write('\t'.join([key,value]) + "\n")
    os.chdir(pwd_dir)

def get_distance_for_tss():
    with open('../CAGE/dominant_ctss-summary','r') as f:
        for line in f:
            exec ("%s={}"%line.strip().split('\t')[0])
    line_num = 0
    with open("../CAGE/Minus.tpm",'r') as fp:
        for line in fp:
            line_num += 1
            if (line_num != 1):
                s=line.strip().split('\t')
                locals()[s[0].split("_")[1]][s[0]] = line
                exec ("%s[s[0]]=line"%s[0].split("_")[1])
    w=open('D_minus','w')
    with open('../Transcript/Raw.bed','r') as f:
        for line in f:
            if '-' in line.strip().split('\t')[5]:
                y=line.strip().split('\t')
                site=int(y[2])
                last_exon_length = int(y[10].split(',')[-2])
                regin =min(500,last_exon_length)
                id=y[3]
                for k in locals()[y[0]]:
                    xx = locals()[y[0]][k].strip().split('\t')
                    site_peak=int(xx[0].split('_')[2])
                    peak_id=xx[0]
                    xxx=site_peak - site
                    if (0<=xxx) & (xxx<= 500):     
                        w.write(id+'\t'+peak_id+'\t'+str(xxx)+'\n')
                    if (xxx< 0) & (abs(xxx) <= regin):
                        w.write(id+'\t'+peak_id+'\t'+str(site_peak - site)+'\n')
    w.close()
    with open('../CAGE/dominant_ctss-summary','r') as f:
        for line in f:
            exec ("%s={}"%line.strip().split('\t')[0])
    line_num = 0
    with open("../CAGE/Plus.tpm",'r') as fp:
        for line in fp:
            line_num += 1
            if (line_num != 1):
                s=line.strip().split('\t')
                locals()[s[0].split("_")[1]][s[0]] = line
                exec ("%s[s[0]]=line"%s[0].split("_")[1])
    w=open('D_plus','w')
    with open('../Transcript/Raw.bed','r') as f:
        for line in f:
            if '+' in line.strip().split('\t')[5]:
                    y=line.strip().split('\t')
                    site=int(y[1]) + 1
                    fisrt_exon_length = int(y[10].split(',')[0])
                    regin =min(500,fisrt_exon_length)
                    id=y[3]
                    for k in locals()[y[0]]:
                        xx = locals()[y[0]][k].strip().split('\t')
                        site_peak=int(xx[0].split('_')[2])
                        peak_id=xx[0]
                        xxx=site-site_peak
                        if (0<=xxx) & (xxx<= 500):   
                            w.write(id+'\t'+peak_id+'\t'+str(xxx)+'\n')
                        if (xxx< 0) & (abs(xxx) <= regin):
                            w.write(id+'\t'+peak_id+'\t'+str(site-site_peak)+'\n')
    return_code = os.system("cat D_plus D_minus > DD" ) >> 8
   
def  get_distance_for_tes():
    with open('../PAS/dominant_ctss-summary','r') as f:
        for line in f:
            exec ("%s={}"%line.strip().split('\t')[0])
    line_num = 0
    with open("../PAS/Minus.tpm",'r') as fp:
        for line in fp:
            line_num += 1
            if (line_num != 1):
                s=line.strip().split('\t')
                locals()[s[0].split("_")[1]][s[0]] = line
                exec ("%s[s[0]]=line"%s[0].split("_")[1])
    w=open('D_minus','w')
    with open('../AD_5end/Raw_ad5end.bed','r') as f:
        for line in f:
            if '-' in line.strip().split('\t')[5]:
                y=line.strip().split('\t')
                site=int(y[1]) + 1
                fisrt_exon_length = int(y[10].split(',')[0])
                regin =min(500,fisrt_exon_length)
                id=y[3]
                for k in locals()[y[0]]:
                    xx = locals()[y[0]][k].strip().split('\t')
                    site_peak=int(xx[0].split('_')[2])
                    peak_id=xx[0]
                    xxx=site-site_peak
                    if (0<=xxx) & (xxx<= 500):  
                        w.write(id+'\t'+peak_id+'\t'+str(xxx)+'\n')
                    if (xxx< 0) & (abs(xxx) <= regin):
                        w.write(id+'\t'+peak_id+'\t'+str(xxx)+'\n')

    w.close()
    with open('../PAS/dominant_ctss-summary','r') as f:
        for line in f:
            exec ("%s={}"%line.strip().split('\t')[0])
    line_num = 0
    with open("../PAS/Plus.tpm",'r') as fp:
        for line in fp:
            line_num += 1
            if (line_num != 1):
                s=line.strip().split('\t')
                locals()[s[0].split("_")[1]][s[0]] = line
                exec ("%s[s[0]]=line"%s[0].split("_")[1])
    w=open('D_plus','w')
    with open('../AD_5end/Raw_ad5end.bed','r') as f:
        for line in f:
            if '+' in line.strip().split('\t')[5]:
                y=line.strip().split('\t')
                site=int(y[2]) 
                last_exon_length = int(y[10].split(',')[-2])
                regin = min(500,last_exon_length)
                id=y[3]
                for k in locals()[y[0]]:
                    xx = locals()[y[0]][k].strip().split('\t')
                    site_peak=int(xx[0].split('_')[2])
                    peak_id=xx[0]
                    xxx=site_peak - site
                    if (0<=xxx) & (xxx<= 500):     
                        w.write(id+'\t'+peak_id+'\t'+str(xxx)+'\n')
                    if (xxx< 0) & (abs(xxx) <= regin):
                        w.write(id+'\t'+peak_id+'\t'+str(xxx)+'\n')
    w.close()
    return_code = os.system("cat D_plus D_minus > DD" ) >> 8

def get_cor_between_txs_and_ngs(fp):
    FPKM =  OrderedDict()
    with open('../Transcript/FPKM.txt','r') as f:
        for line in f:
            s = line.strip().split('\t')	
            if 'sample' not in line:
                FPKM[s[0]] = '\t'.join(s[1:])
            else:
                sample = s
    TPM =  OrderedDict()
    line_num = 0  
    index = []
    with open(fp,'r') as f:
        for line in f:
            line_num += 1
            s = line.strip().split('\t')
            if (line_num == 1):
                for m in range(1,len(sample)):
                    for n in range(0,len(s)):
                        if s[n] == sample[m]:
                            index.append(n+1)
                print(index)
            else :
                data = []
                for i in range(0,len(index)):
                    data.append(s[index[i]])
                tpm = '\t'.join(data)
                TPM[s[0]] = tpm
    w=open('pair-with-dis-and-cor','w')
    with open('DD','r') as f:
        for line in f:
            s=line.strip().split('\t')
            iso_id = s[0]
            cage_id = s[1]
            dis = s[2]
            iso_fpkm = FPKM[iso_id].split('\t')
            for n in range(0,len(iso_fpkm)) :
                iso_fpkm[n] = float(iso_fpkm[n])
            cage_tpm = TPM[cage_id].split('\t')
            cage_tpm_num = 0
            cage_tpm_max = 0
            for n in range(0,len(cage_tpm)) :
                cage_tpm[n] = float(cage_tpm[n])
                if cage_tpm[n] >= 0.5 :
                    cage_tpm_num += 1
                if cage_tpm[n] > cage_tpm_max:
                    cage_tpm_max = cage_tpm[n]
            if (sum(iso_fpkm) == 0) or (sum(cage_tpm) == 0) :
                cor = 'NA'
                w.write(iso_id+'\t'+cage_id+'\t'+dis+'\t'+str(cor) +'\t' +str(cage_tpm_num) +'\t' +str(cage_tpm_max)  +'\n')
            else :
                cor=corrcoef(iso_fpkm,cage_tpm)
                w.write(iso_id+'\t'+cage_id+'\t'+dis+'\t'+str(cor) +'\t' +str(cage_tpm_num) +'\t' +str(cage_tpm_max)  +'\n')

def get_real_iso_pairs():
    ID = OrderedDict()
    w=open('iso_peak_pairs','w')
    with open('pair-with-dis-and-cor','r') as f:
        for line in f:
            s=line.strip().split('\t')
            iso_id, cage_id, dis, cor, peak_num, max_peak_tpm = s[0], s[1], s[2], s[3], s[4], s[5]
            if 'NA' not in cor:
                if (iso_id not in ID) and (float(cor) > -0.1):
                    ID[iso_id] = cage_id +'\t' +dis +'\t' +cor+'\t'+peak_num+'\t' +max_peak_tpm
                else :
                    if float(cor) > 0 :
                        if float(cor) > float(ID[iso_id].split('\t')[2]):
                            ID[iso_id] = cage_id +'\t' +dis +'\t' +cor +'\t'+peak_num+'\t' +max_peak_tpm
                        if float(cor) == float(ID[iso_id].split('\t')[2]):
                            if float(peak_num) > float(ID[iso_id].split('\t')[3]) :
                                ID[iso_id] = cage_id +'\t' +dis +'\t' +cor +'\t'+peak_num+'\t' +max_peak_tpm
                            else :
                                if abs(float(dis)) < abs(float(ID[iso_id].split('\t')[1])) :
                                    ID[iso_id] = cage_id +'\t' +dis +'\t' +cor +'\t'+peak_num+'\t' +max_peak_tpm
    for k in ID.keys():
        w.write(k+'\t'+ID[k] +'\n') 

def get_peak_OrderedDict():
    global PEAK_increse,PEAK_decrese,ID
    ID =OrderedDict()
    PEAK_increse =OrderedDict()
    PEAK_decrese =OrderedDict()
    with open('iso_peak_pairs','r') as f:
        for line in f:
            lines = line.strip().split('\t')
            peak_site=int(lines[1].split('_')[2])
            iso_id=lines[0]
            if '-' in lines[2]:
                PEAK_decrese[iso_id] =  peak_site
            else:
                PEAK_increse[iso_id] =  peak_site

def adjust_5end():
    pwd_dir = str(os.path.abspath('.'))
    return_code = os.system("mkdir -p tmp/AD_5end") >> 8
    next_dir = pwd_dir + "/tmp/AD_5end"
    os.chdir(next_dir)
    get_distance_for_tss()
    return_code = os.system("cat ../CAGE/Minus.tpm > CAGE.tpm" ) >> 8   
    return_code = os.system("sed -n '2,$p' ../CAGE/Plus.tpm >> CAGE.tpm" ) >> 8
    get_cor_between_txs_and_ngs("CAGE.tpm")
    get_real_iso_pairs()
    get_peak_OrderedDict()
    w= open('Raw_ad5end.bed','w')
    with open("../Transcript/Raw.bed",'r') as fp:				
        for line in fp:
                data = line.strip().split('\t')
                id = data[3]
                strand = data[5]
                if (id not in PEAK_increse) and (id not in PEAK_decrese):
                    w.write(line)
                elif (id in PEAK_increse):
                    if (strand == '+') :
                        new_start_site = int(PEAK_increse[id]) - 1
                        dert = int(data[1]) - new_start_site 				
                        exon_len = data[10]
                        exon_len_list = exon_len.split(',')
                        exon_len_list[0] = int(exon_len_list[0]) + dert
                        for i in range(0,len(exon_len_list)):
                            exon_len_list[i] = str(exon_len_list[i])
                        new_exon_len=','.join(exon_len_list)
                        exon_site = data[11]
                        exon_site_list=exon_site.split(',')
                        if len(exon_site_list) == 1:
                            new_exon_site = "0,"
                        else:
                            for i in range(1,(len(exon_site_list))-1):
                                exon_site_list[i] = str(int(exon_site_list[i]) + dert )
                        new_exon_site=','.join(exon_site_list)
                        new_data=[data[0],str(new_start_site),'\t'.join(data[2:6]),str(new_start_site),'\t'.join(data[7:10]),new_exon_len,new_exon_site]  
                        w.write('\t'.join(new_data)+'\n')	
                    elif (strand == '-'):
                        new_end_site = int(PEAK_increse[id]) 
                        dert = new_end_site - int(data[2]) 
                        exon_len = data[10]
                        exon_len_list = exon_len.split(',')
                        exon_len_list[-2]=int(exon_len_list[-2])+dert
                        for i in range(0,len(exon_len_list)):
                            exon_len_list[i] = str(exon_len_list[i])
                        new_exon_len=','.join(exon_len_list)
                        exon_site=data[11]
                        new_exon_site=exon_site
                        new_data=[data[0],data[1],str(new_end_site),'\t'.join(data[3:7]),str(new_end_site),'\t'.join(data[8:10]),new_exon_len,new_exon_site]  
                        w.write('\t'.join(new_data)+'\n')	
                elif (id in PEAK_decrese):	
                    if (strand == '+') :
                        new_start_site = int(PEAK_decrese[id]) - 1
                        dert = new_start_site - int(data[1]) 	
                        exon_len = data[10]
                        exon_len_list = exon_len.split(',')
                        exon_len_list[0] = int(exon_len_list[0]) - dert
                        for i in range(0,len(exon_len_list)):
                            exon_len_list[i] = str(exon_len_list[i])
                        new_exon_len=','.join(exon_len_list)		
                        exon_site = data[11]
                        exon_site_list=exon_site.split(',')
                        if len(exon_site_list) == 1:
                            new_exon_site = "0,"
                        else:
                            for i in range(1,(len(exon_site_list))-1):
                                exon_site_list[i] = str(int(exon_site_list[i]) - dert )
                        new_exon_site=','.join(exon_site_list)												
                        new_data=[data[0],str(new_start_site),'\t'.join(data[2:6]),str(new_start_site),'\t'.join(data[7:10]),new_exon_len,new_exon_site]   
                        w.write('\t'.join(new_data)+'\n')					
                    elif (strand == '-'):
                        new_end_site = int(PEAK_decrese[id]) 
                        dert = int(data[2]) - new_end_site 		
                        exon_len = data[10]
                        exon_len_list = exon_len.split(',')
                        exon_len_list[-2]=int(exon_len_list[-2]) - dert
                        for i in range(0,len(exon_len_list)):
                            exon_len_list[i] = str(exon_len_list[i])
                        new_exon_len=','.join(exon_len_list)
                        exon_site=data[11]
                        new_exon_site=exon_site
                        new_data=[data[0],data[1],str(new_end_site),'\t'.join(data[3:7]),str(new_end_site),'\t'.join(data[8:10]),new_exon_len,new_exon_site]  
                        w.write('\t'.join(new_data)+'\n')		  
    os.chdir(pwd_dir)

def asjust_3end():
    pwd_dir = str(os.path.abspath('.'))
    return_code = os.system("mkdir -p tmp/AD_3end") >> 8
    next_dir = pwd_dir + "/tmp/AD_3end"
    os.chdir(next_dir)
    get_distance_for_tes()
    return_code = os.system("cat ../PAS/Minus.tpm > PAS.tpm" ) >> 8   
    return_code = os.system("sed -n '2,$p' ../PAS/Plus.tpm >> PAS.tpm" ) >> 8
    get_cor_between_txs_and_ngs("PAS.tpm")
    get_real_iso_pairs()
    get_peak_OrderedDict()
    w=open('Raw_ad5end_ad3end.bed','w')
    with open("../AD_5end/Raw_ad5end.bed",'r') as fp:                         
        for line in fp:
                data = line.strip().split('\t')
                id = data[3]
                strand = data[5]
                if (id not in PEAK_increse) and (id not in PEAK_decrese):
                    w.write(line)
                elif (id in PEAK_increse):
                    if (strand == '-') :
                        new_start_site = int(PEAK_increse[id]) -1
                        dert = int(data[1]) - new_start_site 
                        exon_len = data[10]
                        exon_len_list = exon_len.split(',')
                        exon_len_list[0] = int(exon_len_list[0]) + dert
                        for i in range(0,len(exon_len_list)):
                            exon_len_list[i] = str(exon_len_list[i])
                        new_exon_len=','.join(exon_len_list)
                        exon_site = data[11]
                        exon_site_list=exon_site.split(',')
                        if len(exon_site_list) == 1:
                            new_exon_site = "0,"
                        else:
                            for i in range(1,(len(exon_site_list))-1):
                                exon_site_list[i] = str(int(exon_site_list[i]) + dert)
                        new_exon_site=','.join(exon_site_list)
                        new_data=[data[0],str(new_start_site),'\t'.join(data[2:6]),str(new_start_site),'\t'.join(data[7:10]),new_exon_len,new_exon_site]
                        w.write('\t'.join(new_data)+'\n')       
                    elif (strand == '+'):
                        new_end_site = int(PEAK_increse[id]) 
                        dert = new_end_site - int(data[2]) 
                        exon_len = data[10]
                        exon_len_list = exon_len.split(',')
                        exon_len_list[-2]=int(exon_len_list[-2])+dert
                        for i in range(0,len(exon_len_list)):
                            exon_len_list[i] = str(exon_len_list[i])
                        new_exon_len=','.join(exon_len_list)
                        exon_site=data[11]
                        new_exon_site=exon_site
                        new_data=[data[0],data[1],str(new_end_site),'\t'.join(data[3:7]),str(new_end_site),'\t'.join(data[8:10]),new_exon_len,new_exon_site]  
                        w.write('\t'.join(new_data)+'\n')      
                elif (id in PEAK_decrese):    
                    if (strand == '-') :
                        new_start_site = int(PEAK_decrese[id]) - 1
                        dert = new_start_site - int(data[1]) 
                        exon_len = data[10]
                        exon_len_list = exon_len.split(',')
                        exon_len_list[0] = int(exon_len_list[0]) - dert
                        for i in range(0,len(exon_len_list)):
                            exon_len_list[i] = str(exon_len_list[i])
                        new_exon_len=','.join(exon_len_list)
                        exon_site = data[11]
                        exon_site_list=exon_site.split(',')
                        if len(exon_site_list) == 1:
                            new_exon_site = "0,"
                        else:
                            for i in range(1,(len(exon_site_list))-1):
                                exon_site_list[i] = str(int(exon_site_list[i]) - dert )
                        new_exon_site=','.join(exon_site_list)
                        new_data=[data[0],str(new_start_site),'\t'.join(data[2:6]),str(new_start_site),'\t'.join(data[7:10]),new_exon_len,new_exon_site]  
                        w.write('\t'.join(new_data)+'\n')    
                    elif (strand == '+'):
                        new_end_site = int(PEAK_decrese[id]) 
                        dert = int(data[2]) - new_end_site 
                        exon_len = data[10]
                        exon_len_list = exon_len.split(',')
                        exon_len_list[-2]=int(exon_len_list[-2]) - dert
                        for i in range(0,len(exon_len_list)):
                            exon_len_list[i] = str(exon_len_list[i])
                        new_exon_len=','.join(exon_len_list)
                        exon_site=data[11]
                        new_exon_site=exon_site
                        new_data=[data[0],data[1],str(new_end_site),'\t'.join(data[3:7]),str(new_end_site),'\t'.join(data[8:10]),new_exon_len,new_exon_site]  
                        w.write('\t'.join(new_data)+'\n')     
    w=open('Raw_ad5end_ad3end-final.bed','w')
    with open('Raw_ad5end_ad3end.bed','r') as f:
        for line in f:
            s=line.strip().split('\t')
            if (int(s[10].split(',')[0]) == 0) or (int(s[10].split(',')[-2]) == 0 ) or (int(s[2]) - int(s[1]) +1 >= 8000) or (int(s[2]) <= int(s[1])):
                continue
            else:
                w.write(line) 
    os.chdir(pwd_dir)

def main():
    if len(sys.argv) == 1:
        sys.exit(__doc__)
    else:
        options = docopt(__doc__, version=__version__)
        print(options)
        return_code = os.system("mkdir tmp") >> 8
        geno_seq = str(os.path.abspath((options['-g'])))
        geno_ref = str(os.path.abspath((str(options['-r']))))
        ngs_bam_dir = str(os.path.abspath((str(options['--ngs_bam']))))
        #get_BSgenome(geno_seq)
        if options['--ngs'] and options['--pb']:
            combin_ngs_and_pacbio(str(options['--ngs']),str(options['--pb']),geno_ref)
        elif options['--ngs'] and not options['--pb']:
            only_ngs(str(options['--ngs']))
        elif not options['--ngs'] and options['--pb']:
            only_pacbio(str(options['--pb']))
        #calculate_fpkm_for_iso(ngs_bam_dir)

        if options['-c'] : 
            print('1')
            #extract_dominant_ctss(str(options['-c']),"CAGE")
            #get_txs_cluster(str(options['-c']),"CAGE")
            adjust_5end()
        else :
            return_code = os.system("mkdir -p tmp/AD_5end") >> 8
            return_code = os.system("cp tmp/Transcript/Raw.bed tmp/AD_5end/Raw_ad5end.bed") >> 8

        if options['-p'] :
            print('1')
            #extract_dominant_ctss(str(options['-p']),"PAS")
            #get_txs_cluster(str(options['-p']),"PAS")
            asjust_3end()
        else :
            return_code = os.system("mkdir -p tmp/AD_3end") >> 8
            return_code = os.system("cp tmp/AD_5end/Raw_ad5end.bed tmp/AD_3end/Raw_ad5end_ad3end-final.bed") >> 8
        os.chdir("./tmp/")
        return_code = os.system("ln -s " + geno_seq + " GENOME.fa") >> 8
        script_dir = str(os.popen("find / -name Pull.py ").read().strip()).replace("Pull.py","script") + "/"
        return_code = os.system(str(script_dir)+"gtf2bed"+ " " + geno_ref +" > pcg.bed") >> 8
        return_code = os.system("ln -s tmp/tmp/AD_3end/Raw_ad5end_ad3end-final.bed Raw_lnc.bed") >> 8
        return_code = os.system("snakemake -s " + str(script_dir)+"Lnc_identity.sm.rules.py") >>8


if __name__ == '__main__':
    main()	
