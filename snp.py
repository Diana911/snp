# -*- coding: utf-8 -*-
"""
Spyder Editor
import pandas as pd
#import math
#import seaborn as sns
This is a temporary script file.
"""

import io
import os
import pandas as pd
from Bio import SeqIO
import re
from gtfparse import read_gtf
import pybedtools
import numpy as np


def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

def getref(st):
    return re.findall(r'[^:]+$', st)[0].split(',')[0]

def getalt(st):
    return re.findall(r'[^:]+$', st)[0].split(',')[1]

def trim(srr, path, adap):
        os.system("TrimmomaticPE -threads 10 -phred33 "+path+srr+"_1.fastq.gz "+path+srr+"_2.fastq.gz "+\
        path+srr+"_1paired.fastq"+" "+path+srr+"_1unpaired.fastq "+path+srr+"_2paired.fastq "+\
        path+srr+"_2unpaired.fastq ILLUMINACLIP:"+adap+" LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36")
        os.system("rm "+path+"*unpaired.fastq")

def get_bam(srr, path, ind, fl):
    os.system("bowtie2 --no-unal -p30 -x "+ind+" -1 "+path+srr+"_1paired.fastq -2 "+path+srr+"_2paired.fastq \
              | samtools view -q 20 -F 4 -Sb | samtools sort -o "+path+srr+fl+".bam")
    
    

def get_vcf(genome, path, sample, bam, fl):
    
    os.system("bcftools mpileup -f "+genome\
              +" --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR "+\
              bam+" | bcftools call -cv -Oz > "+path+"var_calls"+sample+fl+".vcf.gz")
    os.system("gunzip " +path+ "var_calls"+sample+fl+".vcf.gz")
    
    
def for_new_genome(path, sample, fl):
    datatest = read_vcf(path+"var_calls"+sample+fl+".vcf")
    datatest['ref']=datatest.filter(regex='.bam', axis=1).iloc[:,0].apply(getref)
    datatest['alt']=datatest.filter(regex='.bam', axis=1).iloc[:,0].apply(getalt)
    datatest['ref'] = datatest['ref'].astype(float)
    datatest['alt'] = datatest['alt'].astype(float)
    datatest['log'] = np.log10(datatest['alt']/(datatest['ref']+datatest['alt']))
    c={}
    h=datatest.loc[(datatest.ref > 10) & (datatest.alt > 10) & (datatest.log<1) & (datatest.log>-1) &\
                   ~ (datatest['INFO'].str.contains('INDEL'))]
    ind=h.loc[:,['CHROM','POS']]
    for i, row in h.iterrows():
        if row['CHROM'] not in c:
            c[row['CHROM']]=[[row['POS'],row['ALT']]]
        else:
            c[row['CHROM']].append([row['POS'],row['ALT']])
    return [c, ind]
            
            
 
    
adap = '/usr/share/trimmomatic/TruSeq2-PE.fa:2:30:10'    
genome = '/media/leon/WorkDATA/grch38_snp/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
conf = pd.read_csv("conf.csv")
os.system("bowtie2-build " + genome + " ind")
ind = '/media/leon/DISK2/DataforDiana/ind'


for srr in conf.file_name:
    path = conf.loc[conf['file_name'] == srr, "path"].values[0]
    sample = conf.loc[conf['file_name'] == srr, "biosample"].values[0]
    trim(srr, path, adap)
    get_bam(srr,path,ind,'ref')

for sample in conf.biosample.unique():
    path = conf.loc[conf['biosample'] == sample, "path"].values[0]
    os.system('samtools merge '+path+'allbam'+sample+'.bam '+path+'*.bam')
    get_vcf(genome,path,sample,path+'allbam*','all')


for sample in conf.biosample.unique():
    path = conf.loc[conf['biosample'] == sample, "path"].values[0]
    c = for_new_genome(path, sample,"all")[0]
    with open (path+"newsequence"+sample+".fasta", "w") as newfasta:
        for seq_record in SeqIO.parse(genome, "fasta"):
            if seq_record.id in c:
                lt=list(seq_record.seq)
                for i in c[seq_record.id]:
                    lt[i[0]-1]=i[1]
                st = ''.join(lt)
                newfasta.write(">" + seq_record.description + "\n")
                newfasta.write(st + "\n")
            else:
                newfasta.write(">" + seq_record.description + "\n")
                newfasta.write(str(seq_record.seq) + "\n")
    os.system("bowtie2-build "+path+"newsequence" +sample+ ".fasta " +path+sample+ "ind")
    
for srr in conf.file_name:
    path = conf.loc[conf['file_name'] == srr, "path"].values[0]
    sample = conf.loc[conf['file_name'] == srr, "biosample"].values[0]
    new_ind = path+sample+"ind"
    get_bam(srr,path,new_ind,'alt')
    
for sample in conf.biosample.unique(): 
    path = conf.loc[conf['biosample'] == sample, "path"].values[0]
    new_gen = path+"newsequence"+sample+".fasta"
   
    list1=[]
    list2=[]
    for srr in conf.loc[conf['biosample'] == sample, "file_name"]:  
        alt_srr = path+srr+"alt.bam"
        ref_srr = path+srr+"ref.bam"
        list1.append(alt_srr)
        list2.append(ref_srr)
        string1 = ' '.join([str(item) for item in list1])
        string2 = ' '.join([str(item) for item in list2])
    get_vcf(new_gen,path,sample,string1,'alt')
    get_vcf(genome,path,sample,string2,'ref')
    
homsap = read_gtf('/media/leon/DISK2/DataforDiana/Homo_sapiens.GRCh38.104.gtf')

def chrr(st):
    return 'chr'+re.findall(r'[^.]+$', st)[0]
homsap['seqname']=homsap['seqname'].apply(chrr)

homsap = homsap.loc[:,['seqname','start','end','gene_id']]
homsap['start_prom']=homsap['start']-1500
homsap['end_prom']=homsap['start']+1500
homsap=homsap[['seqname','start_prom','end_prom','start','end','gene_id']]
homsap=homsap.loc[homsap['seqname']!='chrMT']
homsap_pb = pybedtools.BedTool.from_dataframe(homsap)
homsap=homsap.columns.tolist()

def chrrr(num):
    return 'chr'+num

for sample in conf.biosample.unique():
    path = conf.loc[conf['biosample'] == sample, "path"].values[0]
    althet=read_vcf(path + "var_calls"+sample+"alt.vcf")
    althet=althet.rename({'REF':'REF_in_alt', 'ALT':'ALT_in_alt'},axis='columns').drop(['INFO','FORMAT','FILTER','QUAL'], 
                                                                                       axis=1) 
    datatest=read_vcf(path + "var_calls"+sample+"ref.vcf")    
    
    
    for srr in conf.loc[conf['biosample'] == sample, "file_name"]:

        ty=conf.loc[conf['file_name'] == srr, "method"].values[0]
        datatest[ty+srr+'ref']=datatest.filter(regex=srr+'ref.bam', axis=1).iloc[:,0].apply(getref).astype(float) 
        althet[ty+srr+'altt']=althet.filter(regex=srr+'alt.bam', axis=1).iloc[:,0].apply(getref).astype(float)
    
    newall=datatest.merge(althet, how='inner')
    newall=newall.merge(for_new_genome(path, sample,"all")[1]).to_csv(path+'het'+sample+'.vcf',sep='\t',index=False)
    
    
    os.system("java -Xmx12G -jar /media/leon/DISK2/DataforDiana/SnpSift.jar annotate -id /media/leon/DISK2/DataforDiana/00-All.vcf.gz "+path+"het"+sample+".vcf > "+path+sample+
              "_annotated_het.vcf")
    
    
    newall=read_vcf(path+sample+"_annotated_het.vcf")
    newall=newall.filter(regex="^(?:(?!bam).)+$", axis=1)
    newall=newall.loc[newall.ID!='.'].reset_index(drop=True)
    
    newall=newall.drop(['0.0','FILTER','INFO','FORMAT'], axis=1).dropna().reset_index(drop=True)
    newall=newall.rename(columns={"0": "POS"})
    newall['scope1']=newall.POS-200
    newall['scope2']=newall.POS+200
    newall['ro']=0
    
    for l in range(len(newall.POS)):
        for m in range(l-7,l+7):
            if m<0:
                continue
            if m>=len(newall.POS):
                continue
            if newall.POS[m] in range(newall.scope1[l], newall.scope2[l]):
                newall.ro[l]=newall.ro[l]+1
    
    newall=newall.drop(['scope1','scope2'], axis=1).loc[newall['CHROM']!='chrMT']
    newall['CHROM']=newall.CHROM.astype(str).apply(chrrr)
    newall['POS1']=newall['POS']
    cols = newall.columns.tolist()
    cols = cols[0:2]+[cols[-1]]+cols[2:-1]
    newall = newall[cols]
    
    
    chip=newall[(newall.filter(regex='Chip_seq', axis=1)>10).any(axis=1)].filter(regex="^(?:(?!RNA_seq).)+$", axis=1)
    chip_pb = pybedtools.BedTool.from_dataframe(chip)
    my_and_ref = homsap_pb.intersect(chip_pb)
    ref_and_my = chip_pb.intersect(homsap_pb)
    my_and_ref_df=pybedtools.BedTool.to_dataframe(my_and_ref, names=homsap)
    ref_and_my_df=pybedtools.BedTool.to_dataframe(ref_and_my, names=chip.columns.tolist())
    my_and_ref_df['start_prom']=my_and_ref_df['start_prom']+1 #немного костыль
    merge_my_ref = pd.merge(my_and_ref_df, ref_and_my_df, left_on=['seqname','start_prom'], 
                            right_on=['CHROM','POS'], how='left')
    merge_my_ref = merge_my_ref.drop(['POS1'], axis=1)
    merge_my_ref = merge_my_ref.drop_duplicates(subset=['seqname','start_prom'], keep='first')
    
    rna=newall[(newall.filter(regex='RNA_seq', axis=1)>30).any(axis=1)].rename({'ID':'ID_RNA','ro':'ro_RNA',
                                                                                'POS':'POS_RNA','REF':'REF_RNA',
                                                                                'ALT':'ALT_RNA',
                                                                                'REF_in_alt':'REF_in_alt_RNA',
                                                                                'ALT_in_alt':'ALT_in_alt_RNA',
                                                                                'CHROM':'CHROM_RNA'},axis='columns').filter(
                                                                                   regex="^(?:(?!Chip_seq).)+$", axis=1)
    
    cols = merge_my_ref.columns.tolist()
    cols = [cols[0]]+cols[3:5]+cols[1:3]+cols[5:]
    merge_my_ref = merge_my_ref[cols]
    
    gene = pybedtools.BedTool.from_dataframe(merge_my_ref)
    rna_pb = pybedtools.BedTool.from_dataframe(rna)
    rna_and_gene = rna_pb.intersect(gene)
    gene_and_rna = gene.intersect(rna_pb)
    
    rna_and_gene_df=pybedtools.BedTool.to_dataframe(rna_and_gene,names=rna.columns.tolist())
    gene_and_rna_df=pybedtools.BedTool.to_dataframe(gene_and_rna,names=merge_my_ref.columns.tolist())
    gene_and_rna_df['start']=gene_and_rna_df['start']+1 #немного костыль
    
    merge = pd.merge(gene_and_rna_df, rna_and_gene_df, right_on=['CHROM_RNA','POS_RNA'], left_on=['seqname','start'], 
                     how='left')
    merge = merge.drop(['POS1'], axis=1)
    merge.to_csv(path+'table'+sample+'.csv',index=False)