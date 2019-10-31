import os
from Bio import SeqIO
from Bio.Seq import Seq
from collections import OrderedDict

WORK_BASE = "tmp"

rule all :
	input:
		a = os.path.join(WORK_BASE,"blastn_pcg.txt"),
		b = os.path.join(WORK_BASE,"Raw-1.fa"),
		c = os.path.join(WORK_BASE,"blastn_hk.txt"),
		d = os.path.join(WORK_BASE,"Raw-2.fa"),
		e = os.path.join(WORK_BASE,"orf_out"),
		f = os.path.join(WORK_BASE,"Raw-3.fa"),
		g = os.path.join(WORK_BASE,"cpc-out"),
		h = os.path.join(WORK_BASE,"CNCI_out"),
		i = os.path.join(WORK_BASE,"Raw-4.fa"),
		j = os.path.join(WORK_BASE,"Blastx.txt"),
		k = os.path.join(WORK_BASE,"pfam-chuli.txt"),
		q = os.path.join(WORK_BASE,"Raw-5.fa"),
################################################################################
# Exclude known protein coding genes
################################################################################
rule exclude_known_gene :
	input:
        	rawlnc = "Raw_lnc.bed" ,
	        genome_seq = "GENOME.fa" ,
		pcg_annovation = "pcg.bed" ,
	output:
		twobit = os.path.join(WORK_BASE, "GENOME.2bit"),
		lncfa = os.path.join(WORK_BASE,"Raw_lnc.fa"),
		pcgfa = os.path.join(WORK_BASE,"pcg.fa"),
		blastn1 = os.path.join(WORK_BASE,"blastn_pcg.txt"),
		blastn2 = os.path.join(WORK_BASE,"blastn_pcg_reverse.txt"),
	threads:10
	log : 
		"tmp/exclude_known_gene.log",
	params:
		job_name = "exclude_known_gene",
		log = "log",
	shell:"""
faToTwoBit {input.genome_seq} {output.twobit}
twoBitToFa -bed={input.rawlnc} {output.twobit} {output.lncfa}
twoBitToFa -bed={input.pcg_annovation} {output.twobit} {output.pcgfa}
mkdir -p tmp/PCG_DB
makeblastdb -in {output.pcgfa}  -dbtype  nucl -parse_seqids -out tmp/PCG_DB/db
blastn -query {output.lncfa} -strand plus -db tmp/PCG_DB/db -out {output.blastn1} -outfmt "6 qseqid qlen sseqid slen qstart qend evalue qcovs pident"  -max_target_seqs 1 -num_threads 10 -evalue 1e-6 
mkdir -p tmp/LNC_DB
makeblastdb -in {output.lncfa} -dbtype  nucl -parse_seqids -out tmp/LNC_DB/db
blastn -query {output.pcgfa} -strand  plus -db tmp/LNC_DB/db -out {output.blastn2} -outfmt "6 qseqid qlen sseqid slen  qstart qend evalue qcovs pident" -max_target_seqs 1 -num_threads 10 -evalue 1e-6
"""

rule deal1:
	input :
		blastn1 = rules.exclude_known_gene.output.blastn1,
		blastn2 = rules.exclude_known_gene.output.blastn2,
	output:
		raw1 = os.path.join(WORK_BASE,"Raw-1.fa"),
	run :
		ID = OrderedDict()
		w = open(output.raw1,'w') 
		with open(input.blastn1,'r') as f:
			for line in f:
				s = line.strip().split('\t')
				if (int(s[7]) >= 30) and  (float(s[8]) >= 80):
					ID[s[0]] = ''
		with open(input.blastn2,'r') as f:
			for line in f :
				m=line.strip().split('\t')
				if (int(m[7]) >= 30) and  (float(m[8]) >= 80) :
					ID[m[2]] = ''  
		for seq_record in SeqIO.parse(rules.exclude_known_gene.output.lncfa,'fasta'):
	        	if str(seq_record.id) not in ID:
                		w.write('>'+ str(seq_record.id) +'\n' + str(seq_record.seq) +'\n')

################################################################################
# Exclude House-Keeping genes
################################################################################
rule exclude_hk_gene :
	input:
		hkseq = "data/HK.fa",
		raw1 = rules.deal1.output.raw1,
	output:
		blastn3 = os.path.join(WORK_BASE,"blastn_hk.txt"),
	threads: 10
	log :
                "tmp/exclude_hk_gene.log",
	params:
        	job_name = "exclude_hk_gene",
	    	log = "log",
	shell:"""
mkdir -p tmp/HK_DB
makeblastdb -in {input.hkseq} -dbtype  nucl -parse_seqids -out tmp/HK_DB/db
blastn -query {input.raw1} -strand plus -db tmp/HK_DB/db -out {output.blastn3} -outfmt "6 qseqid qlen sseqid slen qstart qend evalue qcovs pident" -evalue 1e-6 -num_threads 10 -max_target_seqs 1
	"""

rule deal2:
	input:	
		blastn3 = rules.exclude_hk_gene.output.blastn3,
	output:
		raw2 = os.path.join(WORK_BASE,"Raw-2.fa"),
	run :
		ID = OrderedDict()
		w=open(output.raw2,'w')
		with open(input.blastn3,'r') as f:
			for line in f:
				s = line.strip().split('\t')	
				if (int(s[7]) >= 30) and (float(s[8]) >= 80) :
					ID[s[0]] = ''
		for seq_record in SeqIO.parse(rules.exclude_hk_gene.input.raw1,'fasta'):
			if str(seq_record.id) not in ID:
				w.write('>'+ str(seq_record.id) +'\n' + str(seq_record.seq) +'\n')
			else :
				continue

################################################################################
# Filter orf and length
################################################################################
rule filter_orf_and_len :
	input :
		raw2 = rules.deal2.output.raw2,
	output:
		orf = os.path.join(WORK_BASE,"orf_out"),
	threads:10
	log :
                "tmp/filter_orf_and_len.log",
	params:
        	job_name = "filter_orf_and_len",
	shell:"""
getorf -sequence {input.raw2}  -outseq {output.orf} -minsize 300 -find 1 -noreverse
"""
rule deal3:
	input :
		orf = rules.filter_orf_and_len.output.orf,
	output:
		raw3 = os.path.join(WORK_BASE,"Raw-3.fa"),
	run:
		ID = OrderedDict()
		w = open(output.raw3,'w')
		for r in SeqIO.parse(input.orf, "fasta"):
			ID[str(r.id).split('_')[0]] = ''
		for r in SeqIO.parse(rules.deal2.output.raw2,"fasta"):
			s=(r.seq)
			l = len(s)
			id =str(r.id)
			if (l < 200) or (id in ID) :
				continue
			else :
				w.write('>'+ str(r.id) +'\n' + str(r.seq) +'\n')

################################################################################
# Filter by cpc2 and cnci
################################################################################
rule cpc_and_cnci :
	input :
		raw3 = rules.deal3.output.raw3,
	output:
		cpc_result = os.path.join(WORK_BASE,"cpc-out"),
		cnci_result = os.path.join(WORK_BASE,"CNCI_out"),
	threads:10
	log :
                "tmp/cpc_cnci.log",
	params:
        	job_name = "cpc_and_cnci",
	shell:"""
python essentials/CPC2-beta/bin/CPC2.py -i {input.raw3} -o {output.cpc_result}
python essentials/CNCI-master/CNCI.py -f {input.raw3} -o {output.cnci_result} -m pl -p 10
"""

rule deal4:
	input :
		cpc_result = rules.cpc_and_cnci.output.cpc_result,
	output:
		raw4 = os.path.join(WORK_BASE,"Raw-4.fa"),
	run:
		ID = OrderedDict()
		ID_all = OrderedDict()
		w=open(output.raw4,'w')
		with open(input,'r') as f:
			for line in f:
				if '#' not in line:
					s=line.strip().split('\t')
					ID_all[s[0]] = ''
					if "coding" == s[7]:
						ID[s[0]] = ''
		with open('tmp/CNCI_out/CNCI.index','rb') as f:
			for line in f:
				lines= line.strip().split('\t')
				if float(lines[2] >= -0.0001) == lines[1] :
					ID[lines[0]] = ''
		with open('tmp/CNCI_out.log','rb') as f:
			for line in f:
				i=line.strip().split(' ')[0].split(">")[1]
				ID[i] = ''
		for seq_record in SeqIO.parse(rules.deal3.output.raw3,'fasta'):
			if seq_record.id not in ID:
				w.write('>'+ str(seq_record.id) +'\n' + str(seq_record.seq) +'\n')

################################################################################
# Exclude hit to known protein sequence against db
################################################################################
rule swiss_prot:
	input:
                proseq = "data/uniprot_sprot.fasta",
		raw4 = rules.deal4.output.raw4,
	output:
		blastx = os.path.join(WORK_BASE,"Blastx.txt"),
	threads:10
	log:
		"tmp/swiss_prot.log",
	params:
		job_name = "swiss_prot",
	shell:"""
mkdir -p tmp/Swisspro
makeblastdb -in {input.proseq} -dbtype prot -parse_seqids -out tmp/Swisspro/db
blastx -db tmp/Swisspro/db -query {input.raw4} -out {output.blastx} -strand plus -evalue 1e-6 -max_target_seqs 1 -num_threads 50 -outfmt "6 qseqid qlen sseqid slen qstart qend evalue qcovs pident"
"""

################################################################################
# Exclude hit to known Protein domain
################################################################################
rule filter_pfam :
	input:
		raw4 = rules.deal4.output.raw4 ,
	output:
		raw4_pep = os.path.join(WORK_BASE,"Raw-4.pep.fa"),
		pfam_out = os.path.join(WORK_BASE,"Pfam_out"),
		chuli = os.path.join(WORK_BASE,"pfam-chuli.txt"),
	log:
                "filter_pfam.log",
	params:
		job_name = "filter_pfam",
	shell:"""
getorf -sequence {input.raw4} -outseq {output.raw4_pep} -minsize 60 -find 1 -noreverse
essentials/hmmer-3.1b2/binaries/hmmscan --notextw --cut_ga --cpu 10 --tblout {output.pfam_out} essentials/DA/Pfam-A.hmm {output.raw4_pep}
awk '{{FS=" "; OFS="\t"}} {{if ($0!~/#/) print $3,$8}}' {output.pfam_out}  > {output.chuli}
"""

rule deal5:
	input :
		blastx = rules.swiss_prot.output.blastx,
		chuli = rules.filter_pfam.output.chuli,
	output:
		raw5 = os.path.join(WORK_BASE,"Raw-5.fa"),
	run:
		ID = OrderedDict()
		I=OrderedDict()
		w=open(output.raw5,'w')
		with open(input.blastx,'r') as f:
			for line in f:
				lines= line.strip().split('\t')
				if (float(lines[8]) >= 50): 
					ID[lines[0]] = '' 

		with open(input.chuli,'rb') as f:
			for line in f:
				if 'name' not in line:
					lines= line.strip().split('\t')
					if float(lines[1]) <= 0.001:
						ID[lines[0].split('.')[0]]= ''
					
		for seq_record in SeqIO.parse(rules.deal4.output.raw4,'fasta'):
        		if str(seq_record.id) not in ID:
                		w.write('>'+str(seq_record.id)+'\n'+str(seq_record.seq) +'\n')

