Plant Full-Length LncRNA (PULL) 

This repository provides a pipline for plant full-length lncRNA identity called PULL. Which take advantages of sequencing for terminal signal to get a relatively complete transcript units collection. Additionally, long read sequencing data can also be added to form a high credible collection. After get the full-length transcriptome, we combined several popular tools to evaluate the coding potential in order to get the widespread non-coding RNAs, especially the long-noncoding RNAs (lncRNAs). The pipline built on python and some other softwares, so we recommend you use anaconda to create a specific work environment. The repository can only work on linux now.

Step1： Downlaod anaconda
wget https://repo.anaconda.com/archive/Anaconda2-2018.12-Linux-x86_64.sh
sh Anaconda2-2018.12-Linux-x86_64.sh

Step2: Creat a new python3 work environment
conda create -n Py3_pull -c bioconda --file requirements_py3.txt python=3

Step3: Installation of PULL
git clone https://github.com/cotton-lab/PULL.git 
Unpack the compressed file, then PULL is ready for use!

Besides, here are also some softwares you are supposed to installed already.
ncbi-blast(v2.7.1+) 
stringtie(v1.3.4d) 
R(v3.4.2) 
CAGEr(v1.22.3)
