#!bin/bash

conda create -y -n ngs -c bioconda fastqc trimmomatic bwa samtools varscan bedtools

#--------- Installation ----------------
# Installe tout les outils nécessaire pour l'analyse

# Téléchargement des fast
wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1DM9g8OulE1ScBk-HoaREfUZs6gurtkBe' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1DM9g8OulE1ScBk-HoaREfUZs6gurtkBe" -O patient7.tar.gz && rm -rf /tmp/cookies.txt

tar -zxvf patient7.tar.gz

# Crée le dossier index si pas déja présent
if [ ! -f index ]; then
	mkdir index
fi

# Téléchargement de l'index du chromosome 16
wget -P index/ http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr16.fa.gz

# Crée le dossier annotation si pas déja présent
if [ ! -d annotation ]; then
	mkdir annotation
fi

# Téléchargement de l'annotation du génome
wget -P annotation/ ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/GRCh37_mapping/gencode.v24lift37.basic.annotation.gtf.gz
