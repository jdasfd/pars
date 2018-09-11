# Processing Yeast PARS Data

[TOC level=1-3]: # " "
- [Processing Yeast PARS Data](#processing-yeast-pars-data)
- [Download reference data](#download-reference-data)
    - [Download PARS10 full site.](#download-pars10-full-site)
    - [Download S288c annotation data from ensembl by rsync](#download-s288c-annotation-data-from-ensembl-by-rsync)
    - [SGD](#sgd)
- [Download strains and outgroups](#download-strains-and-outgroups)
    - [Sanger (NCBI WGS)](#sanger-ncbi-wgs)
    - [Illumina (NCBI ASSEMBLY)](#illumina-ncbi-assembly)
- [Plans of alignments](#plans-of-alignments)
- [AlignDB](#aligndb)
    - [Build alignDB for multiple genomes](#build-aligndb-for-multiple-genomes-n7)
    - [Extract gene-list and snp-codon-list](#extract-gene-list-and-snp-codon-list-n7)
    - [SNPs and indels](#snps-and-indels-n7)
- [Blast](#blast)
- [Features](#features)
- [Real Processing](#real-processing)
- [Download other reference data](#download-other-reference-data)
    - [mRNA levels](#mrna-levels)
    - [ess, rich/minimal and chem](#ess-richminimal-and-chem)
    - [Recombination rates](#recombination-rates)
    - [Protein-protein interactions](#protein-protein-interactions)
- [Phylogeny](#phylogeny)
    - [create protein coding gene list](#create-protein-coding-gene-list)
    - [cut cds alignment](#cut-cds-alignment)
        - [create cds_yml](#create-cds_yml)
        - [cut cds_alignment by cds_yml](#cut-cds_alignment-by-cds_yml)
        - [count cds_alignment proporation in sgd](#count-cds_alignment-proporation-in-sgd)
    - [create gene_phylogeny (n157_nonMosaic)](#create-gene_phylogeny-n157_nonmosaic)
    - [count distance (n157_nonMosaic)](#count-distance-n157_nonmosaic)
- [SNP](#snp)

# Download PARS10 full site.

```bash
mkdir -p ~/data/mrna-structure/PARS10
cd ~/data/mrna-structure/PARS10

perl ~/Scripts/download/list.pl -u http://genie.weizmann.ac.il/pubs/PARS10/
perl ~/Scripts/download/download.pl -i pubs_PARS10.yml

find . -name "*.gz" | xargs gzip -d

```

# Download strains genomes

## Download S288c (soft-masked) from Ensembl

```bash
mkdir -p ~/data/alignment/egaz/download/S288c
cd ~/data/alignment/egaz/download/S288c

aria2c -x 6 -s 3 -c ftp://ftp.ensembl.org/pub/release-82/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz
aria2c -x 6 -s 3 -c ftp://ftp.ensembl.org/pub/release-82/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.82.gff3.gz
find . -name "*.gz" | xargs gzip -t
faops filter -N -s Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz S288c.fa

```

## Download strains from NCBI assembly

**DBVPG6044**
```bash
mkdir ~/data/alignment/egaz/download/DBVPG6044
cd ~/data/alignment/egaz/download/DBVPG6044
aria2c -x 6 -s 3 -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/079/025/GCA_002079025.1_ASM207902v1/GCA_002079025.1_ASM207902v1_genomic.fna.gz
find . -name "*.gz" | xargs gzip -t
faops filter -N -s GCA_002079025.1_ASM207902v1_genomic.fna.gz DBVPG6044.fa
```

**Y12**
```bash
mkdir ~/data/alignment/egaz/download/Y12
cd ~/data/alignment/egaz/download/Y12
aria2c -x 6 -s 3 -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/058/645/GCA_002058645.1_ASM205864v1/GCA_002058645.1_ASM205864v1_genomic.fna.gz
find . -name "*.gz" | xargs gzip -t
faops filter -N -s GCA_002058645.1_ASM205864v1_genomic.fna.gz Y12.fa
```

**SK1**
```bash
mkdir ~/data/alignment/egaz/download/SK1
cd ~/data/alignment/egaz/download/SK1
aria2c -x 6 -s 3 -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/057/885/GCA_002057885.1_ASM205788v1/GCA_002057885.1_ASM205788v1_genomic.fna.gz
find . -name "*.gz" | xargs gzip -t
faops filter -N -s GCA_002057885.1_ASM205788v1_genomic.fna.gz SK1.fa
```

**UWOPS03-461.4**
```bash
mkdir ~/data/alignment/egaz/download/UWOPS03_461_4
cd ~/data/alignment/egaz/download/UWOPS03_461_4
aria2c -x 6 -s 3 -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/058/095/GCA_002058095.1_ASM205809v1/GCA_002058095.1_ASM205809v1_genomic.fna.gz
find . -name "*.gz" | xargs gzip -t
faops filter -N -s GCA_002058095.1_ASM205809v1_genomic.fna.gz UWOPS03_461_4.fa
```

**YPS128**
```bash
mkdir ~/data/alignment/egaz/download/YPS128
cd ~/data/alignment/egaz/download/YPS128
aria2c -x 6 -s 3 -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/057/995/GCA_002057995.1_ASM205799v1/GCA_002057995.1_ASM205799v1_genomic.fna.gz
find . -name "*.gz" | xargs gzip -t
faops filter -N -s GCA_002057995.1_ASM205799v1_genomic.fna.gz YPS128.fa
```

**DBVPG6765**
```bash
mkdir ~/data/alignment/egaz/download/DBVPG6765
cd ~/data/alignment/egaz/download/DBVPG6765
aria2c -x 6 -s 3 -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/057/805/GCA_002057805.1_ASM205780v1/GCA_002057805.1_ASM205780v1_genomic.fna.gz
find . -name "*.gz" | xargs gzip -t
faops filter -N -s GCA_002057805.1_ASM205780v1_genomic.fna.gz DBVPG6765.fa
```

**CBS432**
```bash
mkdir ~/data/alignment/egaz/download/CBS432
cd ~/data/alignment/egaz/download/CBS432
aria2c -x 6 -s 3 -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/079/055/GCA_002079055.1_ASM207905v1/GCA_002079055.1_ASM207905v1_genomic.fna.gz
find . -name "*.gz" | xargs gzip -t
faops filter -N -s GCA_002079055.1_ASM207905v1_genomic.fna.gz CBS432.fa
```

**N44**
```bash
mkdir ~/data/alignment/egaz/download/N44
cd ~/data/alignment/egaz/download/N44
aria2c -x 6 -s 3 -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/079/085/GCA_002079085.1_ASM207908v1/GCA_002079085.1_ASM207908v1_genomic.fna.gz
find . -name "*.gz" | xargs gzip -t
faops filter -N -s GCA_002079085.1_ASM207908v1_genomic.fna.gz N44.fa
```

**YPS138**
```bash
mkdir ~/data/alignment/egaz/download/YPS138
cd ~/data/alignment/egaz/download/YPS138
aria2c -x 6 -s 3 -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/079/115/GCA_002079115.1_ASM207911v1/GCA_002079115.1_ASM207911v1_genomic.fna.gz
find . -name "*.gz" | xargs gzip -t
faops filter -N -s GCA_002079115.1_ASM207911v1_genomic.fna.gz YPS138.fa
```
**UFRJ50816**
```bash
mkdir ~/data/alignment/egaz/download/UFRJ50816
cd ~/data/alignment/egaz/download/UFRJ50816
aria2c -x 6 -s 3 -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/079/145/GCA_002079145.1_ASM207914v1/GCA_002079145.1_ASM207914v1_genomic.fna.gz
find . -name "*.gz" | xargs gzip -t
faops filter -N -s GCA_002079145.1_ASM207914v1_genomic.fna.gz UFRJ50816.fa
```

**UWOPS91-917.1**
```bash
mkdir ~/data/alignment/egaz/download/UWOPS91_917_1
cd ~/data/alignment/egaz/download/UWOPS91_917_1
aria2c -x 6 -s 3 -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/079/175/GCA_002079175.1_ASM207917v1/GCA_002079175.1_ASM207917v1_genomic.fna.gz
find . -name "*.gz" | xargs gzip -t
faops filter -N -s GCA_002079175.1_ASM207917v1_genomic.fna.gz UWOPS91_917_1.fa
```

**EC1118**
```bash
mkdir ~/data/alignment/egaz/download/EC1118
cd ~/data/alignment/egaz/download/EC1118
aria2c -x 6 -s 3 -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/218/975/GCA_000218975.1_ASM21897v1/GCA_000218975.1_ASM21897v1_genomic.fna.gz
find . -name "*.gz" | xargs gzip -t
faops filter -N -s GCA_000218975.1_ASM21897v1_genomic.fna.gz EC1118.fa
```

**Seub**
```bash
mkdir ~/data/alignment/egaz/download/Seub
cd ~/data/alignment/egaz/download/Seub
aria2c -x 6 -s 3 -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/298/625/GCA_001298625.1_SEUB3.0/GCA_001298625.1_SEUB3.0_genomic.fna.gz
find . -name "*.gz" | xargs gzip -t
faops filter -N -s GCA_001298625.1_SEUB3.0_genomic.fna.gz Seub.fa
```

## Download strains from NCBI WGS

**scer**
```bash
cd ~/data/alignment/egaz/download
perl ~/Scripts/withncbi/taxon/wgs_prep.pl \
    -f ~/Scripts/pars/scer_wgs.tsv \
    --fix
bash scer_wgs.rsync.sh
```

**spar**
```bash
cd ~/data/alignment/egaz/download
perl ~/Scripts/withncbi/taxon/wgs_prep.pl \
    -f ~/Scripts/pars/spar_wgs.tsv \
    --fix
bash spar_wgs.rsync.sh
```

## Download strains from 1002genomes project

```bash
cd ~/data/alignment/egaz/download
wget -c http://1002genomes.u-strasbg.fr/files/1011Assemblies.tar.gz
tar -zxvf 1011Assemblies.tar.gz
```

# RepeatMasker

```bash
cd ~/data/alignment/egaz

egaz prepseq download/S288c/S288c.fa -o S288c -v
gzip -d -c download/S288c/Saccharomyces_cerevisiae.R64-1-1.82.gff3.gz > S288c/chr.gff
egaz masked S288c/*.fa -o S288c/repeat.yml

egaz prepseq \
    download/DBVPG6044/DBVPG6044.fa -o DBVPG6044 \
    --repeatmasker '--species Fungi --parallel 8' -v
egaz prepseq \
    download/Y12/Y12.fa -o Y12 \
    --repeatmasker '--species Fungi --parallel 8' -v
egaz prepseq \
    download/SK1/SK1.fa -o SK1 \
    --repeatmasker '--species Fungi --parallel 8' -v
egaz prepseq \
    download/UWOPS03_461_4/UWOPS03_461_4.fa -o UWOPS03_461_4 \
    --repeatmasker '--species Fungi --parallel 8' -v
egaz prepseq \
    download/YPS128/YPS128.fa -o YPS128 \
    --repeatmasker '--species Fungi --parallel 8' -v
egaz prepseq \
    download/DBVPG6765/DBVPG6765.fa -o DBVPG6765 \
    --repeatmasker '--species Fungi --parallel 8' -v
egaz prepseq \
    download/CBS432/CBS432.fa -o CBS432 \
    --repeatmasker '--species Fungi --parallel 8' -v
egaz prepseq \
    download/N44/N44.fa -o N44 \
    --repeatmasker '--species Fungi --parallel 8' -v
egaz prepseq \
    download/YPS138/YPS138.fa -o YPS138 \
    --repeatmasker '--species Fungi --parallel 8' -v
egaz prepseq \
    download/UFRJ50816/UFRJ50816.fa -o UFRJ50816 \
    --repeatmasker '--species Fungi --parallel 8' -v
egaz prepseq \
    download/UWOPS91_917_1/UWOPS91_917_1.fa -o UWOPS91_917_1 \
    --repeatmasker '--species Fungi --parallel 8' -v
egaz prepseq \
    download/EC1118/EC1118.fa -o EC1118 \
    --repeatmasker '--species Fungi --parallel 8' -v

mkdir -p ~/data/alignment/egaz/Seub
cd ~/data/alignment/egaz/Seub
RepeatMasker --species Fungi --parallel 16 -xsmall ../download/Seub/Seub.fa
egaz prepseq \
    ../download/Seub/Seub.fa.masked -v

cd ~/data/alignment/egaz
cat download/scer_wgs.csv \
    | grep -v "^prefix" \
    | cut -d',' -f1,3 \
    | uniq \
    | perl -nl -a -F"," -e 'printf qq{egaz prepseq \\\n   download/%s/%s.*.fsa_nt.gz -o %s \\\n   --about 2000000 --repeatmasker " --species Fungi --parallel 8" --min 1000 --gi -v \n}, $F[1], $F[0], $F[1];' > rm_scer_wgs.sh   
bash rm_scer_wgs.sh

cat download/spar_wgs.csv \
    | grep -v "^prefix" \
    | cut -d',' -f1,3 \
    | uniq \
    | perl -nl -a -F"," -e 'printf qq{egaz prepseq \\\n   download/%s/%s.*.fsa_nt.gz -o %s \\\n   --about 2000000 --repeatmasker " --species Fungi --parallel 8" --min 1000 --gi -v \n}, $F[1], $F[0], $F[1];' > rm_spar_wgs.sh
bash rm_spar_wgs.sh

# 1011
cd ~/data/alignment/egaz/download
for file in $(ls ~/data/alignment/egaz/download/GENOMES_ASSEMBLED/*.re.fa);
do

filename=$(basename $file)
dir=$(echo $filename | perl -p -e 's/^([A-Za-z]+).+/$1/;')

if [ -d ../$dir ];
then echo -n;
else
cd ~/data/alignment/egaz/download
cat $file | perl -nl -e '

if (m/^>([A-Za-z]+)/){

my $dir = $1;
mkdir ("../$dir") unless (-d "../$dir");
last;

}
'
cp -rf $file ../$dir
sed -i".bak" "s/-/_/" ../$dir/*.re.fa

faops filter -a 1000 ../$dir/*.re.fa ../$dir/$dir.fasta
rm -rf ../$dir/*.re.fa
RepeatMasker --species Fungi --parallel 16 -xsmall ../$dir/$dir.fasta

cd ~/data/alignment/egaz/$dir
egaz prepseq \
    ../$dir/$dir.fasta.masked -v
fi
done
```

# Plans of alignments

```bash
# create downloaded genome list
cat ~/Scripts/pars/scer_assembly.csv \
    | grep -v "^#" \
    | cut -d',' -f1,3 \
    | uniq \
    | perl -nl -a -F"," -e 'printf qq{    --download "name=%s;taxon=%s" \\\n}, $F[0], $F[1];'
    
cat ~/Scripts/pars/spar_assembly.csv \
    | grep -v "^#" \
    | cut -d',' -f1,3 \
    | uniq \
    | perl -nl -a -F"," -e 'printf qq{    --download "name=%s;taxon=%s" \\\n}, $F[0], $F[1];'    

mkdir -p ~/data/mrna-structure/alignment/spar_wgs
cd ~/data/mrna-structure/alignment/spar_wgs

perl ~/Scripts/withncbi/pop/gen_pop_conf.pl \
    -i ~/data/mrna-structure/GENOMES/WGS/spar_wgs.data.yml \
    -o spar_wgs.plan.yml \
    -d ~/data/mrna-structure/GENOMES/WGS \
    -m prefix \
    -r '*.fsa_nt.gz' \
    --opt group_name=spar_wgs \
    --opt base_dir='~/data/mrna-structure/alignment' \
    --opt data_dir="~/data/mrna-structure/alignment/spar_wgs" \
    --opt rm_species=Fungi \
    --dd ~/data/mrna-structure/GENOMES/ASSEMBLIES \
    --download "name=CBS432;taxon=100000000007" \
    --download "name=N44;taxon=100000000008" \
    --download "name=UWOPS91_917_1;taxon=100000000009" \
    --download "name=UFRJ50816;taxon=100000000010" \
    --download "name=YPS138;taxon=100000000011" \
    --plan 'name=Spar_n5_pop;t=CBS432;qs=N44,UWOPS91_917_1,UFRJ50816,YPS138' \
    -y

# pop_prep.pl
perl ~/Scripts/withncbi/pop/pop_prep.pl -p 16 -i spar_wgs.plan.yml

bash 01_file.sh
bash 02_rm.sh
bash 03_strain_info.sh

# plan_ALL.sh
bash plan_ALL.sh

bash 1_real_chr.sh
bash 3_pair_cmd.sh
bash 4_rawphylo.sh
bash 5_multi_cmd.sh

# other plans
bash plan_Spar_n5_pop.sh
bash 5_multi_cmd.sh

#consensus
mkdir -p ~/data/mrna-structure/alignment/spar_wgs/consensus
cd ~/data/mrna-structure/alignment/spar_wgs/consensus
cp -rf ~/data/mrna-structure/alignment/spar_wgs/Spar_n5_pop_refined ~/data/mrna-structure/alignment/spar_wgs/consensus
gunzip -rfvc Spar_n5_pop_refined/*.maf.gz.fas.gz > Spar_n5_pop_refined/Spar.fas
fasops consensus Spar_n5_pop_refined/Spar.fas -o Spar_n5_pop_refined/consensus.fas -p 2

mkdir -p ~/data/mrna-structure/GENOMES/ASSEMBLIES/consensus
cd ~/data/mrna-structure/GENOMES/ASSEMBLIES/consensus
cp -rf ~/data/mrna-structure/alignment/spar_wgs/consensus/Spar_n5_pop_refined/consensus.fas ~/data/mrna-structure/GENOMES/ASSEMBLIES/consensus

perl -p -i.bak -e 's/^(\w+)\n\n/$1\n/g','s/\A\s*\Z//' consensus.fas
faops filter -a 1000 consensus.fas consensus.fasta
perl ~/Scripts/pars/program/rename.pl consensus.fasta
rm -rf consensus.fas.bak consensus.fasta.bak consensus.fas

mkdir -p ~/data/mrna-structure/alignment/scer_wgs
cd ~/data/mrna-structure/alignment/scer_wgs

perl ~/Scripts/withncbi/pop/gen_pop_conf.pl \
    -i ~/data/mrna-structure/GENOMES/WGS/scer_wgs.data.yml \
    -o scer_wgs.plan.yml \
    -d ~/data/mrna-structure/GENOMES/WGS \
    -m prefix \
    -r '*.fsa_nt.gz' \
    --opt group_name=scer_wgs \
    --opt base_dir='~/data/mrna-structure/alignment' \
    --opt data_dir="~/data/mrna-structure/alignment/scer_wgs" \
    --opt rm_species=Fungi \
    --dd ~/data/mrna-structure/GENOMES/ASSEMBLIES \
    --download "name=S288c;taxon=559292" \
    --download "name=DBVPG6044;taxon=100000000001" \
    --download "name=UWOPS03_461_4;taxon=100000000002" \
    --download "name=Y12;taxon=100000000003" \
    --download "name=SK1;taxon=100000000004" \
    --download "name=YPS128;taxon=100000000005" \
    --download "name=DBVPG6765;taxon=100000000006" \
    --download "name=EC1118;taxon=643680" \
    --download "name=consensus;taxon=100000000012" \
    --plan 'name=Scer_n7_Spar;t=S288c;qs=EC1118,Kyokai_no_7,RM11_1a,Sigma1278b,T7,YJM789,Spar;o=Spar' \
    --plan 'name=Scer_n7p_Spar;t=S288c;qs=DBVPG6044,UWOPS03_461_4,Y12,SK1,YPS128,DBVPG6765,Spar;o=Spar' \
    --plan 'name=Scer_n157_Spar;t=S288c;qs=beer001,beer002,beer003,beer004,beer005,beer006,beer007,beer008,beer009,beer010,beer011,beer012,beer013,beer014,beer015,beer016,beer017,beer018,beer019,beer020,beer021,beer022,beer023,beer024,beer025,beer026,beer027,beer028,beer029,beer030,beer031,beer032,beer033,beer034,beer035,beer036,beer037,beer038,beer039,beer040,beer041,beer042,beer043,beer044,beer045,beer046,beer047,beer048,beer049,beer050,beer051,beer052,beer053,beer054,beer055,beer056,beer057,beer058,beer059,beer060,beer061,beer062,beer063,beer064,beer065,beer066,beer067,beer068,beer069,beer070,beer071,beer072,beer073,beer074,beer075,beer076,beer077,beer078,beer079,beer080,beer081,beer082,beer083,beer084,beer085,beer086,beer087,beer088,beer089,beer090,beer091,beer092,beer093,beer094,beer095,beer096,beer097,beer098,beer099,beer100,beer101,beer102,bioethanol001,bioethanol002,bioethanol003,bioethanol004,bioethanol005,bread001,bread002,bread003,bread004,laboratory001,laboratory002,sake001,sake002,sake003,sake004,sake005,sake006,sake007,spirits001,spirits002,spirits003,spirits004,spirits005,spirits006,spirits007,spirits008,spirits009,spirits010,spirits011,wine001,wine002,wine003,wine004,wine005,wine006,wine007,wine008,wine009,wine010,wine011,wine012,wine013,wine014,wine015,wine016,wine017,wine018,wine019,wild001,wild002,wild003,wild004,wild005,wild006,wild007,Spar;o=Spar' \
    --plan 'name=Scer_n157_nonMosaic_Spar;t=S288c;qs=beer001,beer002,beer003,beer004,beer005,beer006,beer007,beer008,beer009,beer010,beer011,beer012,beer013,beer014,beer015,beer016,beer020,beer021,beer022,beer023,beer024,beer025,beer026,beer027,beer028,beer029,beer030,beer031,beer032,beer033,beer034,beer036,beer037,beer038,beer039,beer040,beer041,beer043,beer044,beer045,beer046,beer047,beer048,beer049,beer050,beer051,beer052,beer053,beer054,beer055,beer056,beer059,beer061,beer062,beer063,beer064,beer065,beer066,beer067,beer068,beer069,beer070,beer071,beer073,beer075,beer076,beer077,beer078,beer079,beer080,beer081,beer082,beer083,beer084,beer085,beer086,beer087,beer088,beer089,beer090,beer091,beer092,beer094,beer095,beer096,beer097,beer098,beer099,beer100,beer101,beer102,bioethanol001,bioethanol003,bioethanol004,bread001,bread002,bread003,bread004,sake001,sake002,sake003,sake004,sake005,sake006,sake007,spirits001,spirits002,spirits003,spirits004,spirits005,spirits011,wine001,wine003,wine004,wine005,wine006,wine007,wine009,wine010,wine011,wine012,wine013,wine014,wine015,wine017,wine018,wild004,wild005,wild006,wild007,Spar;o=Spar' \
    --plan 'name=Scer_n157_nonMosaic_consensus;t=S288c;qs=beer001,beer002,beer003,beer004,beer005,beer006,beer007,beer008,beer009,beer010,beer011,beer012,beer013,beer014,beer015,beer016,beer020,beer021,beer022,beer023,beer024,beer025,beer026,beer027,beer028,beer029,beer030,beer031,beer032,beer033,beer034,beer036,beer037,beer038,beer039,beer040,beer041,beer043,beer044,beer045,beer046,beer047,beer048,beer049,beer050,beer051,beer052,beer053,beer054,beer055,beer056,beer059,beer061,beer062,beer063,beer064,beer065,beer066,beer067,beer068,beer069,beer070,beer071,beer073,beer075,beer076,beer077,beer078,beer079,beer080,beer081,beer082,beer083,beer084,beer085,beer086,beer087,beer088,beer089,beer090,beer091,beer092,beer094,beer095,beer096,beer097,beer098,beer099,beer100,beer101,beer102,bioethanol001,bioethanol003,bioethanol004,bread001,bread002,bread003,bread004,sake001,sake002,sake003,sake004,sake005,sake006,sake007,spirits001,spirits002,spirits003,spirits004,spirits005,spirits011,wine001,wine003,wine004,wine005,wine006,wine007,wine009,wine010,wine011,wine012,wine013,wine014,wine015,wine017,wine018,wild004,wild005,wild006,wild007,consensus;o=consensus' \
    --plan 'name=Scer_n128_Spar;t=S288c;qs=beer001,beer003,beer004,beer005,beer006,beer007,beer008,beer009,beer010,beer011,beer012,beer013,beer014,beer015,beer016,beer020,beer021,beer022,beer023,beer024,beer025,beer026,beer027,beer028,beer029,beer030,beer031,beer032,beer033,beer034,beer036,beer037,beer038,beer040,beer041,beer043,beer044,beer045,beer046,beer047,beer048,beer049,beer050,beer051,beer052,beer053,beer054,beer055,beer056,beer059,beer061,beer062,beer063,beer064,beer065,beer066,beer067,beer068,beer069,beer070,beer071,beer073,beer075,beer076,beer077,beer078,beer079,beer080,beer081,beer082,beer083,beer084,beer085,beer086,beer087,beer088,beer089,beer090,beer091,beer092,beer094,beer095,beer096,beer097,beer098,beer099,beer100,beer101,beer102,bioethanol001,bioethanol003,bioethanol004,bread001,bread002,bread003,bread004,sake001,sake002,sake003,sake004,sake005,sake006,sake007,spirits001,spirits002,spirits003,spirits004,spirits005,spirits011,wine001,wine003,wine004,wine005,wine006,wine007,wine009,wine010,wine011,wine012,wine013,wine014,wine015,wine017,wine018,wild005,wild006,wild007,Spar;o=Spar' \
    -y

# pop_prep.pl
perl ~/Scripts/withncbi/pop/pop_prep.pl -p 16 -i scer_wgs.plan.yml

bash 01_file.sh
bash 02_rm.sh
bash 03_strain_info.sh

# plan_ALL.sh
bash plan_ALL.sh

bash 1_real_chr.sh
bash 3_pair_cmd.sh
bash 4_rawphylo.sh
bash 5_multi_cmd.sh
bash 7_multi_db_only.sh

# other plans
bash plan_Scer_n7_Spar.sh
bash 5_multi_cmd.sh
bash 7_multi_db_only.sh

# other plans
bash plan_Scer_n7p_Spar.sh
bash 5_multi_cmd.sh
bash 7_multi_db_only.sh

# other plans
bash plan_Scer_n157_Spar.sh
bash 5_multi_cmd.sh
bash 7_multi_db_only.sh

# other plans
bash plan_Scer_n157_nonMosaic_Spar.sh
bash 5_multi_cmd.sh
bash 7_multi_db_only.sh

# other plans
bash plan_Scer_n157_nonMosaic_consensus.sh
bash 5_multi_cmd.sh
bash 7_multi_db_only.sh

# other plans
bash plan_Scer_n128_Spar.sh
bash 5_multi_cmd.sh
bash 7_multi_db_only.sh

```

# AlignDB

## Build alignDB for multiple genomes n7

```bash
mkdir -p ~/data/mrna-structure/xlsx
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/alignDB/alignDB.pl \
    -d Scer_n7_Spar \
    -da ~/data/mrna-structure/alignment/scer_wgs/Scer_n7_Spar_refined \
    -a ~/data/mrna-structure/alignment/scer_wgs/Stats/anno.yml\
    --ensembl saccharomyces_cerevisiae_core_29_82_4 \
    --outgroup \
    --chr ~/data/mrna-structure/alignment/scer_wgs/chr_length.csv \
    -lt 1000 --parallel 8 --batch 5 \
    --run gene

perl ~/Scripts/alignDB/stat/mvar_stat_factory.pl \
    -d Scer_n7_Spar -r 1-60

```

## Build alignDB for multiple genomes n7p

```bash
mkdir -p ~/data/mrna-structure/xlsx
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/alignDB/alignDB.pl \
    -d Scer_n7p_Spar \
    -da ~/data/mrna-structure/alignment/scer_wgs/Scer_n7p_Spar_refined \
    -a ~/data/mrna-structure/alignment/scer_wgs/Stats/anno.yml\
    --ensembl saccharomyces_cerevisiae_core_29_82_4 \
    --outgroup \
    --chr ~/data/mrna-structure/alignment/scer_wgs/chr_length.csv \
    -lt 1000 --parallel 8 --batch 5 \
    --run gene

perl ~/Scripts/alignDB/stat/mvar_stat_factory.pl \
    -d Scer_n7p_Spar -r 1-60

```

## Build alignDB for multiple genomes n157

```bash
mkdir -p ~/data/mrna-structure/xlsx
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/alignDB/alignDB.pl \
    -d Scer_n157_Spar \
    -da ~/data/mrna-structure/alignment/scer_wgs/Scer_n157_Spar_refined \
    -a ~/data/mrna-structure/alignment/scer_wgs/Stats/anno.yml\
    --ensembl saccharomyces_cerevisiae_core_29_82_4 \
    --outgroup \
    --chr ~/data/mrna-structure/alignment/scer_wgs/chr_length.csv \
    -lt 1000 --parallel 8 --batch 5 \
    --run gene

perl ~/Scripts/alignDB/stat/mvar_stat_factory.pl \
    -d Scer_n157_Spar -r 1-60

```

## Build alignDB for multiple genomes n157_nonMosaic

```bash
mkdir -p ~/data/mrna-structure/xlsx
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/alignDB/alignDB.pl \
    -d Scer_n157_nonMosaic_Spar \
    -da ~/data/mrna-structure/alignment/scer_wgs/Scer_n157_nonMosaic_Spar_refined \
    -a ~/data/mrna-structure/alignment/scer_wgs/Stats/anno.yml\
    --ensembl saccharomyces_cerevisiae_core_29_82_4 \
    --outgroup \
    --chr ~/data/mrna-structure/alignment/scer_wgs/chr_length.csv \
    -lt 1000 --parallel 8 --batch 5 \
    --run gene

perl ~/Scripts/alignDB/stat/mvar_stat_factory.pl \
    -d Scer_n157_nonMosaic_Spar -r 1-60

```

## Build alignDB for multiple genomes n157_nonMosaic_consensus

```bash
mkdir -p ~/data/mrna-structure/xlsx
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/alignDB/alignDB.pl \
    -d Scer_n157_nonMosaic_consensus \
    -da ~/data/mrna-structure/alignment/scer_wgs/multi160/Scer_n157_nonMosaic_consensus_refined \
    -a ~/data/mrna-structure/alignment/scer_wgs/multi160/Results/anno.yml\
    --ensembl saccharomyces_cerevisiae_core_29_82_4 \
    --outgroup \
    --chr ~/data/mrna-structure/alignment/scer_wgs/multi160/Results/chr_length.csv \
    -lt 1000 --parallel 8 --batch 5 \
    --run gene

perl ~/Scripts/alignDB/stat/mvar_stat_factory.pl \
    -d Scer_n157_nonMosaic_consensus -r 1-60

```

## Build alignDB for multiple genomes n128_Spar

```bash
mkdir -p ~/data/mrna-structure/xlsx
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/alignDB/alignDB.pl \
    -d Scer_n128_Spar \
    -da ~/data/mrna-structure/alignment/scer_wgs/multi160/Scer_n128_Spar_refined \
    -a ~/data/mrna-structure/alignment/scer_wgs/multi160/Results/anno.yml\
    --ensembl saccharomyces_cerevisiae_core_29_82_4 \
    --outgroup \
    --chr ~/data/mrna-structure/alignment/scer_wgs/multi160/Results/chr_length.csv \
    -lt 1000 --parallel 8 --batch 5 \
    --run gene

perl ~/Scripts/alignDB/stat/mvar_stat_factory.pl \
    -d Scer_n128_Spar -r 1-60

```

## Build alignDB for multiple genomes n128_consensus

```bash
mkdir -p ~/data/mrna-structure/xlsx
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/alignDB/alignDB.pl \
    -d Scer_n128_consensus \
    -da ~/data/mrna-structure/alignment/scer_wgs/multi160/Scer_n128_consensus_refined \
    -a ~/data/mrna-structure/alignment/scer_wgs/multi160/Results/anno.yml\
    --ensembl saccharomyces_cerevisiae_core_29_82_4 \
    --outgroup \
    --chr ~/data/mrna-structure/alignment/scer_wgs/multi160/Results/chr_length.csv \
    -lt 1000 --parallel 8 --batch 5 \
    --run gene

perl ~/Scripts/alignDB/stat/mvar_stat_factory.pl \
    -d Scer_n128_consensus -r 1-60
  
```

## Extract gene-list and snp-codon-list n7

```bash
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n7_Spar.mvar.1-60.xlsx --sheet 'gene_list' \
    > Scer_n7_Spar.mvar.gene_list.csv

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n7_Spar.mvar.1-60.xlsx --sheet 'snp_codon_list' \
    > Scer_n7_Spar.mvar.gene_list.csv
    
```

## Extract gene-list and snp-codon-list n7p

```bash
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n7p_Spar.mvar.1-60.xlsx --sheet 'gene_list' \
    > Scer_n7p_Spar.mvar.gene_list.csv

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n7p_Spar.mvar.1-60.xlsx --sheet 'snp_codon_list' \
    > Scer_n7p_Spar.mvar.gene_list.csv
    
```

## Extract gene-list and snp-codon-list n157

```bash
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n157_Spar.mvar.1-60.xlsx --sheet 'gene_list' \
    > Scer_n157_Spar.mvar.gene_list.csv

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n157_Spar.mvar.1-60.xlsx --sheet 'snp_codon_list' \
    > Scer_n157_Spar.mvar.gene_list.csv
    
```

## Extract gene-list and snp-codon-list n157_nonMosaic

```bash
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n157_nonMosaic_Spar.mvar.1-60.xlsx --sheet 'gene_list' \
    > Scer_n157_nonMosaic_Spar.mvar.gene_list.csv

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n157_nonMosaic_Spar.mvar.1-60.xlsx --sheet 'snp_codon_list' \
    > Scer_n157_nonMosaic_Spar.mvar.gene_list.csv
    
```

## Extract gene-list and snp-codon-list n157_nonMosaic_consensus

```bash
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n157_nonMosaic_consensus.mvar.1-60.xlsx --sheet 'gene_list' \
    > Scer_n157_nonMosaic_consensus.mvar.gene_list.csv

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n157_nonMosaic_consensus.mvar.1-60.xlsx --sheet 'snp_codon_list' \
    > Scer_n157_nonMosaic_consensus.mvar.gene_list.csv
    
```

## Extract gene-list and snp-codon-list n128

```bash
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n128_Spar.mvar.1-60.xlsx --sheet 'gene_list' \
    > Scer_n128_Spar.mvar.gene_list.csv

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n128_Spar.mvar.1-60.xlsx --sheet 'snp_codon_list' \
    > Scer_n128_Spar.mvar.gene_list.csv
    
```

## SNPs and indels n7

Select columns `chr_name,snp_pos` for SNPs.

Select columns `chr_name,indel_start,indel_end` for indels.

```bash
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n7_Spar.mvar.1-60.xlsx --sheet 'snp_list' \
    | perl -nla -F"," -e '
        /^\d/ or next;
        print qq{$F[2]:$F[3]};
    ' \
    > Scer_n7_Spar.snp.pos.txt

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n7_Spar.mvar.1-60.xlsx --sheet 'indel_list' \
    | perl -nla -F"," -e '
        /^\d/ or next;
        if ( $F[3] == $F[4] ) {
            print qq{$F[2]:$F[3]};
        }
        else {
            print qq{$F[2]:$F[3]-$F[4]};
        }
    ' \
    > Scer_n7_Spar.indel.pos.txt

```

## SNPs and indels n7p

Select columns `chr_name,snp_pos` for SNPs.

Select columns `chr_name,indel_start,indel_end` for indels.

```bash
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n7p_Spar.mvar.1-60.xlsx --sheet 'snp_list' \
    | perl -nla -F"," -e '
        /^\d/ or next;
        print qq{$F[2]:$F[3]};
    ' \
    > Scer_n7p_Spar.snp.pos.txt

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n7p_Spar.mvar.1-60.xlsx --sheet 'indel_list' \
    | perl -nla -F"," -e '
        /^\d/ or next;
        if ( $F[3] == $F[4] ) {
            print qq{$F[2]:$F[3]};
        }
        else {
            print qq{$F[2]:$F[3]-$F[4]};
        }
    ' \
    > Scer_n7p_Spar.indel.pos.txt

```

## SNPs and indels n157

Select columns `chr_name,snp_pos` for SNPs.

Select columns `chr_name,indel_start,indel_end` for indels.

```bash
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n157_Spar.mvar.1-60.xlsx --sheet 'snp_list' \
    | perl -nla -F"," -e '
        /^\d/ or next;
        print qq{$F[2]:$F[3]};
    ' \
    > Scer_n157_Spar.snp.pos.txt

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n157_Spar.mvar.1-60.xlsx --sheet 'indel_list' \
    | perl -nla -F"," -e '
        /^\d/ or next;
        if ( $F[3] == $F[4] ) {
            print qq{$F[2]:$F[3]};
        }
        else {
            print qq{$F[2]:$F[3]-$F[4]};
        }
    ' \
    > Scer_n157_Spar.indel.pos.txt

```

## SNPs and indels n157_nonMosaic

Select columns `chr_name,snp_pos` for SNPs.

Select columns `chr_name,indel_start,indel_end` for indels.

```bash
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n157_nonMosaic_Spar.mvar.1-60.xlsx --sheet 'snp_list' \
    | perl -nla -F"," -e '
        /^\d/ or next;
        print qq{$F[2]:$F[3]};
    ' \
    > Scer_n157_nonMosaic_Spar.snp.pos.txt

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n157_nonMosaic_Spar.mvar.1-60.xlsx --sheet 'indel_list' \
    | perl -nla -F"," -e '
        /^\d/ or next;
        if ( $F[3] == $F[4] ) {
            print qq{$F[2]:$F[3]};
        }
        else {
            print qq{$F[2]:$F[3]-$F[4]};
        }
    ' \
    > Scer_n157_nonMosaic_Spar.indel.pos.txt

```

## SNPs and indels n157_nonMosaic_consensus

Select columns `chr_name,snp_pos` for SNPs.

Select columns `chr_name,indel_start,indel_end` for indels.

```bash
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n157_nonMosaic_consensus.mvar.1-60.xlsx --sheet 'snp_list' \
    | perl -nla -F"," -e '
        /^\d/ or next;
        print qq{$F[2]:$F[3]};
    ' \
    > Scer_n157_nonMosaic_consensus.snp.pos.txt

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n157_nonMosaic_consensus.mvar.1-60.xlsx --sheet 'indel_list' \
    | perl -nla -F"," -e '
        /^\d/ or next;
        if ( $F[3] == $F[4] ) {
            print qq{$F[2]:$F[3]};
        }
        else {
            print qq{$F[2]:$F[3]-$F[4]};
        }
    ' \
    > Scer_n157_nonMosaic_consensus.indel.pos.txt

```

## SNPs and indels n128

Select columns `chr_name,snp_pos` for SNPs.

Select columns `chr_name,indel_start,indel_end` for indels.

```bash
cd ~/data/mrna-structure/xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n128_Spar.mvar.1-60.xlsx --sheet 'snp_list' \
    | perl -nla -F"," -e '
        /^\d/ or next;
        print qq{$F[2]:$F[3]};
    ' \
    > Scer_n128_Spar.snp.pos.txt

perl ~/Scripts/fig_table/xlsx2csv.pl -f Scer_n128_Spar.mvar.1-60.xlsx --sheet 'indel_list' \
    | perl -nla -F"," -e '
        /^\d/ or next;
        if ( $F[3] == $F[4] ) {
            print qq{$F[2]:$F[3]};
        }
        else {
            print qq{$F[2]:$F[3]-$F[4]};
        }
    ' \
    > Scer_n128_Spar.indel.pos.txt

```

# Blast

Prepare a combined fasta file of yeast genome and blast genes against
the genome.

```bash
mkdir -p ~/data/mrna-structure/blast
cd ~/data/mrna-structure/blast

cat ~/data/mrna-structure/alignment/scer_wgs/Genomes/S288c/{I,II,III,IV,V,VI,VII,VIII,IX,X,XI,XII,XIII,XIV,XV,XVI,Mito}.fa \
    > S288c.fa

perl -nl -i -e '/^>/ or $_ = uc $_; print'  S288c.fa
faops size S288c.fa > S288c.sizes

# formatdb
~/share/blast/bin/formatdb -p F -o T -i S288c.fa

# blast every transcripts against genome
~/share/blast/bin/blastall -p blastn -F "m D" -m 0 -b 10 -v 10 -e 1e-3 -a 4 \
    -i ../PARS10/pubs/PARS10/data/sce_genes.fasta -d S288C.fa -o sce_genes.blast
```

# Features

```bash
mkdir -p ~/data/mrna-structure/process
cd ~/data/mrna-structure/process

#----------------------------------------------------------#
# gene
#----------------------------------------------------------#
# parse blastn output
perl ~/Scripts/pars/blastn_transcript.pl -f ../blast/sce_genes.blast -m 0

# produce transcript set
# YLR167W	568	chrXII	498888	499455	+
cat sce_genes.blast.tsv \
    | perl -nla -e 'print qq{$F[2]:$F[3]-$F[4]}' \
    | sort \
    > sce_genes.pos.txt
jrunlist cover sce_genes.pos.txt -o sce_genes.yml

#----------------------------------------------------------#
# intergenic
#----------------------------------------------------------#
cat ../sgd/NotFeature.fasta \
    | perl -n -e '
        />/ or next;
        /Chr\s+(\w+)\s+from\s+(\d+)\-(\d+)/ or next;
        $1 eq "Mito" and next;
        print qq{$1:$2-$3\n};
    ' \
    > sce_intergenic.pos.txt
jrunlist cover sce_intergenic.pos.txt -o sce_intergenic.yml

#----------------------------------------------------------#
# intron
#----------------------------------------------------------#
cat ../sgd/orf_coding_all.fasta \
    | perl -n -MAlignDB::IntSpan -e '
        />/ or next;
        /Chr\s+(\w+)\s+from\s+([\d,-]+)/ or next;
        $1 eq "Mito" and next;

        my $chr = $1;
        my $range = $2;
        my @ranges = sort { $a <=> $b } grep {/^\d+$/} split /,|\-/, $range;
        my $intspan = AlignDB::IntSpan->new()->add_range(@ranges);
        my $hole = $intspan->holes;
        
        printf qq{%s:%s\n}, $chr, $hole->as_string if $hole->is_not_empty;
    ' \
    > sce_intron.pos.txt
jrunlist cover sce_intron.pos.txt -o sce_intron.yml

#----------------------------------------------------------#
# utr
#----------------------------------------------------------#
# produce orf_genomic set
cat ../sgd/orf_genomic_all.fasta \
    | perl -n -e '
        />/ or next;
        /Chr\s+(\w+)\s+from\s+(\d+)\-(\d+)/ or next;
        $1 eq "Mito" and next;

        if ($2 == $3) {
            print qq{$1:$2\n};
        }
        elsif ($2 < $3) {
            print qq{$1:$2-$3\n};
        }
        else {
            print qq{$1:$3-$2\n};
        }
    ' \
    > sce_orf_genomic.pos.txt
jrunlist cover sce_orf_genomic.pos.txt -o sce_orf_genomic.yml

jrunlist compare --op diff sce_genes.yml sce_orf_genomic.yml -o sce_utr.yml
runlist convert sce_utr.yml -o sce_utr.pos.txt

jrunlist compare --op diff sce_genes.yml sce_intron.yml -o sce_mRNA.yml
runlist convert sce_mRNA.yml -o sce_mRNA.pos.txt

jrunlist compare --op diff sce_mRNA.yml sce_utr.yml -o sce_cds.yml
runlist convert sce_cds.yml -o sce_cds.pos.txt

# Stats
printf "| %s | %s | %s | %s |\n" \
    "Name" "chrLength" "size" "coverage" \
    > coverage.stat.md
printf "|:--|--:|--:|--:|\n" >> coverage.stat.md

for f in genes intergenic intron orf_genomic utr mRNA cds; do
    printf "| %s | %s | %s | %s |\n" \
        ${f} \
        $( 
            jrunlist stat ../blast/S288c.sizes sce_${f}.yml --all -o stdout \
            | grep -v coverage \
            | sed "s/,/ /g"
        )
done >> coverage.stat.md

cat coverage.stat.md
```
| Name | chrLength | size | coverage |
|:--|--:|--:|--:|
| genes | 12071326 | 4235405 | 0.3509 |
| intergenic | 12071326 | 2864170 | 0.2373 |
| intron | 12071326 | 65144 | 0.0054 |
| orf_genomic | 12071326 | 8895737 | 0.7369 |
| utr | 12071326 | 516569 | 0.0428 |
| mRNA | 12071326 | 4233361 | 0.3507 |
| cds | 12071326 | 3716792 | 0.3079 |


# Real Processing n7

```bash
export NAME=Scer_n7_Spar

cd ~/data/mrna-structure/process

# SNPs within transcripts
runlist position --op superset \
    sce_genes.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.gene.pos.txt

# read gene and snp info file
# produce ${NAME}.gene_variation.yml
perl ~/Scripts/pars/read_fold.pl \
    --pars ../PARS10/pubs/PARS10/data \
    --gene sce_genes.blast.tsv \
    --pos  ${NAME}.snp.gene.pos.txt \
    > fail_pos.txt

# review fail_pos.txt to find SNPs located in overlapped genes

# process ${NAME}.gene_variation.yml
perl ~/Scripts/pars/process_vars_in_fold.pl --file ${NAME}.gene_variation.yml

# SNPs within orf_genomic regions
runlist position --op superset \
    sce_orf_genomic.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.orf_genomic.pos.txt

# SNPs within intergenic regions
runlist position --op superset \
    sce_intergenic.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.intergenic.pos.txt

# SNPs within introns
runlist position --op superset \
    sce_intron.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.intron.pos.txt

# SNPs within utr
runlist position --op superset \
    sce_utr.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.utr.pos.txt

# SNPs within cds
runlist position --op superset \
    sce_cds.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.cds.pos.txt
unset NAME

```

# Real Processing n7p

```bash
export NAME=Scer_n7p_Spar

cd ~/data/mrna-structure/process

# SNPs within transcripts
runlist position --op superset \
    sce_genes.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.gene.pos.txt

# read gene and snp info file
# produce ${NAME}.gene_variation.yml
perl ~/Scripts/pars/read_fold.pl \
    --pars ../PARS10/pubs/PARS10/data \
    --gene sce_genes.blast.tsv \
    --pos  ${NAME}.snp.gene.pos.txt \
    > fail_pos.txt

# review fail_pos.txt to find SNPs located in overlapped genes

# process ${NAME}.gene_variation.yml
perl ~/Scripts/pars/process_vars_in_fold.pl --file ${NAME}.gene_variation.yml

# SNPs within orf_genomic regions
runlist position --op superset \
    sce_orf_genomic.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.orf_genomic.pos.txt

# SNPs within intergenic regions
runlist position --op superset \
    sce_intergenic.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.intergenic.pos.txt

# SNPs within introns
runlist position --op superset \
    sce_intron.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.intron.pos.txt

# SNPs within utr
runlist position --op superset \
    sce_utr.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.utr.pos.txt

# SNPs within cds
runlist position --op superset \
    sce_cds.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.cds.pos.txt
unset NAME

```

# Real Processing n157

```bash
export NAME=Scer_n157_Spar

cd ~/data/mrna-structure/process

# SNPs within transcripts
runlist position --op superset \
    sce_genes.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.gene.pos.txt

# read gene and snp info file
# produce ${NAME}.gene_variation.yml
perl ~/Scripts/pars/read_fold.pl \
    --pars ../PARS10/pubs/PARS10/data \
    --gene sce_genes.blast.tsv \
    --pos  ${NAME}.snp.gene.pos.txt \
    > fail_pos.txt

# review fail_pos.txt to find SNPs located in overlapped genes

# process ${NAME}.gene_variation.yml
perl ~/Scripts/pars/process_vars_in_fold.pl --file ${NAME}.gene_variation.yml

# SNPs within orf_genomic regions
runlist position --op superset \
    sce_orf_genomic.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.orf_genomic.pos.txt

# SNPs within intergenic regions
runlist position --op superset \
    sce_intergenic.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.intergenic.pos.txt

# SNPs within introns
runlist position --op superset \
    sce_intron.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.intron.pos.txt

# SNPs within utr
runlist position --op superset \
    sce_utr.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.utr.pos.txt

# SNPs within cds
runlist position --op superset \
    sce_cds.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.cds.pos.txt
unset NAME

```

# Real Processing n157_nonMosaic

```bash
export NAME=Scer_n157_nonMosaic_Spar

cd ~/data/mrna-structure/process

# SNPs within transcripts
runlist position --op superset \
    sce_genes.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.gene.pos.txt

# read gene and snp info file
# produce ${NAME}.gene_variation.yml
perl ~/Scripts/pars/read_fold.pl \
    --pars ../PARS10/pubs/PARS10/data \
    --gene sce_genes.blast.tsv \
    --pos  ${NAME}.snp.gene.pos.txt \
    > fail_pos.txt

# review fail_pos.txt to find SNPs located in overlapped genes

# process ${NAME}.gene_variation.yml
perl ~/Scripts/pars/process_vars_in_fold.pl --file ${NAME}.gene_variation.yml

# SNPs within orf_genomic regions
runlist position --op superset \
    sce_orf_genomic.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.orf_genomic.pos.txt

# SNPs within intergenic regions
runlist position --op superset \
    sce_intergenic.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.intergenic.pos.txt

# SNPs within introns
runlist position --op superset \
    sce_intron.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.intron.pos.txt

# SNPs within utr
runlist position --op superset \
    sce_utr.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.utr.pos.txt

# SNPs within cds
runlist position --op superset \
    sce_cds.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.cds.pos.txt
unset NAME

```

# Real Processing n157_nonMosaic_consensus

```bash
export NAME=Scer_n157_nonMosaic_consensus

cd ~/data/mrna-structure/process

# SNPs within transcripts
runlist position --op superset \
    sce_genes.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.gene.pos.txt

# read gene and snp info file
# produce ${NAME}.gene_variation.yml
perl ~/Scripts/pars/read_fold.pl \
    --pars ../PARS10/pubs/PARS10/data \
    --gene sce_genes.blast.tsv \
    --pos  ${NAME}.snp.gene.pos.txt \
    > fail_pos.txt

# review fail_pos.txt to find SNPs located in overlapped genes

# process ${NAME}.gene_variation.yml
perl ~/Scripts/pars/process_vars_in_fold.pl --file ${NAME}.gene_variation.yml

# SNPs within orf_genomic regions
runlist position --op superset \
    sce_orf_genomic.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.orf_genomic.pos.txt

# SNPs within intergenic regions
runlist position --op superset \
    sce_intergenic.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.intergenic.pos.txt

# SNPs within introns
runlist position --op superset \
    sce_intron.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.intron.pos.txt

# SNPs within utr
runlist position --op superset \
    sce_utr.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.utr.pos.txt

# SNPs within cds
runlist position --op superset \
    sce_cds.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.cds.pos.txt
unset NAME

```

# Real Processing n128

```bash
export NAME=Scer_n128_Spar

cd ~/data/mrna-structure/process

# SNPs within transcripts
runlist position --op superset \
    sce_genes.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.gene.pos.txt

# read gene and snp info file
# produce ${NAME}.gene_variation.yml
perl ~/Scripts/pars/read_fold.pl \
    --pars ../PARS10/pubs/PARS10/data \
    --gene sce_genes.blast.tsv \
    --pos  ${NAME}.snp.gene.pos.txt \
    > fail_pos.txt

# review fail_pos.txt to find SNPs located in overlapped genes

# process ${NAME}.gene_variation.yml
perl ~/Scripts/pars/process_vars_in_fold.pl --file ${NAME}.gene_variation.yml

# SNPs within orf_genomic regions
runlist position --op superset \
    sce_orf_genomic.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.orf_genomic.pos.txt

# SNPs within intergenic regions
runlist position --op superset \
    sce_intergenic.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.intergenic.pos.txt

# SNPs within introns
runlist position --op superset \
    sce_intron.yml ../xlsx/${NAME}.snp.pos.txt \
    -o ${NAME}.snp.intron.pos.txt

# SNPs within utr
runlist position --op superset \
    sce_utr.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.utr.pos.txt

# SNPs within cds
runlist position --op superset \
    sce_cds.yml ${NAME}.snp.gene.pos.txt \
    -o ${NAME}.snp.cds.pos.txt
unset NAME

```

# Download other reference data

## mRNA levels

* Data source: Quantification of the yeast transcriptome by
  single-molecule sequencing. Lipson, D. et al. Nature Biotechnology 27,
  652-658 (2009) doi:10.1038/nbt.1551

```bash
mkdir -p ~/data/mrna-structure/gene-traits
cd ~/data/mrna-structure/gene-traits

wget -N http://www.nature.com/nbt/journal/v27/n7/extref/nbt.1551-S2.xls

perl ~/Scripts/fig_table/xlsx2csv.pl -f nbt.1551-S2.xls --sheet 'counts' \
    | perl -nla -F/,/ -e '
    if ( /^#/ ) {
        print qq{#ORF\tGene\tAvg};
    }
    elsif ( /^\d/ ) {
        print qq{$F[1]\t$F[2]\t$F[9]};
    }
    ' \
    > mrna_levels.tsv
    
```

## ess, rich/minimal and chem

* ess: Giaever, G., et al. Functional Profiling of theSaccharomyces
  cerevisiae Genome. Nature 418, 387-391. (2002)
* rich/minimal: Mechanisms of Haploinsufficiency Revealed by Genome-Wide
  Profiling in Yeast Deutschbauer, AM. et al. GENETICS April 1, 2005
  vol. 169 no. 4 1915-1925; 10.1534/genetics.104.036871
* chem: The Chemical Genomic Portrait of Yeast: Uncovering a Phenotype
  for All Genes. Hillenmeyer, M.E. et al. Science 18 Apr 2008: Vol. 320,
  Issue 5874, pp. 362-365 DOI: 10.1126/science.1150021

```bash
mkdir -p ~/data/mrna-structure/gene-traits
cd ~/data/mrna-structure/gene-traits

wget -N http://www-sequence.stanford.edu/group/yeast_deletion_project/Essential_ORFs.txt

cat Essential_ORFs.txt \
    | perl -nl -e '
        next unless /^\d/;
        my $orf = ( split /\s+/ )[1];
        print qq{$orf};
    ' \
    > ess_orf.tsv

wget -N http://www-sequence.stanford.edu/group/research/HIP_HOP/supplements/01yfh/files/OrfGeneData.txt

cat OrfGeneData.txt \
    | perl -nla -F"\t" -e '
        printf q{#} if /^orf/;
        print qq{$F[0]\t$F[1]\t$F[5]\t$F[13]\t$F[17]};
    ' \
    > rich_orf.tsv

# http://chemogenomics.stanford.edu/supplements/global/

wget -N http://chemogenomics.stanford.edu/supplements/global/download/data/hom.z_tdist_pval_nm.counts.smallmol.cutoff.01.xls

perl ~/Scripts/fig_table/xlsx2csv.pl -f hom.z_tdist_pval_nm.counts.smallmol.cutoff.01.xls --sheet 'hom.z_tdist_pval_nm.smallmol.co' \
    | perl -nla -F"," -MText::CSV_XS -e '
    BEGIN { 
        our $csv = Text::CSV_XS->new(); 
        print qq{#ORF\tGene\tCount};
    }
    
    if ( /^Y/ ) {
        if ($csv->parse($_)) {
            my @fields = $csv->fields();
            print qq{$fields[0]\t$fields[1]\t$fields[6]};
        }
    }
    ' \
    > chem_orf.tsv

```

## Recombination rates

* Data source: Global mapping of meiotic recombination hotspots and
  coldspots in the yeast Saccharomyces cerevisiae. vol. 97 no. 21 PNAS
  Jennifer L. Gerton, 11383–11390

```bash
mkdir -p ~/data/mrna-structure/gene-traits
cd ~/data/mrna-structure/gene-traits

wget -N http://derisilab.ucsf.edu/data/hotspots/forWebORFs.txt

cat forWebORFs.txt \
    | perl -nla -F"\t" -MStatistics::Lite -e '
        next unless /^Y/;    # ORF stable id start with a "Y"
        next if @F < 2;
        my $rec_rate  = Statistics::Lite::median(grep {defined} @F[1 .. 7]);
        print qq{$F[0]\t$rec_rate};
    ' \
    > rec_rate.tsv
    
```

## Protein-protein interactions

* Data source:

```bash
mkdir -p ~/data/mrna-structure/gene-traits
cd ~/data/mrna-structure/gene-traits

wget -N http://drygin.ccbr.utoronto.ca/%7Ecostanzo2009/sgadata_costanzo2009_stringentCutoff_101120.txt.gz

gzip -d -c sgadata_costanzo2009_stringentCutoff_101120.txt.gz \
    | perl -nla -F"\t" -e '
        BEGIN { our %interact_of; }
        
        next unless /^Y/;    # ORF stable id start with a "Y"
        next if @F != 7;
        $interact_of{$F[0]}++;
        
        END {
            for my $key ( sort keys %interact_of ) {
                print qq{$key\t$interact_of{$key}};
            }
        }
    ' \
    > interact_count.tsv
    
```

# Phylogeny

## create protein coding gene list

```bash
mkdir -p ~/data/mrna-structure/phylogeny
cd ~/data/mrna-structure/phylogeny

# sgd/saccharomyces_cerevisiae.gff → protein coding gene list

perl ~/Scripts/pars/program/protein_coding_list.pl --file ~/data/mrna-structure/sgd/saccharomyces_cerevisiae.gff --output protein_coding_list.csv

perl ~/Scripts/pars/program/protein_coding_list_range.pl --file ~/data/mrna-structure/sgd/saccharomyces_cerevisiae.gff --output protein_coding_list_range.csv

perl ~/Scripts/pars/program/protein_coding_list_range_chr.pl --file ~/data/mrna-structure/sgd/saccharomyces_cerevisiae.gff --output protein_coding_list_range_chr.csv

```

## cut cds alignment

### create cds_yml

```bash
cd ~/data/mrna-structure/phylogeny
mkdir -p ~/data/mrna-structure/phylogeny/gene_cds_yml

perl ~/Scripts/pars/program/cut_cds_yml.pl --file protein_coding_list_range_chr.csv --output gene_cds_yml

```

### cut cds_alignment by cds_yml

```bash

export NAME=Scer_n7_Spar
cp -rf ~/data/mrna-structure/alignment/scer_wgs/${NAME}_refined ~/data/mrna-structure/phylogeny/${NAME}_refined
cd ~/data/mrna-structure/phylogeny/${NAME}_refined
gunzip -rfvc *.maf.gz.fas.gz > species.fas
mkdir -p ~/data/mrna-structure/phylogeny/${NAME}_gene_alignment_cds
cd ~/data/mrna-structure/phylogeny/${NAME}_gene_alignment_cds
cat ../protein_coding_list.csv |
   parallel --line-buffer -j 8 '
   	   fasops slice ../${NAME}_refined/species.fas ../gene_cds_yml/{}.yml -n S288c -o {}.fas.fas
   '
unset NAME

export NAME=Scer_n7p_Spar
cp -rf ~/data/mrna-structure/alignment/scer_wgs/${NAME}_refined ~/data/mrna-structure/phylogeny/${NAME}_refined
cd ~/data/mrna-structure/phylogeny/${NAME}_refined
gunzip -rfvc *.maf.gz.fas.gz > species.fas
mkdir -p ~/data/mrna-structure/phylogeny/${NAME}_gene_alignment_cds
cd ~/data/mrna-structure/phylogeny/${NAME}_gene_alignment_cds
cat ../protein_coding_list.csv |
   parallel --line-buffer -j 8 '
   	   fasops slice ../${NAME}_refined/species.fas ../gene_cds_yml/{}.yml -n S288c -o {}.fas.fas
   '
unset NAME

export NAME=Scer_n157_Spar
cp -rf ~/data/mrna-structure/alignment/scer_wgs/${NAME}_refined ~/data/mrna-structure/phylogeny/${NAME}_refined
cd ~/data/mrna-structure/phylogeny/${NAME}_refined
gunzip -rfvc *.maf.gz.fas.gz > species.fas
mkdir -p ~/data/mrna-structure/phylogeny/${NAME}_gene_alignment_cds
cd ~/data/mrna-structure/phylogeny/${NAME}_gene_alignment_cds
cat ../protein_coding_list.csv |
   parallel --line-buffer -j 8 '
   	   fasops slice ../${NAME}_refined/species.fas ../gene_cds_yml/{}.yml -n S288c -o {}.fas.fas
   '
unset NAME

export NAME=Scer_n157_nonMosaic_Spar
cp -rf ~/data/mrna-structure/alignment/scer_wgs/${NAME}_refined ~/data/mrna-structure/phylogeny/${NAME}_refined
cd ~/data/mrna-structure/phylogeny/${NAME}_refined
gunzip -rfvc *.maf.gz.fas.gz > species.fas
mkdir -p ~/data/mrna-structure/phylogeny/${NAME}_gene_alignment_cds
cd ~/data/mrna-structure/phylogeny/${NAME}_gene_alignment_cds
cat ../protein_coding_list.csv |
   parallel --line-buffer -j 8 '
   	   fasops slice ../${NAME}_refined/species.fas ../gene_cds_yml/{}.yml -n S288c -o {}.fas.fas
   '
unset NAME

export NAME=Scer_n157_nonMosaic_consensus
cp -rf ~/data/mrna-structure/alignment/scer_wgs/${NAME}_refined ~/data/mrna-structure/phylogeny/${NAME}_refined
cd ~/data/mrna-structure/phylogeny/${NAME}_refined
gunzip -rfvc *.maf.gz.fas.gz > species.fas
mkdir -p ~/data/mrna-structure/phylogeny/${NAME}_gene_alignment_cds
cd ~/data/mrna-structure/phylogeny/${NAME}_gene_alignment_cds
cat ../protein_coding_list.csv |
   parallel --line-buffer -j 8 '
   	   fasops slice ../${NAME}_refined/species.fas ../gene_cds_yml/{}.yml -n S288c -o {}.fas.fas
   '
unset NAME

export NAME=Scer_n128_Spar
cp -rf ~/data/mrna-structure/alignment/scer_wgs/multi128_Spar/${NAME}_refined ~/data/mrna-structure/phylogeny/${NAME}_refined
cd ~/data/mrna-structure/phylogeny/${NAME}_refined
gunzip -rfvc *.maf.gz.fas.gz > species.fas
mkdir -p ~/data/mrna-structure/phylogeny/${NAME}_gene_alignment_cds
cd ~/data/mrna-structure/phylogeny/${NAME}_gene_alignment_cds
cat ../protein_coding_list.csv |
   parallel --line-buffer -j 16 '
   	   fasops slice ../${NAME}_refined/species.fas ../gene_cds_yml/{}.yml -n S288c -o {}.fas.fas
   '
unset NAME

```

### count cds_alignment proporation in sgd

```bash

export NAME=Scer_n7_Spar
cd ~/data/mrna-structure/phylogeny
perl ~/Scripts/pars/program/count_gene_range.pl --file protein_coding_list_range.csv --dir ${NAME}_gene_alignment_cds --output ${NAME}_gene_range.csv
unset NAME

export NAME=Scer_n7p_Spar
cd ~/data/mrna-structure/phylogeny
perl ~/Scripts/pars/program/count_gene_range.pl --file protein_coding_list_range.csv --dir ${NAME}_gene_alignment_cds --output ${NAME}_gene_range.csv
unset NAME

export NAME=Scer_n157_Spar
cd ~/data/mrna-structure/phylogeny
perl ~/Scripts/pars/program/count_gene_range.pl --file protein_coding_list_range.csv --dir ${NAME}_gene_alignment_cds --output ${NAME}_gene_range.csv
unset NAME

export NAME=Scer_n157_nonMosaic_Spar
cd ~/data/mrna-structure/phylogeny
perl ~/Scripts/pars/program/count_gene_range.pl --file protein_coding_list_range.csv --dir ${NAME}_gene_alignment_cds --output ${NAME}_gene_range.csv
unset NAME

export NAME=Scer_n157_nonMosaic_consensus
cd ~/data/mrna-structure/phylogeny
perl ~/Scripts/pars/program/count_gene_range.pl --file protein_coding_list_range.csv --dir ${NAME}_gene_alignment_cds --output ${NAME}_gene_range.csv
unset NAME

export NAME=Scer_n128_Spar
cd ~/data/mrna-structure/phylogeny
perl ~/Scripts/pars/program/count_gene_range.pl --file protein_coding_list_range.csv --dir ${NAME}_gene_alignment_cds --output ${NAME}_gene_range.csv
unset NAME
```

## create gene_phylogeny (n157_nonMosaic)

```bash
export NAME=Scer_n157_nonMosaic_Spar

mkdir -p ~/data/mrna-structure/phylogeny/${NAME}_newick
cd ~/data/mrna-structure/phylogeny/${NAME}_newick

cat ../protein_coding_list.csv |
  parallel --line-buffer -j 1 '
      if [ -e {}.nwk ]; then
          echo >&2 '    {}.nwk already presents  '
          exit;
      fi
      egaz raxml ../${NAME}_gene_alignment_cds/{}.fas.fas --seed 999 --tmp . --parallel 8 --outgroup Spar -o {}.nwk
  '
unset NAME
```

## count distance (n157_nonMosaic)

```bash
export NAME=Scer_n157_nonMosaic_Spar

cd ~/data/mrna-structure/phylogeny
echo 'S288c,beer001,beer002,beer003,beer004,beer005,beer006,beer007,beer008,beer009,beer010,beer011,beer012,beer013,beer014,beer015,beer016,beer020,beer021,beer022,beer023,beer024,beer025,beer026,beer027,beer028,beer029,beer030,beer031,beer032,beer033,beer034,beer036,beer037,beer038,beer039,beer040,beer041,beer043,beer044,beer045,beer046,beer047,beer048,beer049,beer050,beer051,beer052,beer053,beer054,beer055,beer056,beer059,beer061,beer062,beer063,beer064,beer065,beer066,beer067,beer068,beer069,beer070,beer071,beer073,beer075,beer076,beer077,beer078,beer079,beer080,beer081,beer082,beer083,beer084,beer085,beer086,beer087,beer088,beer089,beer090,beer091,beer092,beer094,beer095,beer096,beer097,beer098,beer099,beer100,beer101,beer102,bioethanol001,bioethanol003,bioethanol004,bread001,bread002,bread003,bread004,sake001,sake002,sake003,sake004,sake005,sake006,sake007,spirits001,spirits002,spirits003,spirits004,spirits005,spirits011,wine001,wine003,wine004,wine005,wine006,wine007,wine009,wine010,wine011,wine012,wine013,wine014,wine015,wine017,wine018,wild004,wild005,wild006,wild007,Spar' \
| tr "," "\n" \
> ${NAME}_strain_name.list

mkdir ~/data/mrna-structure/phylogeny/${NAME}_distance
cat protein_coding_list.csv |
   parallel --line-buffer -j 8 '
    if [ -e "${NAME}_newick/{}.nwk" ]; then
       perl ~/Scripts/pars/program/count_distance.pl --file ~/data/mrna-structure/phylogeny/${NAME}_newick/{}.nwk --list ~/data/mrna-structure/phylogeny/${NAME}_strain_name.list --output ${NAME}_distance/{}.csv
    fi
    '
cd ~/data/mrna-structure/phylogeny
echo 'gene,asian,wine,beer1,beer2,mixed,beer1_beer2,beer1_beer2_mixed,wine_beer1_beer2_mixed' > ${NAME}_mean_distance.csv
cat protein_coding_list.csv |
    parallel --line-buffer -j 8 '
        if [ -e "${NAME}_distance/{}.csv" ]; then
            perl ~/Scripts/pars/program/count_distance_mean.pl --file ${NAME}_distance/{}.csv
        fi
    ' \
>> ${NAME}_mean_distance.csv

# generate domestication gene （asian → nodom_gene, wine → limdom_gene, beer1+beer2 → strdom_gene)
Rscript ~/Scripts/pars/program/${NAME}_distance_processed.R

unset NAME

```

```bash
#  生成alignment_proporation_1.list

export NAME=Scer_n7_Spar
Rscript ~/Scripts/pars/program/${NAME}_distance_processed.R
unset NAME

export NAME=Scer_n7p_Spar
Rscript ~/Scripts/pars/program/${NAME}_distance_processed.R
unset NAME

export NAME=Scer_n157_Spar
Rscript ~/Scripts/pars/program/${NAME}_distance_processed.R
unset NAME

export NAME=Scer_n128_Spar
Rscript ~/Scripts/pars/program/${NAME}_distance_processed.R
unset NAME

```

# SNP

## count per gene GC content

```bash

export NAME=Scer_n7_Spar
mkdir -p ~/data/mrna-structure/result/$NAME
cd ~/data/mrna-structure/result/$NAME
perl ~/Scripts/pars/program/count_ACGT_percent.pl --file ~/data/mrna-structure/process/$NAME.gene_variation.process.yml --varfold ~/data/mrna-structure/process/$NAME.gene_variation.fold_class.tsv --output $NAME.gene_variation.fold_class.csv
unset NAME

export NAME=Scer_n7p_Spar
mkdir -p ~/data/mrna-structure/result/$NAME
cd ~/data/mrna-structure/result/$NAME
perl ~/Scripts/pars/program/count_ACGT_percent.pl --file ~/data/mrna-structure/process/$NAME.gene_variation.process.yml --varfold ~/data/mrna-structure/process/$NAME.gene_variation.fold_class.tsv --output $NAME.gene_variation.fold_class.csv
unset NAME

export NAME=Scer_n157_Spar
mkdir -p ~/data/mrna-structure/result/$NAME
cd ~/data/mrna-structure/result/$NAME
perl ~/Scripts/pars/program/count_ACGT_percent.pl --file ~/data/mrna-structure/process/$NAME.gene_variation.process.yml --varfold ~/data/mrna-structure/process/$NAME.gene_variation.fold_class.tsv --output $NAME.gene_variation.fold_class.csv
unset NAME

export NAME=Scer_n157_nonMosaic_Spar
mkdir -p ~/data/mrna-structure/result/$NAME
cd ~/data/mrna-structure/result/$NAME
perl ~/Scripts/pars/program/count_ACGT_percent.pl --file ~/data/mrna-structure/process/$NAME.gene_variation.process.yml --varfold ~/data/mrna-structure/process/$NAME.gene_variation.fold_class.tsv --output $NAME.gene_variation.fold_class.csv
unset NAME

export NAME=Scer_n128_Spar
mkdir -p ~/data/mrna-structure/result/$NAME
cd ~/data/mrna-structure/result/$NAME
perl ~/Scripts/pars/program/count_ACGT_percent.pl --file ~/data/mrna-structure/process/$NAME.gene_variation.process.yml --varfold ~/data/mrna-structure/process/$NAME.gene_variation.fold_class.tsv --output $NAME.gene_variation.fold_class.csv
unset NAME

```

## count SNPs and gene 

```bash

Rscript -e 'install.packages("getopt", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN")'
Rscript -e 'install.packages("ape", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN")'
Rscript -e 'install.packages("ggplot2", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN")'
Rscript -e 'install.packages("scales", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN")'
Rscript -e 'install.packages("reshape", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN")'
Rscript -e 'install.packages("pander", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN")'
Rscript -e 'install.packages("gridExtra", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN")'
Rscript -e 'install.packages("plyr", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN")'
Rscript -e 'install.packages("dplyr", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN")'
Rscript -e 'install.packages("proto", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN")'
Rscript -e 'install.packages("gsubfn", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN")'
Rscript -e 'install.packages("RSQLite", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN")'
Rscript -e 'install.packages("sqldf", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN")'

export NAME=Scer_n7_Spar
cd ~/data/mrna-structure/result/$NAME
Rscript ~/Scripts/pars/program/${NAME}_stat_SNPs.R
sed -i "" "s/-&gt;/->/g" data_SNPs_PARS_*.csv  # debug "->"
unset NAME

export NAME=Scer_n7p_Spar
cd ~/data/mrna-structure/result/$NAME
Rscript ~/Scripts/pars/program/${NAME}_stat_SNPs.R
sed -i "" "s/-&gt;/->/g" data_SNPs_PARS_*.csv  # debug "->"
unset NAME

export NAME=Scer_n157_Spar
cd ~/data/mrna-structure/result/$NAME
Rscript ~/Scripts/pars/program/${NAME}_stat_SNPs.R
sed -i "" "s/-&gt;/->/g" data_SNPs_PARS_*.csv  # debug "->"
unset NAME

export NAME=Scer_n157_nonMosaic_Spar
cd ~/data/mrna-structure/result/$NAME
Rscript ~/Scripts/pars/program/${NAME}_stat_SNPs.R
sed -i "" "s/-&gt;/->/g" data_SNPs_PARS_*.csv  # debug "->"
unset NAME

export NAME=Scer_n128_Spar
cd ~/data/mrna-structure/result/$NAME
Rscript ~/Scripts/pars/program/${NAME}_stat_SNPs.R
sed -i "" "s/-&gt;/->/g" data_SNPs_PARS_*.csv  # debug "->"
unset NAME

```

## count A/T <->G/C

```bash

export NAME=Scer_n7_Spar
cd ~/data/mrna-structure/result/$NAME
mkdir -p ~/data/mrna-structure/result/$NAME/freq_each
Rscript ~/Scripts/pars/program/${NAME}_count_AT_GC.R
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_each/PARS_cds_stat.csv --output freq_each/PARS_cds_stat_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_each/PARS_utr_stat.csv --output freq_each/PARS_utr_stat_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_each/PARS_syn_stat.csv --output freq_each/PARS_syn_stat_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_each/PARS_nsy_stat.csv --output freq_each/PARS_nsy_stat_chi_square.csv
unset NAME

export NAME=Scer_n7p_Spar
cd ~/data/mrna-structure/result/$NAME
mkdir -p ~/data/mrna-structure/result/$NAME/freq_each
Rscript ~/Scripts/pars/program/${NAME}_count_AT_GC.R
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_each/PARS_cds_stat.csv --output freq_each/PARS_cds_stat_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_each/PARS_utr_stat.csv --output freq_each/PARS_utr_stat_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_each/PARS_syn_stat.csv --output freq_each/PARS_syn_stat_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_each/PARS_nsy_stat.csv --output freq_each/PARS_nsy_stat_chi_square.csv
unset NAME

export NAME=Scer_n157_Spar
cd ~/data/mrna-structure/result/$NAME
mkdir -p ~/data/mrna-structure/result/$NAME/freq_each
mkdir -p ~/data/mrna-structure/result/$NAME/freq_10
Rscript ~/Scripts/pars/program/${NAME}_count_AT_GC.R
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_10/PARS_cds_stat_freq_10.csv --output freq_10/PARS_cds_stat_freq_10_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_10/PARS_utr_stat_freq_10.csv --output freq_10/PARS_utr_stat_freq_10_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_10/PARS_syn_stat_freq_10.csv --output freq_10/PARS_syn_stat_freq_10_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_10/PARS_nsy_stat_freq_10.csv --output freq_10/PARS_nsy_stat_freq_10_chi_square.csv
unset NAME

export NAME=Scer_n157_nonMosaic_Spar
cd ~/data/mrna-structure/result/$NAME
mkdir -p ~/data/mrna-structure/result/$NAME/freq_each
mkdir -p ~/data/mrna-structure/result/$NAME/freq_10
Rscript ~/Scripts/pars/program/${NAME}_count_AT_GC.R
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_10/PARS_cds_stat_freq_10.csv --output freq_10/PARS_cds_stat_freq_10_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_10/PARS_utr_stat_freq_10.csv --output freq_10/PARS_utr_stat_freq_10_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_10/PARS_syn_stat_freq_10.csv --output freq_10/PARS_syn_stat_freq_10_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_10/PARS_nsy_stat_freq_10.csv --output freq_10/PARS_nsy_stat_freq_10_chi_square.csv
unset NAME


export NAME=Scer_n128_Spar
cd ~/data/mrna-structure/result/$NAME
mkdir -p ~/data/mrna-structure/result/$NAME/freq_each
mkdir -p ~/data/mrna-structure/result/$NAME/freq_10
Rscript ~/Scripts/pars/program/${NAME}_count_AT_GC.R
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_10/PARS_cds_stat_freq_10.csv --output freq_10/PARS_cds_stat_freq_10_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_10/PARS_utr_stat_freq_10.csv --output freq_10/PARS_utr_stat_freq_10_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_10/PARS_syn_stat_freq_10.csv --output freq_10/PARS_syn_stat_freq_10_chi_square.csv
perl ~/Scripts/pars/program/count_stem_loop_chi_square.pl --file freq_10/PARS_nsy_stat_freq_10.csv --output freq_10/PARS_nsy_stat_freq_10_chi_square.csv
unset NAME

```

## count stem length selection
```bash
export NAME=Scer_n157_nonMosaic_Spar
cd ~/data/mrna-structure/result/$NAME 
perl ~/Scripts/pars/program/count_position_gene.pl --file ~/data/mrna-structure/process/$NAME.gene_variation.process.yml --origin data_SNPs_PARS_cds.csv --output data_SNPs_PARS_cds_pos.csv
Rscript ~/Scripts/pars/program/$NAME_count_AT_GC_gene_trait.R
unset NAME

```

## count per gene cds_utr
```bash
export NAME=Scer_n157_nonMosaic_Spar
perl ~/Scripts/pars/program/count_cut_range.pl --file ~/data/mrna-structure/process/$NAME.gene_variation.process.yml --cut ~/data/mrna-structure/process/sce_cds.yml --output stem_loop_cds_length.csv 
perl ~/Scripts/pars/program/count_cut_range.pl --file ~/data/mrna-structure/process/$NAME.gene_variation.process.yml --cut ~/data/mrna-structure/process/sce_utr.yml --output stem_loop_utr_length.csv
perl ~/Scripts/pars/program/count_per_gene_ACGT_percent.pl --file data_SNPs_PARS_cds.csv --output data_SNPs_PARS_cds_per_gene_ATGC.csv
perl ~/Scripts/pars/program/count_per_gene_ACGT_percent.pl --file data_SNPs_PARS_utr.csv --output data_SNPs_PARS_utr_per_gene_ATGC.csv
Rscript ~/Scripts/pars/program/$NAME_cds_utr.R
unset NAME

```
