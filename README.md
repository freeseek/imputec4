imputec4
========

A protocol to impute C4 alleles from MHC genotypes computed from genotype array or whole genome sequence data using the HapMap3 CEU reference panel from:
```
Sekar, A., McCarroll S., et al. Schizophrenia risk from complex
variation of complement component 4. Nature 530, 177–183 (2016)
```
This reference panel was generated using droplet digital PCR as explained <a href="http://mccarrolllab.org/wp-content/uploads/2014/12/Molecular-genetic-analysis-of-C4-structural-variation.pdf">here</a>. For non-European populations and lower-frequency alleles efforts are underway to create a more advanced reference panel from whole genome sequencing data: check <a href="http://mccarrolllab.org/resources/">here</a> for updates. For any feedback, send an email to giulio.genovese@gmail.com or mccarroll@genetics.med.harvard.edu

![](http://mccarrolllab.org/wp-content/uploads/2014/12/C4haplotypes-300x183.png)

Installation
============

Install basic tools (Debian/Ubuntu specific):
```
sudo apt install wget gzip unzip samtools openjdk-11-jre-headless
```

Optionally, you can install these libraries to activate further HTSlib features:
```
sudo apt install libbz2-dev libssl-dev liblzma-dev libgsl0-dev
```

Preparation steps
```
mkdir -p $HOME/bin $HOME/res && cd /tmp
```

Download latest version of <a href="https://github.com/samtools/htslib">HTSlib</a> and <a href="https://github.com/samtools/bcftools">BCFtools</a> (if not downloaded already)
```
git clone --branch=develop git://github.com/samtools/htslib.git
git clone --branch=develop git://github.com/samtools/bcftools.git
```

Compile latest version of HTSlib (optionally disable bz2, gcs, and lzma) and BCFtools (make sure you are using gcc version 5 or newer)
```
cd htslib && autoheader && (autoconf || autoconf) && ./configure --disable-bz2 --disable-gcs --disable-lzma && make && cd ..
cd bcftools && make && cd .
/bin/cp bcftools/{bcftools,plugins/fixref.so} $HOME/bin/
```

Make sure the directory with the plugins is available to bcftools
```
export PATH="$HOME/bin:$PATH"
export BCFTOOLS_PLUGINS="$HOME/bin"
```

Download Beagle binary
```
wget -P $HOME/res/ https://faculty.washington.edu/browning/beagle/beagle.r1399.jar
```

Install latest development version of plink 1.9
```
wget https://www.cog-genomics.org/static/bin/plink/plink_linux_x86_64_dev.zip
unzip -od $HOME/bin/ plink_linux_x86_64_dev.zip plink
```

Download GRCh37 human genome reference
```
wget -O- ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz | \
  gzip -d > $HOME/res/human_g1k_v37.fasta
samtools faidx $HOME/res/human_g1k_v37.fasta
```

Download GRCh38 human genome reference
```
wget -O- ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz | \
  gzip -d > $HOME/res/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
samtools faidx $HOME/res/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
```

Download C4 reference panel in Beagle 3 format and an additional file with marker positions for the GRCh37 human genome reference
```
wget -P $HOME/res/ http://mccarrolllab.com/wp-content/uploads/2014/12/MHC_haplotypes_CEU_HapMap3_ref_panel.bgl
wget -P $HOME/res/ https://raw.githubusercontent.com/freeseek/imputec4/master/MHC_haplotypes_CEU_HapMap3_ref_panel.GRCh37
```

Liftover marker positions for the GRCh38 human genome reference
```
wget -P $HOME/bin/ http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
chmod a+x $HOME/bin/liftOver
wget -P $HOME/res/ http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
awk 'NR>1 {print "chr6\t"$2-1"\t"$2"\t"$1}' $HOME/res/MHC_haplotypes_CEU_HapMap3_ref_panel.GRCh37 | \
  $HOME/bin/liftOver /dev/stdin $HOME/res/hg19ToHg38.over.chain.gz /dev/stdout /dev/stderr | \
  awk 'BEGIN {print "id\tpos"} {print $4"\t"$3}' > $HOME/res/MHC_haplotypes_CEU_HapMap3_ref_panel.GRCh38
```

Generate C4 reference panels in VCF format for both the GRCh37 and GRCh38 human genome references
```
declare -A fasta=( ["37"]="$HOME/res/human_g1k_v37.fasta" ["38"]="$HOME/res/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna" )
declare -A chrom=( ["37"]="6" ["38"]="chr6" )
declare -A length=( ["37"]="171115067" ["38"]="170805979" )

for build in 37 38; do
  (echo "##fileformat=VCFv4.2"; \
  echo "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"; \
  echo "##contig=<ID=${chrom[$build]},length=${length[$build]}>"; \
  tr '\r' '\n' < $HOME/res/MHC_haplotypes_CEU_HapMap3_ref_panel.bgl | \
    grep C4 | cut -f3- | tr '\t' '\n' | sort | uniq | awk '{print "##ALT=<ID="$0">"}'; \
  echo "##reference=${fasta[$build]}"; \
  echo -en "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"; \
  tr '\r' '\n' < $HOME/res/MHC_haplotypes_CEU_HapMap3_ref_panel.bgl | \
    paste $HOME/res/MHC_haplotypes_CEU_HapMap3_ref_panel.GRCh$build - | \
    awk -v chr="${chrom[$build]}" 'NR==1 {for (i=5; i<NF; i+=2) printf "\t"$i; printf "\n"} \
    NR>1 {printf chr"\t"$2"\t"$4; delete x; delete y; j=0; \
    for (i=5; i<=NF; i++) if (!($i in x)) {x[$i]=j; y[j++]=$i} \
    if ($4=="C4") {printf "\tG\t<"y[0]">"; for (i=1; i<j; i++) printf ",<"y[i]">" } \
    else { printf "\t"y[0]"\t"y[1]; for (i=2; i<j; i++) printf ","y[i] } \
    printf "\t.\t.\t.\tGT"; \
    if ($4=="C4") for (i=5; i<NF; i+=2) printf "\t"1+x[$i]"|"1+x[$(i+1)]; \
    else for (i=5; i<NF; i+=2) printf "\t"x[$i]"|"x[$(i+1)]; \
    printf "\n"}') | \
    bcftools +fixref --no-version -Ov \
    -o $HOME/res/MHC_haplotypes_CEU_HapMap3_ref_panel.GRCh$build.vcf -- \
    --fasta-ref ${fasta[$build]} --mode flip

  bcftools norm --no-version --check-ref w --fasta-ref ${fasta[$build]} \
    $HOME/res/MHC_haplotypes_CEU_HapMap3_ref_panel.GRCh$build.vcf -o /dev/null 2>&1 | \
    grep REF_MISMATCH | cut -f2,3 | awk -F"\t" -v OFS="\t" 'NR==FNR {x[$2]++} \
    NR>FNR && $2 in x {alt=$4; ref=$5; $4=ref; $5=alt; \
    for (i=10; i<=NF; i++) $i=1-substr($i,1,1)"|"1-substr($i,3,1)} NR>FNR' \
    - $HOME/res/MHC_haplotypes_CEU_HapMap3_ref_panel.GRCh$build.vcf | \
    bcftools view --no-version -Oz \
    -o $HOME/res/MHC_haplotypes_CEU_HapMap3_ref_panel.GRCh$build.vcf.gz
done
```

Run C4 imputation from an input VCF
===================================

Run imputation using Beagle
```
vcf="..."
out="..."
build=38 # build=37
declare -A chr=( ["37"]="6" ["38"]="chr6" )
declare -A reg=( ["37"]="6:24894177-33890574" ["38"]="chr6:24893949-33922797" )

bcftools view --no-version "$vcf" -r ${reg[$build]} | \
  java -Xmx8g -jar $HOME/res/beagle.r1399.jar gt=/dev/stdin \
  ref=$HOME/res/MHC_haplotypes_CEU_HapMap3_ref_panel.GRCh$build.vcf.gz out="$out" \
  map=<(awk -v chr=${chr[$build]} 'NR>1 {print chr"\t.\t"$2/1e7"\t"$2}' \
  $HOME/res/MHC_haplotypes_CEU_HapMap3_ref_panel.GRCh$build)
```
Notice that, due to the special nature of the MHC region, we select a genetic map with almost no recombination allowed

Extract imputed C4 alleles into a table
```
out="..."
build=38 # build=37
declare -A reg=( ["37"]="6:31948000-31948000" ["38"]="chr6:31980223-31980223" )

bcftools index "$out.vcf.gz" && \
bcftools query -f "[%SAMPLE\t%ALT\t%GT\n]" "$out.vcf.gz" -r ${reg[$build]} | tr -d '[<>]' | \
  awk -F"\t" -v OFS="\t" '{split($2,a,","); a["0"]="NA"; split($3,b,"|"); \
  print $1,a[b[1]],a[b[2]]}' > "$out.tsv"
```
Notice that due to the location of the C4 gene in the MHC locus, the imputed C4 alleles are particularly susceptible to potential confounding due to: (i) linkage disequilibrium mediated correlation with genotypes at other MHC variants, including HLA variants; (ii) population stratification as MHC haplotypes are prone to high allele frequency differences across ethnic groups due to strong natural selection at the MHC locus.

Check for consistency of the C4 reference panel
===============================================

Convert the C4 and 1000 Genomes project reference panels to plink and then merge to compute consistency
```
url="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/integrated_call_sets/ALL.chr6.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz"

bcftools annotate --no-version -x ID -I +'%CHROM:%POS:%REF:%ALT' \
  $HOME/res/MHC_haplotypes_CEU_HapMap3_ref_panel.GRCh37.vcf.gz | \
  $HOME/bin/plink --vcf /dev/stdin --biallelic-only --keep-allele-order --const-fid --make-bed \
  --out MHC_haplotypes_CEU_HapMap3_ref_panel.GRCh37

bcftools view --no-version $url -r 6:24894177-33890574 | \
  awk 'NF==2 {print "##contig=<ID=6,length=171115067>"} {print}' | \
  bcftools view --no-version -v snps | \
  bcftools annotate --no-version -x ID -I +'%CHROM:%POS:%REF:%ALT' | \
  $HOME/bin/plink --vcf /dev/stdin --keep-allele-order --const-fid --make-bed \
  --out ALL.chr6.integrated_phase1_v3.20101123.snps_indels_svs.genotypes

$HOME/bin/plink --bfile MHC_haplotypes_CEU_HapMap3_ref_panel.GRCh37 \
  --bmerge ALL.chr6.integrated_phase1_v3.20101123.snps_indels_svs.genotypes --merge-mode 6
```

You should get the following result
```
577275 overlapping calls, 577275 nonmissing in both filesets.
575472 concordant, for a concordance rate of 0.996877.
```
