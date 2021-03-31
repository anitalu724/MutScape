export VCF2MAF_URL=`curl -sL https://api.github.com/repos/mskcc/vcf2maf/releases | grep -m1 tarball_url | cut -d\" -f4`
curl -L -o mskcc-vcf2maf.tar.gz $VCF2MAF_URL; tar -zxf mskcc-vcf2maf.tar.gz; cd mskcc-vcf2maf-*
perl vcf2maf.pl --man
perl maf2maf.pl --man

curl -sL https://repo.anaconda.com/miniconda/Miniconda3-py37_4.8.3-Linux-x86_64.sh -o /tmp/miniconda.sh
sh /tmp/miniconda.sh -bfp $HOME/miniconda3

export PATH=$HOME/miniconda3/bin:$PATH

conda install -qy -c conda-forge -c bioconda -c defaults ensembl-vep==100.2 samtools==1.9 bcftools==1.9 ucsc-liftover==377

vep_install --AUTO cf --SPECIES homo_sapiens --ASSEMBLY GRCh38 --CACHEDIR $HOME/.vep

curl -sLO https://raw.githubusercontent.com/Ensembl/ensembl-vep/release/100/examples/homo_sapiens_GRCh38.vcf
vep --species homo_sapiens --assembly GRCh38 --offline --no_progress --no_stats --sift b --ccds --uniprot --hgvs --symbol --numbers --domains --gene_phenotype --canonical --protein --biotype --uniprot --tsl --pubmed --variant_class --shift_hgvs 1 --check_existing --total_length --allele_number --no_escape --xref_refseq --failed 1 --vcf --minimal --flag_pick_allele --pick_order canonical,tsl,biotype,rank,ccds,length --dir $HOME/.vep --fasta $HOME/.vep/homo_sapiens/100_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz --input_file homo_sapiens_GRCh38.vcf --output_file homo_sapiens_GRCh38.vep.vcf --polyphen b --af --af_1kg --af_esp --regulatory
```