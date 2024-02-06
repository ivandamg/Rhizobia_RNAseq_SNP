# Rhizobia_RNAseq_SNP

Pipeline for SNP identification from RNAseq reads


## 0. Download and prepare files


      module load SRA-Toolkit
      prefetch SRR27534840
      fasterq-dump --split-files SRR27534843

      sbatch --partition=pshort_el8 --job-name=SRA --time=0-05:30:00 --mem-per-cpu=12G --ntasks=1 --cpus-per-task=1 --output=index.out --error=index.error --mail-type=END,FAIL --wrap "cd /data/projects/p495_SinorhizobiumMeliloti/08_OtherRNAseqs/01_Fcasuarinae/02_PublishedData/01_RawData; module load SRA-Toolkit; prefetch SRX2582554;  prefetch SRX2582555;  prefetch SRX2582556;  prefetch SRX2582557;  prefetch SRX2582558;  prefetch SRX2582559;  "
 
sbatch --partition=pshort_el8 --job-name=SRA2 --time=0-05:30:00 --mem-per-cpu=12G --ntasks=1 --cpus-per-task=1 --output=index.out --error=index.error --mail-type=END,FAIL --wrap "cd /data/projects/p495_SinorhizobiumMeliloti/08_OtherRNAseqs/01_Fcasuarinae/02_PublishedData/01_RawData; module load SRA-Toolkit; fasterq-dump --split-files SRR5278714;fasterq-dump --split-files SRR5278715;fasterq-dump --split-files SRR5278716;fasterq-dump --split-files SRR5278717;fasterq-dump --split-files SRR5278718;fasterq-dump --split-files SRR5278719;fasterq-dump --split-files SRR5278720;fasterq-dump --split-files SRR5278721;fasterq-dump --split-files SRR5278722;fasterq-dump --split-files SRR5278723;fasterq-dump --split-files SRR5278724;fasterq-dump --split-files SRR5278725"
 



## 1. compress files


      for FILE in $(ls *.fastq); do echo $FILE ;sbatch --partition=pshort_el8 --job-name=GZIP --time=0-03:30:00 --mem-per-cpu=12G --ntasks=1 --cpus-per-task=1 --output=index.out --error=index.error --mail-type=END,FAIL --wrap "cd /data/projects/p495_SinorhizobiumMeliloti/08_OtherRNAseqs/01_Fcasuarinae/02_PublishedData/01_RawData ; gzip $FILE;"; done

## 2. Change names



      mv SRR27534822_2.fastq.gz Gansu_Rhizobia_root_Hgexpose_Rep2_R-8-2_2.fastq.gz
      mv SRR27534823_2.fastq.gz Gansu_Rhizobia_root_Hgexpose_Rep1_R-8-1_2.fastq.gz
      mv SRR27534824_2.fastq.gz Gansu_Rhizobia_root_Control_Rep3_R-7-3_2.fastq.gz
      mv SRR27534825_2.fastq.gz Gansu_Rhizobia_root_Control_Rep2_R-7-2_2.fastq.gz
      mv SRR27534826_2.fastq.gz Gansu_Rhizobia_root_Control_Rep1_R-7-1_2.fastq.gz
    
    
## 3. Trimming and quality filter with Fastp


                     
           for FILE in $(ls *1.fastq.gz); do echo $FILE; sbatch --partition=pshort_el8 --job-name=$(echo $FILE | cut -d'_' -f1,2,3)fastp --time=0-01:00:00 --mem-per-cpu=12G --ntasks=1 --cpus-per-task=1 --output=$(echo $FILE | cut -d'_' -f1,2,3)_fastp.out --error=$(echo $FILE | cut -d'_' -f1,2,3)_fastp.error --mail-type=END,FAIL --wrap " cd /data/projects/p495_SinorhizobiumMeliloti/08_OtherRNAseqs/01_Fcasuarinae/02_PublishedData/01_RawData; module load FastQC; ~/00_Software/fastp --in1 $FILE --in2 $(echo $FILE | cut -d'_' -f1,2,3,4)_2.fastq.gz --out1 ../02_TrimmedData/$(echo $FILE | cut -d'_' -f1,2,3)_1_trimmed.fastq.gz --out2 ../02_TrimmedData/$(echo $FILE | cut -d'_' -f1,2,3)_2_trimmed.fastq.gz -h ../02_TrimmedData/$(echo $FILE | cut -d',' -f1,2,3)_fastp.html --thread 4; fastqc -t 4 ../02_TrimmedData/$(echo $FILE | cut -d'_' -f1,2,3)_1_trimmed.fastq.gz; fastqc -t 4 ../02_TrimmedData/$(echo $FILE | cut -d'_' -f1,2,3)_2_trimmed.fastq.gz"; sleep  1; done



# 4. Mapping to Fcasuarinae with STAR https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

### a. Index reference 

                  sbatch --partition=pshort_el8 --job-name=StarIndex --time=0-01:00:00 --mem-per-cpu=12G --ntasks=1 --cpus-per-task=1 --output=StarIndex.out --error=StarIndex.error --mail-type=END,FAIL --wrap "cd /data/projects/p495_SinorhizobiumMeliloti/08_OtherRNAseqs/01_Fcasuarinae/01_RawData/ncbi_dataset/data/GCF_000013345/; module load STAR/2.7.10a_alpha_220601-GCC-10.3.0; STAR --runThreadN 1 --runMode genomeGenerate --genomeDir /data/projects/p495_SinorhizobiumMeliloti/08_OtherRNAseqs/01_Fcasuarinae/01_RawData/ncbi_dataset/data/GCF_000013345 --genomeFastaFiles GCF_000013345.1_ASM1334v1_genomic.fna --sjdbGTFfile genomic.gtf --sjdbOverhang 99 --genomeSAindexNbases 10"

### b. Map reads 

                  for FILE in $(ls *_1_trimmed.fastq.gz ); do echo $FILE; sbatch --partition=pshort_el8 --job-name=$(echo $FILE | cut -d'_' -f1,2,3)_STAR --time=0-06:00:00 --mem-per-cpu=128G --ntasks=8 --cpus-per-task=1 --output=$(echo $FILE | cut -d'_' -f1,2)_STAR.out --error=$(echo $FILE | cut -d'_' -f1,2)_STAR.error --mail-type=END,FAIL --wrap "module load STAR/2.7.10a_alpha_220601-GCC-10.3.0; cd /data/projects/p495_SinorhizobiumMeliloti/08_OtherRNAseqs/01_Fcasuarinae/02_PublishedData/02_TrimmedData; STAR --runThreadN 8 --genomeDir /data/projects/p495_SinorhizobiumMeliloti/08_OtherRNAseqs/01_Fcasuarinae/01_RawData/ncbi_dataset/data/GCF_000013345 --readFilesIn $FILE $(echo $FILE | cut -d'_' -f1,2,3)_2_trimmed.fastq.gz --readFilesCommand zcat --outFileNamePrefix ../$(echo $FILE | cut -d'_' -f1,2,3)_Mapped --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 1065539232"; sleep  1; done

# 5. Call snps

                for FILE in $(ls *_MappedAligned.sortedByCoord.out.bam); do echo $FILE; sbatch --partition=pshort_el8 --job-name=$(echo $FILE | cut -d'_' -f1,2,3)mpileup --time=0-03:00:00 --mem-per-cpu=64G --ntasks=8 --cpus-per-task=1 --output=$(echo $FILE | cut -d'_' -f1,2,3)_mpileup.out --error=$(echo $FILE | cut -d'_' -f1,2,3)_mpileup.error --mail-type=END,FAIL --wrap "cd /data/projects/p495_SinorhizobiumMeliloti/08_OtherRNAseqs/01_Fcasuarinae/02_PublishedData/03_Mapped; module load SAMtools ; module load BCFtools; bcftools mpileup --threads 8 -a AD,DP,SP -f /data/projects/p495_SinorhizobiumMeliloti/08_OtherRNAseqs/01_Fcasuarinae/01_RawData/ncbi_dataset/data/GCF_000013345/GCF_000013345.1_ASM1334v1_genomic.fna $FILE | bcftools call --threads 8 -mv -Ov -o $(echo $FILE | cut -d'_' -f1,2,3).vcf; bcftools view --threads 8 --exclude 'QUAL <= 30 ' $(echo $FILE | cut -d'_' -f1,2,3).vcf -Oz -o $(echo $FILE | cut -d'_' -f1,2,3)_bcftoolsV1_Q30.vcf.gz"   ; sleep 1; done










######

# END




Sorting and indexing

            for FILE in $(ls *_L*_bwa-mem2_Medicago.bam); do echo $FILE; sbatch --partition=pall --job-name=$(echo $FILE | cut -d'_' -f1,2)ST --time=0-03:00:00 --mem-per-cpu=64G --ntasks=2 --cpus-per-task=1 --output=$(echo $FILE | cut -d'_' -f1,2)_ST.out --error=$(echo $FILE | cut -d'_' -f1,2)_ST.error --mail-type=END,FAIL --wrap "module load UHTS/Analysis/samtools/1.10; cd /data/projects/p495_SinorhizobiumMeliloti/02_DuplexSeq/03_Mapped_Medicago/; samtools sort $FILE -o $(echo $FILE | cut -d'.' -f1)_Sorted.bam; samtools index $(echo $FILE | cut -d'.' -f1)_Sorted.bam; "; sleep 1; done


  
6.  Map reads with bbmap


        Sinteractive -c 16 -m 128G -t 02:00:00
        module load gcc bbmap/38.63
        bbmap.sh in1=C3-4_Rirregularis_1_val_1.fq in2=C3-4_Rirregularis_2_val_2.fq out=C3-4_mappedSuperZ17p.sam ref=../00_GenomeAssemblies/SuperAssembly_Z17ALL_primary_ctg.fasta
        
        # low quality repetitive assembly
        
        for i in $(ls Unmapped_Mesculenta_*.1.fq); do echo $i; bbmap.sh in1=$i in2=$(echo $i | cut -d'.' -f1).2.fq out=$(echo $i | cut -d'.' -f1)_mapped_MergedDAOM197198B1.sam ref=../../00_GenomeAssemblies/Merged_DAOM197198_B1.fna vslow k=8 maxindel=200 minratio=0.1 ; done


7. transform to bam sort and Qfilter
            
            module load gcc samtools
            samtools sort -m4G -@4 -o Z17Hifireads_RefSuperZ17p_sorted.bam Z17Hifireads_RefSuperZ17p.sam
            samtools view -bq 30 Z17Hifireads_RefSuperZ17p_sorted.bam > Z17Hifireads_RefSuperZ17p_sorted_Q30.bam
            samtools index Z17Hifireads_RefSuperZ17p_sorted.bam


8. FeatureCounts on chimer.

1. merge GTF file of both isolates.

            cat DAOM197198_Rhiir2_v2.gff > Merged_DAOM197198_B1.gff
            cat Rhizophagus_irregularis_B1.gff3 | grep -v "#" >> Merged_DAOM197198_B1.gff
            
2. run feature counts on both mapped/reference B1 and DAOM independently

            /work/FAC/FBM/DEE/isanders/popgen_to_var/IM/ZZ_Soft/subread-2.0.3-source/bin/featureCounts -T 8 -a /work/FAC/FBM/DEE/isanders/popgen_to_var/IM/01_Coinoc_v2/03_MappedDAOM/DAOM197198_Rhiir2_v2.gff -t CDS -g ID -o Counts_REFDAOM_COL2215.txt  Unmapped_Mesculenta_COL2215_B1DAOM197198_3_mapped_DAOM197198_Sorted_Q30.bam Unmapped_Mesculenta_COL2215_B1_1_mapped_DAOM197198_Sorted_Q30.bam Unmapped_Mesculenta_COL2215_B1_2_mapped_DAOM197198_Sorted_Q30.bam Unmapped_Mesculenta_COL2215_B1_3_mapped_DAOM197198_Sorted_Q30.bam Unmapped_Mesculenta_COL2215_DAOM197198_1_mapped_DAOM197198_Sorted_Q30.bam Unmapped_Mesculenta_COL2215_DAOM197198_2_mapped_DAOM197198_Sorted_Q30.bam Unmapped_Mesculenta_COL2215_DAOM197198_3_mapped_DAOM197198_Sorted_Q30.bam Unmapped_Mesculenta_COL2215_Mock_1_mapped_DAOM197198_Sorted_Q30.bam Unmapped_Mesculenta_COL2215_Mock_2_mapped_DAOM197198_Sorted_Q30.bam -p --countReadPairs
          


6. Call SNPs from RNA-seq BAM files

             # Call SNP
            module load Bioinformatics/Software/vital-it; module add UHTS/Analysis/freebayes/1.2.0; mkdir 03_SNP_Rhiir;
            for i in $(ls *.bam); do echo $i; freebayes --fasta-reference DAOM.fa --report-monomorphic -C 10 -p 1 $i > 03_SNP_Rhiir_$(echo $i | cut -d'_' -f1,2,3,4,5).vcf; done

            # Filter quality SNP Call

6. Orthologs between B1 and DAOM annotation.

Run orthologs analysis with orthofinder. then in excel identify wich are single-copy orthologs. 
Do normal DE analysis and then look if they are orthologs. always compare amf1 amf2 vs Co-inoculation


7. Then look DE genes with one reference. look at ortholog in other reference and check if they are also DE. If yes, it means that gene is activated in both isolates. If not it means that it is activated in a single-isolate.

8. Can confirm this by looking at SNP in the co-inoculation. We should expect presence of both variants if expressed in both isolates.
9. Coverage analysis 

            module load gcc bedtools2
            for i in $(ls *Q30.bam); do echo $i; bedtools genomecov -ibam $i -bga -split > $(echo $i | cut -d'_' -f3,4,5)_Covdetail.txt ; done


## TEST 2. Mapp reads individually to each genome assembly. and then with variant calling. identify to which isolate the counts are from. 

1. Map cassava unmapped reads to DAOM197198. 
2. Map cassava unmapped reads to B1.
3. Differential expression independently.
4. Variant calling to determine from wich isolate the counts come from: is the gene activated in one isolate or both? This can be done by comparing SNP.
5. chimera not a good aproach because a lot of genes will map together.

6. Review feature counts command to make it better.
