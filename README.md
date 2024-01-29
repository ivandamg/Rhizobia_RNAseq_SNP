# Rhizobia_RNAseq_SNP

Pipeline for SNP identification from RNAseq reads

0. Download and prepare files
1. 
      module load SRA-Toolkit
      prefetch SRR27534840
      fasterq-dump --split-files SRR27534843
      # compress files
      for FILE in $(ls *.fastq); do echo $FILE; sbatch --partition=pshort_el8 --job-name=$(echo $FILE | cut -d'_' -f1,2)GZIP --time=0-03:30:00 --mem-per-cpu=12G --ntasks=1 --cpus-per-task=1 --output=index.out --error=index.error --mail-type=END,FAIL --wrap "cd /data/projects/p495_SinorhizobiumMeliloti/08_OtherRNAseqs/01_Hg_PRJNA1063170/02_Ecotype_Dongbei/01_RawData ; gzip $FILE;"; done
      
      mv SRR27534822_2.fastq.gz Gansu_Rhizobia_root_Hgexpose_Rep2_R-8-2_2.fastq.gz
      mv SRR27534823_2.fastq.gz Gansu_Rhizobia_root_Hgexpose_Rep1_R-8-1_2.fastq.gz
      mv SRR27534824_2.fastq.gz Gansu_Rhizobia_root_Control_Rep3_R-7-3_2.fastq.gz
      mv SRR27534825_2.fastq.gz Gansu_Rhizobia_root_Control_Rep2_R-7-2_2.fastq.gz
      mv SRR27534826_2.fastq.gz Gansu_Rhizobia_root_Control_Rep1_R-7-1_2.fastq.gz
    
    
2. Trimming and quality filter with Fastp


                     
     for FILE in $(ls *1.fastq.gz); do echo $FILE; sbatch --partition=pshort_el8 --job-name=$(echo $FILE | cut -d'_' -f6)fastp --time=0-02:00:00 --mem-per-cpu=64G --ntasks=4 --cpus-per-task=1 --output=$(echo $FILE | cut -d'_' -f6)_fastp.out --error=$(echo $FILE | cut -d'_' -f6)_fastp.error --mail-type=END,FAIL --wrap " cd /data/projects/p495_SinorhizobiumMeliloti/08_OtherRNAseqs/01_Hg_PRJNA1063170/01_Ecotype_Gansu/01_RawData; module load FastQC; ~/00_Software/fastp --in1 $FILE --in2 $(echo $FILE | cut -d'_' -f1,2,3,4,5,6)_2.fastq.gz --out1 ../02_TrimmedData/$(echo $FILE | cut -d'_' -f1,2,3,4,5,6)_1_trimmed.fastq.gz --out2 ../02_TrimmedData/$(echo $FILE | cut -d'_' -f1,2,3,4,5,6)_2_trimmed.fastq.gz -h ../02_TrimmedData/$(echo $FILE | cut -d'.' -f1)_fastp.html --thread 4; fastqc -t 4 ../02_TrimmedData/$(echo $FILE | cut -d'_' -f1,2,3,4,5,6)_1_trimmed.fastq.gz; fastqc -t 4 ../02_TrimmedData/$(echo $FILE | cut -d'_' -f1,2,3,4,5,6)_2_trimmed.fastq.gz"; sleep  1; done

  
  
3.  Map reads with bbmap


        Sinteractive -c 16 -m 128G -t 02:00:00
        module load gcc bbmap/38.63
        bbmap.sh in1=C3-4_Rirregularis_1_val_1.fq in2=C3-4_Rirregularis_2_val_2.fq out=C3-4_mappedSuperZ17p.sam ref=../00_GenomeAssemblies/SuperAssembly_Z17ALL_primary_ctg.fasta
        
        # low quality repetitive assembly
        
        for i in $(ls Unmapped_Mesculenta_*.1.fq); do echo $i; bbmap.sh in1=$i in2=$(echo $i | cut -d'.' -f1).2.fq out=$(echo $i | cut -d'.' -f1)_mapped_MergedDAOM197198B1.sam ref=../../00_GenomeAssemblies/Merged_DAOM197198_B1.fna vslow k=8 maxindel=200 minratio=0.1 ; done


4. transform to bam sort and Qfilter
            
            module load gcc samtools
            samtools sort -m4G -@4 -o Z17Hifireads_RefSuperZ17p_sorted.bam Z17Hifireads_RefSuperZ17p.sam
            samtools view -bq 30 Z17Hifireads_RefSuperZ17p_sorted.bam > Z17Hifireads_RefSuperZ17p_sorted_Q30.bam
            samtools index Z17Hifireads_RefSuperZ17p_sorted.bam


5. FeatureCounts on chimer.

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
