# Example analysis of exome data from IGSR/1000 Genomes

This document describes utilization of this workflow on **exome** example data downloaded from IGSR/1000 Genomes using reference genome version **GRCh38**.

## 1. Download data
### 1.1. Download git repository

Download git repository

	git clone 


### 1.2. Download human genome reference GRCh38/HG38

Download human genome reference GRCh38/HG38 bundle.

Read more:
- [Download human genome reference GRCh38/HG38](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle)

### 1.3. Select sample(s)

Browse and select wanted sample(s) from the IGSR [data portal](https://www.internationalgenome.org/data-portal/sample)

In this example we will choose sample **HG00188**.

### 1.4. Download sequencing data for sample HG00188 - Raw sequencing data / FASTQ files

We have selected sample HG00188 from IGSR data portal.

Under data collections choose **1000 Genomes on GRCh38**. On the left under **Data types** choose **Sequence**. Under **Technologies** choose **Exome**. Now under File you have seeminlgy paired-end sequencing 

Using linux terminal let's download the samples:

	mkdir -p samples/HG00188
	cd samples/HG00188

	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR070/SRR070795/SRR070795_2.fastq.gz
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR070/SRR070504/SRR070504_2.fastq.gz
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR070/SRR070795/SRR070795_1.fastq.gz
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR070/SRR070504/SRR070504_1.fastq.gz

Read more:
- [About FASTQ sequence read files](https://www.internationalgenome.org/faq/about-fastq-sequence-read-files/)

### 1.5. Exome / Exon targetted intervals

Since we are going to analyze exome sequencing data which is targeted data to specific intervals in the genome, we need those target/bait intervals for our analysis.

Let's check the IGSR [FAQ](https://www.internationalgenome.org/faq)
- [How was exome and exon targetted sequencing used?](https://www.internationalgenome.org/faq/how-was-exome-and-exon-targetted-sequencing-used/)

Download the provided targets.

	cd ../..
	mkdir intervals
	cd intervals
	wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/exome_pull_down_targets/20130108.exome.targets.bed.README
	wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/exome_pull_down_targets/20130108.exome.targets.bed

Unfortunately 20130108.exome.targets.bed.README says the target intervals are in GRCh37 reference genome version coordinates. We want GRCh38 reference genome version. One option would be to use liftOver-tool to lift over the GRCh37 to GRCh38. Let's continue to dig around.

#### 1.5.1. Search for GRCh38 exome target coordinates

Let's check the articles related to IGSR and the 1000 Genomes Project [here](https://www.internationalgenome.org/about). Read the [IGSR’s Nucleic Acids Research publication](https://academic.oup.com/nar/article/48/D1/D941/5580898)

> 1000 Genomes Project data updated to GRCh38
> Sequence data produced by the 1000 Genomes Project was realigned to the GRCh38 assembly (7). This comprised low coverage WGS and exome data for 2548 unrelated samples and an additional 150 related samples. The alignments have been used in calling variants on GRCh38, in a process that used BCFtools, GATK and FreeBayes for site discovery and produced a joint-genotyped, integrated, phased, biallelic SNV call set in late 2018 (manuscript submitted, https://doi.org/10.12688/wellcomeopenres.15126.1). Since then, an updated call set, adding INDELs, has been released. All data from this work is available on the IGSR FTP site and all variants have been submitted to the European Variation Archive (EVA).

Alright let's check out the manuscript 

> "manuscript submitted, https://doi.org/10.12688/wellcomeopenres.15126.1". 

The site actually prompts the visitor of [newer version](https://wellcomeopenresearch.org/articles/4-50/v2). In **Methods** section under **Quality control of alignment files** there is a link to exome target coordinates, [ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/working/20190125_coords_exon_target/](ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/working/20190125_coords_exon_target/)

Let's download the files from the ftp-link:

	wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/working/20190125_coords_exon_target/20190125_exon_coord_MANIFEST.txt
	wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/working/20190125_coords_exon_target/20190125_exon_coord_README.txt
	wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/working/20190125_coords_exon_target/output_1000G_Exome.v1.bed

File **20190125_exon_coord_README.txt** says that this interval file is a result of actual liftOver of GRCh37 version. Well, alright, let's use this new interval file anyway.

### 1.6. Search for additional information

We need some extra information to do the analysis.

Per sample we need following information for the **input_fofn** file

	file(cols[2]), // fastq_1
    file(cols[3]), // fastq_2
    cols[7], // run_date
    cols[1], // sample_name
    // NEW
    cols[4], // lane
    cols[5], // library_name
    cols[8], // platform_name
    cols[6], // platform_unit
    cols[0], // readgroup_name
    cols[9] // sequencing_center

We have the fastq_1 and fastq_2 files. We are missing **run_date**, **lane**, **library_name**, **platform_name**, **platform_unit**, **readgroup_name**, **sequencing_center** information.

Browse around the ftp-server...
Check [http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/1000genomes.sequence.index](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/1000genomes.sequence.index)...

Download

	cd ..
	wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/1000genomes.sequence.index

It's basically a big table, let's read it into R and select only relevant information for us:

```r
IGSR_seq <- read.csv(
  file = "path_where_you_downloaded_this_file/1000genomes.sequence.index",
  header = TRUE,
  sep = "\t",
  # comment.char = "##"
  skip = 28
)

#View(IGSR_seq)

#install.packages("dplyr")
library(dplyr)

# Let's filter by our sample and partial filenames
specific_samples <- IGSR_seq %>% 
  filter(SAMPLE_NAME == "HG00188") %>%
  filter(grepl('SRR070/SRR070504/SRR070504',`X.FASTQ_ENA_PATH`) | 
           grepl('SRR070/SRR070795/SRR070795',`X.FASTQ_ENA_PATH`)
         )

# discard few unneeded/uninformative columns
discarded_cols <- c("MD5","STUDY_NAME","POPULATION", "RUNBLOCK_NAME", "INSERT_SIZE", "LIBRARY_LAYOUT", "WITHDRAWN", "WITHDRAWN_DATE", "COMMENT", "ANALYSIS_GROUP")
subset_df <- specific_samples[ , -which(names(specific_samples) %in% discarded_cols)]

# output to markdown type table
knitr::kable(subset_df, "pipe")

```

Output is:

|X.FASTQ_ENA_PATH                                                         |RUN_ID    |STUDY_ID  |CENTER_NAME |SUBMISSION_ID |SUBMISSION_DATE     |SAMPLE_ID |SAMPLE_NAME |EXPERIMENT_ID |INSTRUMENT_PLATFORM |INSTRUMENT_MODEL            |LIBRARY_NAME |RUN_NAME |RUN_BLOCK_NAME |PAIRED_FASTQ                                                             |READ_COUNT |BASE_COUNT |
|:------------------------------------------------------------------------|:---------|:---------|:-----------|:-------------|:-------------------|:---------|:-----------|:-------------|:-------------------|:---------------------------|:------------|:--------|:--------------|:------------------------------------------------------------------------|:----------|:----------|
|ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR070/SRR070504/SRR070504_1.fastq.gz |SRR070504 |SRP004058 |WUGSC       |SRA025233     |2010-10-28 00:00:00 |SRS006906 |HG00188     |SRX029520     |ILLUMINA            |Illumina Genome Analyzer II |2862178100   |65845    |NA             |ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR070/SRR070504/SRR070504_2.fastq.gz |27855659   |5571131800 |
|ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR070/SRR070504/SRR070504_2.fastq.gz |SRR070504 |SRP004058 |WUGSC       |SRA025233     |2010-10-28 00:00:00 |SRS006906 |HG00188     |SRX029520     |ILLUMINA            |Illumina Genome Analyzer II |2862178100   |65845    |NA             |ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR070/SRR070504/SRR070504_1.fastq.gz |27855659   |5571131800 |
|ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR070/SRR070795/SRR070795_1.fastq.gz |SRR070795 |SRP004058 |WUGSC       |SRA025589     |2010-10-28 00:00:00 |SRS006906 |HG00188     |SRX029733     |ILLUMINA            |Illumina Genome Analyzer II |2862178100   |65853    |NA             |ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR070/SRR070795/SRR070795_2.fastq.gz |27986688   |5597337600 |
|ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR070/SRR070795/SRR070795_2.fastq.gz |SRR070795 |SRP004058 |WUGSC       |SRA025589     |2010-10-28 00:00:00 |SRS006906 |HG00188     |SRX029733     |ILLUMINA            |Illumina Genome Analyzer II |2862178100   |65853    |NA             |ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR070/SRR070795/SRR070795_1.fastq.gz |27986688   |5597337600 |

Inspect the table by filenames SRR070504/SRR070504 ... and SRR070795/SRR070795. 

Read more here:
- [https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups](https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups)

## 2. Prepare main sample information containing input file for analysis

Next we need to input your sample data to a text file where values are separted by tab. Information needs to placed in the following order.

| readgroup_name | sample_name | fastq_1 | fastq_2 | lane | library_name | platform_unit | run_date | platform_name | sequencing_center |

Read more about different fields here:
- [https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups](https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups)

From our previous section with we select appropriate information from 1000genomes.sequence.index per our sample and more specifically sample run.

**readgroup_name** can basically be anything that is unique to sequencing run, for example RUN_ID, SUBMISSION_ID, EXPERIMENT_ID or RUN_NAME for example. For this we'll choose RUN_NAME: 65853 and 65845 (as per specific filename).

**sample_name** is HG00188.

**fastq_1** and **fastq_2** are as per RUN_NAME. Give full path to your downloaded files.

**lane** since we don't exactly have real lanes we'll use the RUN_ID as lane designation. SRR070504 and SRR070795.

**library_name** let's use the LIBRARY_NAME values. 2862178100 for both, so the different sequencing runs are from the same library prepration.

**platform_unit** is kind of tricky. We can check from our FASTQ-files if they have any usabale information:

	zcat <path_to_your_file>/samples/HG00188/ftp.sra.ebi.ac.uk/vol1/fastq/SRR070/SRR070504/SRR070504_1.fastq.gz | head

	@SRR070504.1 HWUSI-EAS729_104869904:8:1:1576:1011/1
	NTAAGGAGCTCAGGGTTGTTTTCTGAAGCGAAAATGCAGGCAGATGAGCATAGGCTGAGCCAGGTTCCCAGAAGAGTAACAGTGGGAGCTGGTCTCCAGC
	+
	!",,)21/01@@@@@@@@@@22222@@@@@@@@@@@@@@@@@@@@@@@@@@@:@@#############################################
	@SRR070504.2 HWUSI-EAS729_104869904:8:1:1666:1008/1
	NGAAAATAAACTTGTCCAAGGTATAAAACACATAATTCTGGAGTGTTACAATGACAGTAGGGAGTGTTCTGGATAGTTCCAACTTTAATAATAAAAGCCA
	+
	!,,-*770002@C@@C@@@C@@@@@@@2222@@C@@@@CCC@C@C@@@@@@CCC@@@@@@<<<<<@@@@@@@@@@@CC@@@@@@@C@8C@@@@@222@8:
	@SRR070504.3 HWUSI-EAS729_104869904:8:1:1747:1007/1
	NCGCCTGCCCTGAGGACCGGTGGGACGCCCTTTCTCTCCTTGCGCCCAGACTGCAGAAGAGGAAGGACACGCAGAGACGGAGGCGCAGGAAGAATAGAGG

Here you can read about FASTQ-format:
- [https://help.basespace.illumina.com/articles/descriptive/fastq-files/](https://help.basespace.illumina.com/articles/descriptive/fastq-files/)
- [https://en.wikipedia.org/wiki/FASTQ_format](https://en.wikipedia.org/wiki/FASTQ_format)

First part @SRR070504 is SRA accession ID, number after "." is likely a read number.  
HWUSI-EAS729_104869904 is likely combination of unique instrument name (HWUSI-EAS729) with flowcell barcode/name (104869904).  
Next after : , we have flowcell lane, 8.
Next tile number within the flowcell.
Next 'x'-coordinate of the cluster within the tile.
Next 'y'-coordinate of the cluster within the tile.
And /1 is the member of a pair.

According to GATK read groups article, platform_unit value consists of:

>{FLOWCELL_BARCODE}.{LANE}.{SAMPLE_BARCODE}.

So in this case flowcell_barcode would be 104869904. Lane would be 8. sample_barcode is unavailable maybe leave it out. So what we get ofr platform_unit value is HWUSI-EAS729_104869904:8 and HWUSI-EAS729_104869904:7. The HWUSI-EAS729_ part is probably unnecessary.

**run_date** is 2010-10-28. Adding the specifc time 00:00:00, unfortunately ended up with an error in one of the workflow processes/tools.

**platform_name** is ILLUMINA.

**sequencing_center** is WUGSC.

Let's compile our table

readgroup_name | sample_name | fastq_1 | fastq_2 | lane | library_name | platform_unit | run_date | platform_name | sequencing_center
--- | --- | --- | --- | --- | --- | --- | --- | --- | ---
65845 | HG00188 | path_to_my_files/HG00188/ftp.sra.ebi.ac.uk/vol1/fastq/SRR070/SRR070504/SRR070504_1.fastq.gz | path_to_my_files/HG00188/ftp.sra.ebi.ac.uk/vol1/fastq/SRR070/SRR070504/SRR070504_2.fastq.gz | SRR070504 | 2862178100 | HWUSI-EAS729_104869904:8 | 2010-10-28 | ILLUMINA | WUGSC
65853 | HG00188 | ppath_to_my_files/HG00188/ftp.sra.ebi.ac.uk/vol1/fastq/SRR070/SRR070795/SRR070795_1.fastq.gz | path_to_my_files/HG00188/ftp.sra.ebi.ac.uk/vol1/fastq/SRR070/SRR070795/SRR070795_2.fastq.gz | SRR070795 | 2862178100 | HWUSI-EAS729_104869904:7 | 2010-10-28 | ILLUMINA | WUGSC

Save file as tab separated text file. (Note! above table is formatted to markdown (github)).

Let's automate the previous steps a little bit in R as far as we can.

```r
# select specific columns in specific order
selected_df <- subset_df[,c("RUN_NAME", "SAMPLE_NAME", "X.FASTQ_ENA_PATH", "PAIRED_FASTQ", "RUN_ID", "LIBRARY_NAME", "SUBMISSION_DATE", "INSTRUMENT_PLATFORM", "CENTER_NAME")]
selected_df

# add platform_unit column in right order
library(tibble)
selected_df <- add_column(selected_df, plaform_unit = NA, .after = "LIBRARY_NAME")

# change the name of the columns
column_names <- c("readgroup_name", "sample_name", "fastq_1", "fastq_2", "lane", "library_name", "platform_unit", "run_date", "platform_name", "sequencing_center")
names(selected_df) <-  column_names

# save the file
write.table(selected_df, file = "<path_to_my_files>/manifest-1-sample.txt", quote = FALSE, sep = "\t", col.names = TRUE)
```

Now you need to fill in the platform_unit column values and modify the actual FASTQ-file columns fastq_1, fastq_2.

### 2.1. Lanes or just a separate sequencing run with same library?

As in the previous section we checked information for our sample per its' files/sequencing runs. It seems that the same sample library (identifcal library designation) was fully run multiple times. So the sample was simply sequenced in two separate occasions. When sample is split across lanes it might be to in case of one lane fails.

Lanes or not, I reckon that combining sequence data from same library source would be similar. Some might just combine the fastq-files together, but in this workflow we'll follow Broad Institute's advice.

Read more here:
- [https://gatk.broadinstitute.org/hc/en-us/articles/360035889471-How-should-I-pre-process-data-from-multiplexed-sequencing-and-multi-library-designs-](https://gatk.broadinstitute.org/hc/en-us/articles/360035889471-How-should-I-pre-process-data-from-multiplexed-sequencing-and-multi-library-designs-)

### 2.2. Convert bed to interval_list, scatter the interval_list to multiple files

We need to convert our .bed formatted interval list file to interval_list format.

	cd intervals
	picard BedToIntervalList \
	I=<path_to_my_files>/intervals/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/working/20190125_coords_exon_target/output_1000G_Exome.v1.bed \
	O=<path_to_my_files>/output_1000G_Exome.v1.interval_list \
	SD=<path_to_my_files>/hg38/Homo_sapiens_assembly38.dict
	
Further more we need to scatter the interval_list to multiple files for our workflow.

	gatk SplitIntervals \
	-R <path_to_my_files>/hg38/v0/v0/Homo_sapiens_assembly38.fasta \
	-L <path_to_my_files>/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/working/20190125_coords_exon_target/output_1000G_Exome.v1.interval_list \
	--scatter-count 30 \
	--subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION \
	-O <path_to_my_files>/intervals/output_1000G_Exome_scattered_intervals

Unfortunately the scattered file names are renamed as .intervals which causes the gatk to error out. We need to rename the files to .interval_list. We also need to generate single file which has the paths of all the scattered .interval_list files.

Create a small renaming script:

	cat #!/bin/bash \
	for file in $(find . -name "*$1"); do \
		mv "$file" "${file%$1}$2" \
	done > rename.sh

Now enter change the extensions with the script:

	./rename.sh .intervals .interval_list

Generate .txt file with full paths of splitted interval list

	cd output_1000G_Exome_scattered_intervals/
	readlink -f 00* > ../scattered_20190125_coords_exon_targets.txt

### 2.3. Fill in your inputs to params.yaml

**unmapped_bams_list:** - you can comment this out with #

**input_fofn:** - table with sample information we assembled earlier, updated manifest-1-sample.txt.  
**outdir** - output directory where we want to save the end result(s).  
**scattered_calling_interval** - path to file that contains scattered interval paths. Same as in the previous section (scattered_20190125_coords_exon_targets.txt).  

Reference files:
**fasta**  
**dbSNP_vcf:**  
**known_indels_mills:**  
**known_indels_dbSNP:**  
- fill in matching reference files you have downloaded.  
**sequence_grouping:**  
**sequence_grouping_unmapped:**  
- you can find sequence_grouping files from the repository (gatk4-germline-snps-indels-1/resources/original_public_sources/sequence_grouping ...")  


Additional input:  
**temp_dir:** - path to temporary folder, helps when executing locally or on a shared server.

## 3. Run the workflow
### 3.1. Run analysis on local computer

Go to folder where you want the intermediate files to be stored or alternatively specify -workdir <your_chosen_path> argument for nextflow command. Singularity is used in this example instead of docker.

	nextflow run <path_to_my_files>/gatk4-germline-snps-indels-1/main.nf \
	-params-file <path_to_my_files>/gatk4-germline-snps-indels-1/params.local.yaml 
	-c <path_to_my_files>/gatk4-germline-snps-indels-1/nextflow_singularity-local.config -profile singularity,local

You can check the config for run options, such as limiting cpus or memory values for specific executors or profiles as in the for example local.


### 3.2. Run analysis on slurm cluster

To run analysis on slurm cluster, login to slurm login-node. Singularity is used instead of docker. Then we just execute nextflow with slurm specific configuration options.

	nextflow run <path_to_my_files>/gatk4-germline-snps-indels-1/main.nf \
	-params-file <path_to_my_files>/gatk4-germline-snps-indels-1/params.local.yaml 
	-c <path_to_my_files>/gatk4-germline-snps-indels-1/nextflow_singularity-local.config -profile singularity,slurm
