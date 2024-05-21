# VDJbase repertoire annotation and downstream analysis


The pipeline performs anotation and downstrean analysis of AIRR-seq.

The pipeline can be devided into seven main componenet:

**1. Initial repertoire alignment and annotation based reference set**

> In this section, the repertoire sequences are annotated using IgBlast and MakeDb (presto) against the supplied reference set and collapse tham.

**2. Undocumented allele inference**

> In this section, inference for undocumented V allele (not found in the initial reference set) is performed using TIgGER.

**3. Second repertoire alignment and annotation with the discovered undocumented alleles.**

> In this section, in case undocumented alleles were inferred the repertoire is re-aligned and annotated with the additional alleles.

**4. Clonal inference and selection of colonal representative**

> In this section, clones are infered for the annotated repertoire, and a single representative for each clone whith the least number of mutation is chosen.

**5. Genotype inference**

> In this section, genotypes are inferred for each of the IG calls: V, D, and J using TIgGER Bayesian inferernce tool, and a personal reference set is created.

**6. Third repertoire alignment and annotation with the personal reference set.**

> In this section, the repertoire sequences are annotated using IgBlast and MakeDb (presto) against the personal reference set from step 5.

**7. OGRDB statistics.**

> the reperotire statistics are deduced using ogrdbstats.


### Input files:

1. An AIRR-seq in fasta format.
2. heavy_chain (yes/no)
3. chain - (IGHV/IGLV/IGKV)
4. Reference set files for IG*V, IG*D, and IG*J alleles in fasta format for the right chain.
5. auxiliary_data and custom_internal_data for the right chain.

### Output:

1. The aligned repertoire with the personal reference set.
2. The aligned repertoire with the init reference set.
3. A genotype reports
4. A log files
5. An ogrdb statistics report for the init aligned repertoire
6. An ogrdb statistics report for the final aligned repertoire
7. pipeline_statistic.csv - table of pass and fail reads for some of the steps.

### Docker images: 

The pipeline uses two docker images:

1. peresay/suite
2. williamlees/ogrdbstats



