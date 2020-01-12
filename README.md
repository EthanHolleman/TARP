# TARP: Transposable Element Assembly Remapping Pipeline

## PAG Poster Link
[Poster Link](https://drive.google.com/open?id=1H8wGdqWAdffG0Yw6Ii8iBKpZ3ta6S2VX)

## PAG Abstract

Transposable elements are prolific throughout plant genomes and while originally regarded as nothing more than genetic parasites increasing number of studies have shown these elements to contribute significantly to the phenotypes of many significant agricultural species. Therefore, accurate and up to date information on the content of transposable elements in plant, genomes is an essential piece of information for many researchers.  


However, as plant genome assemblies continue to be improved and updates to old versions released the public libraries of transposable elements often remain the same. The repetitive nature of transposable elements means they can easily be left in unassembled scaffolds in early genome versions and then over time make their way into the chromosomal sequences as sequencing data becomes more robust. There is, therefore, a need for computational methods to reliably and accurately locate and remap transposable elements between outdated and current assemblies.  


Here, we present Transposable Element Assembly Remap Pipeline or TARP. TARP utilizes iterative global alignments of transposable element consensus sequences to locate known and novel elements in updated plant genome assemblies without searching for either one of these categories explicitly. This allows researchers interested in the transposable element content of updated plant assemblies to utilize previous transposable element libraries and data in their search instead of starting from scratch.


Additionally, using flanking sequence comparisons TART can distinguish between novel and known elements in a new assembly even if there have been significant locational changes to known elements due to the incorporation of scaffolds into updated assemblies. This allows researchers to compare where elements of interest are located between two or more assembly versions.

## Quick Start Guide

### Gather Assembly Files
Navigate to NCBI and find the assembly entry for the organism of interest. Then
via the FTP site (using FTP is important to insure correct chromosome to accession
number conversion) download the assembled chromsome files and chr2acc file for
the assembly version your elements are mapped to and the updated assembly. Store
these in seperate directories. I recommend using [filezilla](https://filezilla-project.org/) for this.  
To make things easier you can concatinate the individual chromosome files into
a single fasta using the command below.
```
cat *fna > [Complete Genome File Name]
```

### Create Bowtie2 Index
Assuming you have already installed Bowtie2 create a new index using the command
```
bowtie2-build [path to most recent assembly file] [Name of your index]
```
This will build a bowtie2 index in your working directory.

### Build BLAST Databases
Navigate to the directories the current and outdated assemblies are now stored.
In each, run the command below to build the searchable BLAST databases.
```
makeblastdb -in [assembly file name] -dbtype nucl -title [your title] -out
[your db name] -parse_seqids
```

### Gather Backmapping Files
If you wish to use the backmapping functionality you will need a few additional
files. To run standard backmapping you will need a summary file that contains
data on your outdated element library. This should be a tab delimited file; an
example of the format is shown below.
(This requirement will soon be phased out)
```
Element ID	Reference	Class	Subclass	Order	Super Family	Family	Description	Unanchored Scaffold	Start Position	End Position
RLG_Gmr3_Gm1-1	Du et al. 2010 BMC Genomics 2010, 11:113	I	I	LTR	Gypsy	Gmr3	INTACT	Gm01		5622921	5628239
```
Pass this file into TARP with the argument
`-sum [Your file path]`  

Using just the summary file you will be able to run either flanking sequence
based remapping or distance minimization remapping. To run a feature based remap
supply either a file containing genes or vcf formatted SNP file.
Pass this file into TARP with the argument `-F [Your file path]`.
You must also supply an argument for `-sum` to run feature based backmapping.  

Once you decide what method of backmapping to utilize, if any, identify it using
the argument `-M` and refer to the chart below to select method.

| Argument | Method                |
|-----------|-----------------------|
| 1         | Distance Minimization |
| 2         | Flank Based           |
| 3         | Feature Based         |

### Example Command

A generalized TARP command is shown below. Fill in the [] with paths that make
sense for your machine.
```
./Tarp -I [Intact Elements] -S [Solo Elements] -P [Outdated Assembly BLAST DB] 
-C [Current Assembly BLAST DB] -B [Bowtie2 Index] -acc_c [Current acc2chr file]
-acc_o [Outdated acc2chr file] -M [Backmap integer] -sum [Summary File] 
-F [Feature file (if M = 3)]
```
