# PLEASE NOTE - THIS FORK IS PROVIDED AS-IS

This fork is provided as is, I'm not the maintainer, just a bioinformatician that wanted to use
[PMGA (PubMLST Genome Annotator v2.0)](https://github.com/CDCgov/BMGAP/tree/master/pipeline/PMGA)
on the command-line. The folks at the [CDC](https://github.com/CDCgov), in particular
[Nadav Topaz](https://github.com/ntopaz) deserves all the kudos and credit for
creating PMGA!

Given I only wanted to run this on the command-line, please keep that in mind if any issues arise.
I might be able to assist with technical issues, but I cannot help interpret the results
(e.g. serotype/serogroup). For help with the biological significance of the results, it would
probably be best to reach out to [Bacterial Meningitis Genome Analysis Platform (BMGAP)](https://github.com/CDCgov/BMGAP).

I'm sure they would be willing to help!

To review the original docs I recommend you visit: [PMGA (PubMLST Genome Annotator v2.0)](https://github.com/CDCgov/BMGAP/tree/master/pipeline/PMGA)

---

## PMGA - PubMLST Genome Annotator

PMGA is a tool for serotyping/serogrouping all *Neisseria* species and *Haemophilus influenzae*. PMGA requires you to
build a set of BLAST databases using data that is available from [PubMLST](https://pubmlst.org/). You can then query
your uncompressed FASTA files against the BLAST databases. If you do not provide a species using the `--species` option,
Mash will be used to determine the species. Once complete, you will get output files (discussed below).

## Installation

Currently `pmga` is available on my personal Conda channel, but is expected to become available on BioConda soon.

```{bash}
mamba create -n pmga -c pmga -c conda-forge -c bioconda pmga
```

## Building BLAST Databases

Prior to running `pmga`, you will need to build the BLAST databases with `pmga-build`. During this process, `pmga-build`
will download each loci for all schemes associated with *Neisseria* scpecies and *Haemophilus influenzae*. Building the
databases can take upwards to 90 minutes, but this is a one time step. After `pmga-build` is completed, by defualt 
all the outputs will be written to `./pubmlst_dbs_all`.

### `pmga-build` Usage

```{bash}
pmga-build --help
usage: pmga-build [-h] [--outdir STR] [--force] [--silent] [--debug] [--version]

pmga-build - Script for creating local BlastDBs from PubMLST alleles

options:
  -h, --help    show this help message and exit
  --outdir STR  Directory to save BLAST databases to (Default: ./pubmlst_dbs_all)
  --force       Overwrite existing directories.
  --silent      Only critical errors will be printed.
  --debug       Print debug related text.
  --version     show program's version number and exit
```

## Running `pmga`

After you have sucessfully build all the necessary BLAST databases, you are now ready to start serotyping/serogrouping
your *Neisseria* and *H. influenzae* samples! `pmga` can be executed on an **uncompressed** FASTA file and will output
a number of files.

### `pmga` Usage

```{bash}
pmga --help
usage: pmga [-h] [--prefix STR] [--blastdir STR] [--species STR] [-t INT] [-o STR] [--force] [--verbose]
            [--silent] [--version]
            FASTA

pmga - Serotyping, serotyping and MLST of all Neisseria species and Haemophilus influenzae

positional arguments:
  FASTA                 Input FASTA file to analyze

options:
  -h, --help            show this help message and exit
  --prefix STR          Prefix for outputs (Default: Use basename of input FASTA file)
  --blastdir STR        Directory containing BLAST DBs built by pmga-build (Default: ./pubmlst_dbs_all

Additional Options:
  --species STR         Use this as the input species (Default: use Mash distance). Available Choices:
                        neisseria, hinfluenzae
  -t INT, --threads INT
                        Number of cores to use (default=1)
  -o STR, --outdir STR  Directory to output results to (Default: ./pmga)
  --force               Force overwrite existing output file
  --verbose             Print debug related text.
  --silent              Only critical errors will be printed.
  --version             show program's version number and exit
```

### `pmga` Output Files

Below are the expected outputs from each `pmga` run. The files will be output to the value set by `--outdir`, which
defaults to `./pmga`. Each of the output files will use the value of `--prefix` for filenames. By default the prefix
is set to the basename of the input file.

```{bash}
<OUTDIR>
├── <PREFIX>-allele-matrix.txt
├── <PREFIX>-final-blast-results.json.gz
├── <PREFIX>-loci-counts.txt
├── <PREFIX>-raw-blast-results.json.gz
├── <PREFIX>.gff.gz
└── <PREFIX>.txt
```

| Extension                       | Description                                                                  |
|---------------------------------|------------------------------------------------------------------------------|
| `*-allele-matrix.txt`           | A tab-delimitted file with the allele ID for each loci with a hit            |
| `*-blast-final-results.json.gz` | Filtered BLAST results in JSON format                                        |
| `*-blast-raw-results.json.gz`   | Unfiltered BLAST results in JSON format                                      |
| `*-loci-counts.txt`             | A tab-delimitted file with the number of hits per loci                       |
| `*.gff.gz`                      | A GFF3 file annotated with the BALST results                                 |
| `*.txt`                         | A tab-delimitted file with final predicted serotype/serogroup for the sample |

#### Example Serotype/Serogroup Output

```{bash}
sample	species	prediction	genes_present	notes
GCF_003355215	neisseria_serogroup	B	csb,cssA,cssB,cssC,ctrA,ctrB,ctrC,ctrD,ctrE,ctrF,tex	B backbone: All essential capsule genes intact and present
```

Above an example of the output predictions, which contains 5 columns. The 5 columns are

| Column        | Description                                                                      |
|---------------|----------------------------------------------------------------------------------|
| sample        | The input sample name                                                            |
| species       | The species of the sample either `neisseria_serogroup` or `hinfluenzae_serotype` |
| prediction    | The predicted serotype or serogroup                                              |
| genes_present | A list of genes present in the sample                                            |
| notes         | Any notes associated with the prediction                                         |

## Citations

If you make use of this tool, please cite the following:

* **[Bacterial Meningitis Genome Analysis Platform (BMGAP)](https://github.com/CDCgov/BMGAP)**  
An analysis pipeline, ExpressJS API, and ReactJS webapp for the analysis and characterization of bacterial meningitis samples  
*Buono SA, Kelly RJ, Topaz N, Retchless AC, Silva H, Chen A, Ramos E, Doho G, Khan AN, Okomo-Adhiambo MA, Hu F, Marasini D, Wang X. [Web-Based Genome Analysis of Bacterial Meningitis Pathogens for Public Health Applications Using the Bacterial Meningitis Genomic Analysis Platform (BMGAP).](https://doi.org/10.3389/fgene.2020.601870) Front Genet. 2020 Nov 26;11:601870.*  

* **[BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)**  
Basic Local Alignment Search Tool  
*Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL [BLAST+: architecture and applications](http://dx.doi.org/10.1186/1471-2105-10-421). BMC Bioinformatics 10, 421 (2009)*  

* **[Mash](https://github.com/marbl/Mash)**  
Fast genome and metagenome distance estimation using MinHash  
*Ondov BD, Treangen TJ, Melsted P, Mallonee AB, Bergman NH, Koren S, Phillippy AM [Mash: fast genome and metagenome distance estimation using MinHash](http://dx.doi.org/10.1186/s13059-016-0997-x). Genome Biol 17, 132 (2016)*  

* **[Pigz](https://zlib.net/pigz/)**  
A parallel implementation of gzip for modern multi-processor, multi-core machines.  
*Adler, M. [pigz: A parallel implementation of gzip for modern multi-processor, multi-core machines.](https://zlib.net/pigz/) Jet Propulsion Laboratory (2015).*  

* **[PubMLST.org](https://pubmlst.org/)**  
A database housing MLST shemes for many bacterial species.  
*Jolley KA, Bray JE, Maiden MCJ [Open-access bacterial population genomics: BIGSdb software, the PubMLST.org website and their applications.](http://dx.doi.org/10.12688/wellcomeopenres.14826.1) Wellcome Open Res 3, 124 (2018)*  
