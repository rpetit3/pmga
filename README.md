# PLEASE NOTE - THIS FORK IS PROVIDED AS-IS

This fork is provided as is, I'm not the maintainer, just a bioinformatician that wanted to use
[PMGA (PubMLST Genome Annotator v2.0)](https://github.com/CDCgov/BMGAP/tree/master/pipeline/PMGA)
on the command-line. The folks at the [CDC](https://github.com/CDCgov), in particular
[Nadav Topaz](https://github.com/ntopaz) deserves all the kudos and credit for
creating PMGA!

Given I only wanted to run this on the command-line, please keep that in mind if any issues arise.
I might be able to assist with technical issues, but I cannot help interpret the results
(e.g. serotype/serogroup). For help with the biological significance of the results, it would
probably be best to reach out to  [BMGAP](https://github.com/CDCgov/BMGAP).

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
your *Neisseria* and *H. influenzae* samples! For this you will use `pmga`

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

| Extension                       | Description                                           |
|---------------------------------|-------------------------------------------------------|
| `*-allele-matrix.txt`           | A file with the allele ID for each loci with a hit    |
| `*-blast-final-results.json.gz` | Filtered BLAST results in JSON format                 |
| `*-blast-raw-results.json.gz`   | Unfiltered BLAST results in JSON format               |
| `*-loci-counts.txt`             | A file with the number of hits per loci               |
| `*.gff.gz`                      | A GFF3 file annotated with the BALST results          |
| `*.txt`                         | The final predicted serotype/serogroup for the sample |
