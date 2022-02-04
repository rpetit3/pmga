#!/usr/bin/env python3
"""
usage: pmga-build [-h] [--outdir STR] [--force] [--silent] [--debug] [--version]

pmga-build - Script for creating local BlastDBs from PubMLST alleles

optional arguments:
  -h, --help    show this help message and exit
  --outdir STR  Directory to save BLAST databases to
  --force       Overwrite existing directories.
  --silent      Only critical errors will be printed.
  --debug       Print debug related text.
  --version     show program's version number and exit
"""
import logging
import os
import sys

VERSION = "3.0.1"
PROGRAM = "pmga-build"
DESCRIPTION = "Script for creating local BlastDBs from PubMLST alleles"
STDOUT = 11
STDERR = 12
logging.addLevelName(STDOUT, "STDOUT")
logging.addLevelName(STDERR, "STDERR")

# Not sure why these are skipped
# should look into if more can be skipped considering >5000 loci are downloaded
# Currently takes ~90 minutes to complete the building
SKIP_ALLELES = [
    "HAEM1147", "HAEM1148", "HAEM1149", "HAEM1150", "HAEM1151", "HAEM1152", "HAEM1153", "HAEM1154", "HAEM1155",
    "HAEM1156", "HAEM1157", "HAEM1158", "HAEM1159", "HAEM1160", "HAEM1161", "HAEM1162", "HAEM1163", "HAEM1164"
]
SHARED_ACROSS_SPECIES = ["NEIS2210"]


def set_log_level(error, debug):
    """Set the output log level."""
    return logging.ERROR if error else logging.DEBUG if debug else logging.INFO


def get_log_level():
    """Return logging level name."""
    return logging.getLevelName(logging.getLogger().getEffectiveLevel())


def execute(cmd, directory=os.getcwd(), capture=False, stdout_file=None, stderr_file=None):
    """A simple wrapper around executor."""
    from executor import ExternalCommand, ExternalCommandFailed
    try:
        command = ExternalCommand(
            cmd, directory=directory, capture=True, capture_stderr=True,
            stdout_file=stdout_file, stderr_file=stderr_file
        )
        command.start()
        if get_log_level() == 'DEBUG':
            logging.log(STDOUT, command.decoded_stdout)
            logging.log(STDERR, command.decoded_stderr)

        if capture:
            return [command.decoded_stdout, command.decoded_stderr]
        return True
    except ExternalCommandFailed as error:
        raise error


def get_neisseria_schemes():
    logging.info("Obtaining pubMLST schemes")
    ### Obtain list of schemes and then each scheme individually from pubMLST for Neisseria sp. and load into dict ###
    schemes = {}
    scheme_url = "http://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/schemes"
    scheme_request_data = requests.get(scheme_url).json()
    for scheme in scheme_request_data["schemes"]:
        url = scheme["scheme"]
        request_data = requests.get(url).json()
        description = request_data["description"]
        if "loci" in request_data:
            schemes[description] = []
            for loci in request_data["loci"]:
                allele = loci.split("/")[-1]
                schemes[description].append(allele)
    logging.info("Schemes obtained")
    return schemes


if __name__ == "__main__":
    import argparse as ap
    import json
    import requests
    from Bio import SeqIO

    parser = ap.ArgumentParser(
        prog=PROGRAM,
        conflict_handler='resolve',
        description=(f'{PROGRAM} (v{VERSION}) - {DESCRIPTION}'),
        formatter_class=ap.RawDescriptionHelpFormatter
    )
    parser.add_argument('--outdir', metavar="STR", type=str, default="pubmlst_dbs_all",
                        help='Directory to save BLAST databases to')
    parser.add_argument('--force', action='store_true', help='Overwrite existing directories.')
    parser.add_argument('--silent', action='store_true',  help='Only critical errors will be printed.')
    parser.add_argument('--debug', action='store_true', help='Print debug related text.')
    parser.add_argument('--version', action='version', version=f'{PROGRAM} {VERSION}')

    args = parser.parse_args()

    # Setup logs
    FORMAT = '%(asctime)s:%(name)s:%(levelname)s - %(message)s'
    logging.basicConfig(format=FORMAT, datefmt='%Y-%m-%d %H:%M:%S',)
    logging.getLogger().setLevel(set_log_level(args.silent, args.debug))

    # Check if outdir exists
    outdir = os.path.abspath(os.path.expanduser(args.outdir))
    if os.path.isdir(args.outdir):
        if args.force:
            logging.info(f"Found --force, removing existing {outdir}")
            execute(f"rm -rf {outdir}")
        else:
            logging.error(f"Output Directory {outdir} aleady exists, please use --force to overwrite")
            sys.exit(1)

    # Setup initial directories
    custom_allele_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "custom_allele_sets")
    logging.info(f"Created Output Directory: {outdir}")
    execute(f"mkdir -p {outdir}")
    execute(f"cp -r {custom_allele_path}/hinfluenzae {outdir}")
    execute(f"cp -r {custom_allele_path}/neisseria {outdir}")

    # Start downloading and building blast databases
    for db in ["neisseria", "hinfluenzae"]:
        request_json = requests.get(f"http://rest.pubmlst.org/db/pubmlst_{db}_isolates/loci?return_all=1").json()

        # Parse through each of the loci
        logging.info(f"Found {len(request_json['loci'])} loci for {db}, database builds might take a while")
        for i, locus in enumerate(request_json["loci"]):
            logging.debug(f"Building database for locus: {locus}")
            locus_url = locus.replace("isolates", "seqdef")
            allele_name = locus_url.split("/")[-1]
            if allele_name in SKIP_ALLELES:
                continue

            # Verify it looks like a FASTA
            fasta_data = requests.get(f"{locus_url}/alleles_fasta").text
            if "{" in fasta_data:
                # no results found
                continue 

            allele_name = allele_name.replace("'", "#")
            allele_db_path = f"{outdir}/{db}/{allele_name}"
            allele_fasta = f"{allele_db_path}/{allele_name}.fasta"
            execute(f"mkdir -p {allele_db_path}")

            # Write the FASTA
            with open(allele_fasta, "w") as fh:
                fh.write(fasta_data)

            # Determine if Nucleotide or Protein
            mol_type = "nucl"
            for letter in next(SeqIO.parse(allele_fasta, "fasta")).seq:
                if letter not in ["A", "G", "C", "T", "N"]:
                    mol_type = "prot"
                    break

            execute(f"makeblastdb -in {allele_fasta} -input_type fasta -title '{allele_name}' -dbtype {mol_type} -out {allele_db_path}/{allele_name}")
            if allele_name in SHARED_ACROSS_SPECIES:
                allele_db_path_hi = allele_db_path.replace("neisseria","hinfluenzae")
                execute(f"cp -r {allele_db_path} {allele_db_path_hi}")

            if i % 1000 == 0:
                logging.info(f"Built database for {i} of {len(request_json['loci'])} loci")

    neisseria_schemes = get_neisseria_schemes()
    with open(f"{outdir}/neisseria_schemes.json", "wt") as json_fh:
        json.dump(neisseria_schemes, json_fh)
