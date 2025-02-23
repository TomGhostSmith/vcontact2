"""vConTACT 2 - Copyright 2018 Benjamin Bolduc, Guilhem Doulcier.

vConTACT2 (viral Contig Automatic Cluster Taxonomy) is tool to perform
"Guilt-by-contig-association" automatic classification of viral
contigs.

This program is distributed under the term of the GNU General Public
Licence v3 (or later) with ABSOLUTELY NO WARRANTY. This is free
software, and you are welcome to redistribute it.
"""

import os
import logging
import argparse
import subprocess
import shutil
import psutil
import pandas as pd
import numpy as np
import pkg_resources  # type: ignore
import sys
sys.path.append(".")
import vcontact2
import vcontact2.protein_clusters
import vcontact2.pcprofiles
import vcontact2.contig_clusters
import vcontact2.evaluations
import vcontact2.cluster_refinements
import vcontact2.modules
import vcontact2.exports.csv
import vcontact2.exports.krona
import vcontact2.exports.cytoscape
import vcontact2.exports.profiles
import vcontact2.exports.summaries

import fcntl
import time

GPULock = "/tmp/GPULock"
CPULock = "/tmp/CPULock"

def acquire_lock():
    """Acquire a lock using a file-based locking mechanism."""
    lock1 = "/tmp/CPULock1"
    lock2 = "/tmp/CPULock2"
    while True:
        try:
            # Attempt to open the file in write mode
            lock_file = open(lock1, 'w')
            # Try to acquire a non-blocking lock (using fcntl)
            fcntl.flock(lock_file, fcntl.LOCK_EX | fcntl.LOCK_NB)
            print(f"got the CPU lock1")
            # return lock_file
            return lock_file
        except BlockingIOError:
            # If the lock file is already locked, wait for a while before retrying
            time.sleep(1)
        try:
            # Attempt to open the file in write mode
            lock_file = open(lock2, 'w')
            # Try to acquire a non-blocking lock (using fcntl)
            fcntl.flock(lock_file, fcntl.LOCK_EX | fcntl.LOCK_NB)
            print(f"got the CPU lock2")
            # return lock_file
            return lock_file
        except BlockingIOError:
            time.sleep(1)


def release_lock(lock_file):
    """Release the lock on the given lock file."""
    print(f"release the CPU lock")
    fcntl.flock(lock_file, fcntl.LOCK_UN)
    lock_file.close()


db_root = pkg_resources.resource_filename("vcontact2", "data/")

ref_dbs = {
    "ProkaryoticViralRefSeq85-ICTV": os.path.join(
        db_root, "ViralRefSeq-prokaryotes-v85.faa.gz"
    ),
    "ArchaeaViralRefSeq85-Merged": os.path.join(
        db_root, "ViralRefSeq-archaea-v85.faa.gz"
    ),
    "ProkaryoticViralRefSeq85-Merged": os.path.join(
        db_root, "ViralRefSeq-prokaryotes-v85.faa.gz"
    ),
    "ProkaryoticViralRefSeq88-Merged": os.path.join(
        db_root, "ViralRefSeq-prokaryotes-v88.faa.gz"
    ),
    "ArchaeaViralRefSeq94-Merged": os.path.join(
        db_root, "ViralRefSeq-archaea-v94.faa.gz"
    ),
    "ProkaryoticViralRefSeq94-Merged": os.path.join(
        db_root, "ViralRefSeq-prokaryotes-v94.faa.gz"
    ),
    "ArchaeaViralRefSeq97-Merged": os.path.join(
        db_root, "ViralRefSeq-archaea-v97.faa.gz"
    ),
    "ProkaryoticViralRefSeq97-Merged": os.path.join(
        db_root, "ViralRefSeq-prokaryotes-v97.faa.gz"
    ),
    "ArchaeaViralRefSeq99-Merged": os.path.join(
        db_root, "ViralRefSeq-archaea-v99.faa.gz"
    ),
    "ProkaryoticViralRefSeq99-Merged": os.path.join(
        db_root, "ViralRefSeq-prokaryotes-v99.faa.gz"
    ),
    "ArchaeaViralRefSeq201-Merged": os.path.join(
        db_root, "ViralRefSeq-archaea-v201.faa.gz"
    ),
    "ProkaryoticViralRefSeq201-Merged": os.path.join(
        db_root, "ViralRefSeq-prokaryotes-v201.faa.gz"
    ),
    "ArchaeaViralRefSeq207-Merged": os.path.join(
        db_root, "ViralRefSeq-archaea-v207.faa.gz"
    ),
    "ArchaeaViralRefSeq211-Merged": os.path.join(
        db_root, "ViralRefSeq-archaea-v211.faa.gz"
    ),
    "ProkaryoticViralRefSeq207-Merged": os.path.join(
        db_root, "ViralRefSeq-prokaryotes-v207.faa.gz"
    ),
    "ProkaryoticViralRefSeq211-Merged": os.path.join(
        db_root, "ViralRefSeq-prokaryotes-v211.faa.gz"
    ),
}
ref_tax = {
    "ProkaryoticViralRefSeq85-ICTV": os.path.join(
        db_root, "ViralRefSeq-prokaryotes-v85.ICTV-reference.csv"
    ),
    "ArchaeaViralRefSeq85-Merged": os.path.join(
        db_root, "ViralRefSeq-archaea-v85.Merged-reference.csv"
    ),
    "ProkaryoticViralRefSeq85-Merged": os.path.join(
        db_root, "ViralRefSeq-prokaryotes-v85.Merged-reference.csv"
    ),
    "ProkaryoticViralRefSeq88-Merged": os.path.join(
        db_root, "ViralRefSeq-prokaryotes-v88.Merged-reference.csv"
    ),
    "ArchaeaViralRefSeq94-Merged": os.path.join(
        db_root, "ViralRefSeq-archaea-v94.Merged-reference.csv"
    ),
    "ProkaryoticViralRefSeq94-Merged": os.path.join(
        db_root, "ViralRefSeq-prokaryotes-v94.Merged-reference.csv"
    ),
    "ArchaeaViralRefSeq97-Merged": os.path.join(
        db_root, "ViralRefSeq-archaea-v97.Merged-reference.csv"
    ),
    "ProkaryoticViralRefSeq97-Merged": os.path.join(
        db_root, "ViralRefSeq-prokaryotes-v97.Merged-reference.csv"
    ),
    "ArchaeaViralRefSeq99-Merged": os.path.join(
        db_root, "ViralRefSeq-archaea-v99.Merged-reference.csv"
    ),
    "ProkaryoticViralRefSeq99-Merged": os.path.join(
        db_root, "ViralRefSeq-prokaryotes-v99.Merged-reference.csv"
    ),
    "ArchaeaViralRefSeq201-Merged": os.path.join(
        db_root, "ViralRefSeq-archaea-v201.Merged-reference.csv"
    ),
    "ProkaryoticViralRefSeq201-Merged": os.path.join(
        db_root, "ViralRefSeq-prokaryotes-v201.Merged-reference.csv"
    ),
    "ArchaeaViralRefSeq207-Merged": os.path.join(
        db_root, "ViralRefSeq-archaea-v207.Merged-reference.csv"
    ),
    "ProkaryoticViralRefSeq207-Merged": os.path.join(
        db_root, "ViralRefSeq-prokaryotes-v207.Merged-reference.csv"
    ),
    "ArchaeaViralRefSeq211-Merged": os.path.join(
        db_root, "ViralRefSeq-archaea-v211.Merged-reference.csv"
    ),
    "ProkaryoticViralRefSeq211-Merged": os.path.join(
        db_root, "ViralRefSeq-prokaryotes-v211.Merged-reference.csv"
    ),
}
ref_g2c = {
    "ProkaryoticViralRefSeq85-ICTV": os.path.join(
        db_root, "ViralRefSeq-prokaryotes-v85.protein2contig.csv"
    ),
    "ArchaeaViralRefSeq85-Merged": os.path.join(
        db_root, "ViralRefSeq-archaea-v85.protein2contig.csv"
    ),
    "ProkaryoticViralRefSeq85-Merged": os.path.join(
        db_root, "ViralRefSeq-prokaryotes-v85.protein2contig.csv"
    ),
    "ProkaryoticViralRefSeq88-Merged": os.path.join(
        db_root, "ViralRefSeq-prokaryotes-v88.protein2contig.csv"
    ),
    "ArchaeaViralRefSeq94-Merged": os.path.join(
        db_root, "ViralRefSeq-archaea-v94.protein2contig.csv"
    ),
    "ProkaryoticViralRefSeq94-Merged": os.path.join(
        db_root, "ViralRefSeq-prokaryotes-v94.protein2contig.csv"
    ),
    "ArchaeaViralRefSeq97-Merged": os.path.join(
        db_root, "ViralRefSeq-archaea-v97.protein2contig.csv"
    ),
    "ProkaryoticViralRefSeq97-Merged": os.path.join(
        db_root, "ViralRefSeq-prokaryotes-v97.protein2contig.csv"
    ),
    "ArchaeaViralRefSeq99-Merged": os.path.join(
        db_root, "ViralRefSeq-archaea-v99.protein2contig.csv"
    ),
    "ProkaryoticViralRefSeq99-Merged": os.path.join(
        db_root, "ViralRefSeq-prokaryotes-v99.protein2contig.csv"
    ),
    "ArchaeaViralRefSeq201-Merged": os.path.join(
        db_root, "ViralRefSeq-archaea-v201.protein2contig.csv"
    ),
    "ProkaryoticViralRefSeq201-Merged": os.path.join(
        db_root, "ViralRefSeq-prokaryotes-v201.protein2contig.csv"
    ),
    "ArchaeaViralRefSeq207-Merged": os.path.join(
        db_root, "ViralRefSeq-archaea-v207.protein2contig.csv"
    ),
    "ArchaeaViralRefSeq211-Merged": os.path.join(
        db_root, "ViralRefSeq-archaea-v211.protein2contig.csv"
    ),
    "ProkaryoticViralRefSeq207-Merged": os.path.join(
        db_root, "ViralRefSeq-prokaryotes-v207.protein2contig.csv"
    ),
    "ProkaryoticViralRefSeq211-Merged": os.path.join(
        db_root, "ViralRefSeq-prokaryotes-v211.protein2contig.csv"
    ),
}


logging.addLevelName(
    logging.WARNING, "\033[1;31m%s\033[1;0m" % logging.getLevelName(logging.WARNING)
)
logging.addLevelName(
    logging.ERROR, "\033[1;41m%s\033[1;0m" % logging.getLevelName(logging.ERROR)
)
logging.addLevelName(
    logging.INFO, "\033[1;42m%s\033[1;0m" % logging.getLevelName(logging.INFO)
)
logging.addLevelName(
    logging.DEBUG, "\033[1;43m%s\033[1;0m" % logging.getLevelName(logging.DEBUG)
)
logging.getLogger().addFilter(logging.Filter("vcontact2"))
logger = logging.getLogger("vcontact2")


# Argparse multiple inheritance
class CustomFormatter(
    argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter
):
    pass


parser = argparse.ArgumentParser(description=__doc__, formatter_class=CustomFormatter)

inputs = parser.add_argument_group("Main Arguments")

inputs.add_argument(
    "-r",
    "--raw-proteins",
    type=str,
    default="test_data/VIRSorter_genomes.faa",
    dest="raw_proteins",
    help="FASTA-formatted proteins file. If provided alongside --proteins-fn, vConTACT will start prior"
    " to PC generation.",
)
inputs.add_argument(
    "--rel-mode",
    type=str,
    choices=["BLASTP", "Diamond", "MMSeqs2"],
    default="Diamond",
    dest="rel_mode",
    help="Method to use to create the protein-protein similarity edge file. This "
    "is what the PC clustering is applied against.",
)
inputs.add_argument(
    "-b",
    "--blast-fp",
    type=str,
    dest="blast_fp",
    help="Blast results file in CSV or TSV format. Used for generating the protein clusters. If "
    "provided alongside --proteins-fn, vConTACT will start from the PC-generation step. If raw "
    "proteins are provided, THIS WILL BE SKIPPED. Reference DBs CANNOT BE USED IF THIS OPTION "
    "IS ENABLED!!",
)
inputs.add_argument(
    "-p",
    "--proteins-fp",
    type=str,
    default="test_data/VIRSorter_genomes_g2g.csv",
    dest="proteins_fp",
    help="A file linking the protein name (as in the blast file) and the genome names (csv or tsv)."
    " If provided alongside --blast-fp, vConTACT will start from the PC-generation step. If "
    "provided alongside --raw-proteins, vConTACT will start from creating the all-verses-all "
    "protein comparison and then generate protein clusters.",
)
inputs.add_argument(
    "--db",
    type=str,
    choices=["None", *ref_dbs],
    # default="ProkaryoticViralRefSeq85-ICTV",
    default="ProkaryoticViralRefSeq211-Merged",
    dest="db",
    help="Select a reference database to de novo cluster proteins against. If 'None' selected, be "
    "aware that there will be no references included in the analysis.",
)
inputs.add_argument(
    "--pcs-mode",
    type=str,
    choices=["ClusterONE", "MCL"],
    default="MCL",
    dest="pcs_mode",
    help="Whether to use ClusterONE or MCL for Protein Cluster (PC) generation.",
)
inputs.add_argument(
    "--vcs-mode",
    type=str,
    choices=["ClusterONE", "MCL"],
    default="ClusterONE",
    dest="vc_mode",
    help="Whether to use ClusterONE or MCL for Viral Cluster (VC) generation.",
)
inputs.add_argument(
    "--c1-bin",
    default="cluster_one-1.0.jar",
    type=str,
    dest="cluster_one",
    help="Location for clusterONE file. Path only used if vConTACT cant find in $PATH.",
)
inputs.add_argument(
    "--blastp-bin",
    type=str,
    dest="blastp_fp",
    help="Location for BLASTP executable. Path only used if vConTACT cant find in $PATH.",
)
inputs.add_argument(
    "--diamond-bin",
    type=str,
    dest="diamond_fp",
    help="Location for DIAMOND executable. Path only used if vConTACT cant find in $PATH.",
)
inputs.add_argument(
    "-o",
    "--output-dir",
    type=str,
    dest="output_dir",
    default="vContact_Output",
    help="Output directory",
)
inputs.add_argument(
    "-t",
    "--threads",
    type=int,
    dest="threads",
    default=psutil.cpu_count(logical=False),
    help="Number of CPUs to use. If nothing is provided, vConTACT will attempt to identify the number"
    " of available CPUs.",
)

pcs = parser.add_argument_group("Protein Cluster (PC) Arguments")
pcs.add_argument(
    "--pc-evalue",
    default="0.0001",
    type=float,
    dest="evalue",
    help="E-value used by BLASTP or Diamond when creating the protein-protein similarity network.",
)
pcs.add_argument(
    "--reported-alignments",
    default=25,
    type=int,
    dest="pc_alignments",
    help="Maximum number of target sequences per query to report alignments for in Diamond.",
)
pcs.add_argument(
    "--max-overlap",
    default=0.8,
    type=float,
    dest="pc_overlap",
    help="Specifies the maximum allowed overlap between two clusters. (ClusterONE only)",
)
pcs.add_argument(
    "--penalty",
    default=2.0,
    type=float,
    dest="pc_penalty",
    help="Sets a penalty value for the inclusion of each node. It can be used to model the possibility of "
    "uncharted connections for each node, so nodes with only a single weak connection to a cluster "
    "will not be added to the cluster as the penalty value will outweigh the benefits of adding the "
    "node. (ClusterONE only)",
)
pcs.add_argument(
    "--haircut",
    default=0.1,
    type=float,
    dest="pc_haircut",
    help="Apply a haircut transformation as a post-processing step on the detected clusters. A haircut "
    "transformation removes dangling nodes from a cluster. (ClusterONE only)",
)
pcs.add_argument(
    "--pc-inflation",
    default=2.0,
    type=float,
    dest="pc_inflation",
    help="Inflation value to define proteins clusters with MCL. (default: 2.0) (MCL only)",
)

vcs = parser.add_argument_group("Viral Cluster (VC) Arguments")
vcs.add_argument(
    "--vc-inflation",
    default=2.0,
    type=float,
    dest="vc_inflation",
    help="Inflation parameter to define contig clusters with MCL. (MCL only)",
)
vcs.add_argument(
    "--min-density",
    default=0.3,
    type=float,
    dest="vc_density",
    help="Sets the minimum density of predicted complexes. (ClusterONE only)",
)
vcs.add_argument(
    "--min-size",
    default=2,
    type=int,
    dest="vc_size",
    help="Minimum size for the Viral Cluster.",
)
vcs.add_argument(
    "--vc-overlap",
    default=0.9,
    type=float,
    dest="vc_overlap",
    help="Specifies the maximum allowed overlap between two clusters. (ClusterONE only)",
)
vcs.add_argument(
    "--vc-penalty",
    default=2.0,
    type=int,
    dest="vc_penalty",
    help="Sets a penalty value for the inclusion of each node. It can be used to model the possibility of "
    "uncharted connections for each node, so nodes with only a single weak connection to a cluster "
    "will not be added to the cluster as the penalty value will outweigh the benefits of adding the "
    "node. (ClusterONE only)",
)
vcs.add_argument(
    "--vc-haircut",
    default=0.55,
    type=float,
    dest="vc_haircut",
    help="Apply a haircut transformation as a post-processing step on the detected clusters. A haircut "
    "transformation removes dangling nodes from a cluster. (ClusterONE only)",
)
vcs.add_argument(
    "--merge-method",
    default="single",
    choices=["single", "multi"],
    dest="vc_method",
    help="Specifies the method to be used to merge highly overlapping complexes. (ClusterONE only)",
)
vcs.add_argument(
    "--similarity",
    default="match",
    choices=["match", "simpson", "jaccard", "dice"],
    dest="vc_similarity",
    help="Specifies the similarity function to be used in the merging step. (ClusterONE only)",
)
vcs.add_argument(
    "--seed-method",
    default="nodes",
    dest="vc_seed",
    choices=["unused_nodes", "nodes", "edges", "cliques"],
    help="Specifies the seed generation method to use. (ClusterONE only)",
)
vcs.add_argument(
    "--optimize",
    dest="optimize",
    action="store_true",
    help="Optimize hierarchical distances during second-pass of the viral clusters",
)

network = parser.add_argument_group("Similarity Network and Module Options")
network.add_argument(
    "--sig",
    default=1.0,
    type=float,
    help="Significance threshold in the contig similarity network.",
)
network.add_argument(
    "--max-sig",
    dest="max_sig",
    default=300,
    type=int,
    help="Significance threshold in the contig similarity network.",
)
network.add_argument(
    "--permissive",
    action="store_true",
    help="Use permissive affiliation for associating VCs with reference sequences.",
)
network.add_argument(
    "--mod-inflation",
    default=5.0,
    type=float,
    dest="mod_inflation",
    help="Inflation parameter to define protein modules with MCL.",
)
network.add_argument(
    "--mod-sig",
    default=1.0,
    type=float,
    dest="mod_sig",
    help="Significance threshold in the protein cluster similarity network.",
)
network.add_argument(
    "--mod-shared-min",
    default=3,
    type=int,
    dest="mod_shared_min",
    help="Minimal number (inclusive) of contigs a PC must appear into to be taken into account in the "
    "modules computing.",
)
network.add_argument(
    "--link-sig",
    default=1.0,
    type=float,
    dest="link_sig",
    help="Significitaivity threshold to link a cluster and a module",
)
network.add_argument(
    "--link-prop",
    default=0.5,
    type=float,
    dest="link_prop",
    help="Proportion of a module's PC a contig must have to be considered as displaying this module.",
)

outputs = parser.add_argument_group("Output Options")
outputs.add_argument(
    "-e",
    "--exports",
    nargs="*",
    default=["csv"],
    help='Export backend. Suported values are "csv", "krona" and "cytoscape"',
)
outputs.add_argument(
    "--cluster-filter",
    nargs="*",
    default=[0],
    help="Id of the clusters to export (Cytoscape export only).",
)
outputs.add_argument(
    "--criterion",
    default="predicted_family",
    help="Pooling criterion for cluster export (Cytoscape export only).",
)

misc = parser.add_argument_group("Misc. Options")
misc.add_argument(
    "-v",
    "--verbose",
    action="count",
    default=-2,
    help="Verbosity level : -v warning, -vv info, -vvv debug, (default info)",
)
misc.add_argument(
    "-f", "--force-overwrite", action="store_true", help="Overwrite existing files"
)


def init_logger(verbose: int):
    log_levels = [logging.WARNING, logging.INFO, logging.DEBUG]
    # create logger with 'spam_application'
    logger.setLevel(log_levels[verbose])
    ch = logging.StreamHandler()
    ch.setLevel(log_levels[verbose])
    ch.setFormatter(logging.Formatter("%(levelname)s:%(name)s: %(message)s"))
    logger.addHandler(ch)
    print("\n{:=^80}\n".format("This is vConTACT2 {}".format(vcontact2.__version__)))


def check_deps(args):
    # Checks
    cluster_one_fp = args.cluster_one
    cluster_one_java = False  # Need to keep track if it's jarfile or not
    if ("ClusterONE" in args.pcs_mode) or ("ClusterONE" in args.vc_mode):
        if not cluster_one_fp:  # If user hasn't provided a path
            cluster_one_fp = shutil.which("cluster_one-1.0.jar")
            if cluster_one_fp is None:  # Try other install method
                cluster_one_fp = shutil.which("clusterone")
                if cluster_one_fp is None:
                    logger.error("Could not find ClusterONE java file.")
                    raise FileNotFoundError("Could not find ClusterONE java file.")
            else:
                cluster_one_java = True
                logger.info("Found ClusterONE: {}".format(cluster_one_fp))

            # Ensure that ClusterONE is working
            if cluster_one_java:
                cluster_one_cmd = "java -jar {} --help".format(cluster_one_fp)
            else:
                cluster_one_cmd = "{} --help".format(cluster_one_fp)

            logger.info("Verifying clusterONE works...")
            res = subprocess.run(
                cluster_one_cmd,
                shell=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE,
            )
            if res.returncode == 0:
                logger.info("It appears that clusterONE is functional!")
            else:
                logger.error(
                    "It appears something went wrong with trying to run ClusterONE:\n{}".format(
                        res.stderr.decode("utf-8")
                    )
                )
                raise FileNotFoundError("Could not run ClusterONE.")

    diamond_fp = args.diamond_fp
    if args.rel_mode == "Diamond":
        if not diamond_fp:
            diamond_fp = shutil.which("diamond")
            if diamond_fp is None:
                logger.error("Could not find Diamond in $PATH.")
                raise FileNotFoundError("Could not find Diamond.")
            else:
                logger.info("Found Diamond: {}".format(diamond_fp))

    blastp_fp = args.blastp_fp
    if args.rel_mode == "BLASTP":
        if not blastp_fp:
            blastp_fp = shutil.which("blastp")
            if blastp_fp is None:
                logger.error("Could not find BLASTp in $PATH.")
                raise FileNotFoundError("Could not find BLASTp.")
            else:
                logger.info("Found BLASTP: {}".format(blastp_fp))

    if "MCL" in (args.pcs_mode or args.vc_mode):
        mcl_fp = shutil.which("mcxload")
        if shutil.which("mcxload") is None:
            logger.error("Could not find MCL (specifically mcxload) in $PATH.")
            raise FileNotFoundError("Could not find MCL.")
        else:
            logger.info("Found MCL: {}".format(mcl_fp))
    return cluster_one_fp


def read_dfs(contigs_fp, pcs_fp, profiles_fp, proteins_fp=None):
    """
    pcs_df = pos, id, pc_id, nb_proteins
    contigs_csv_df = pos, id (contig), proteins
    profiles = dict w/ matrix and singletons
    Refactor later
    """
    print("\n\n" + "{:-^80}".format("Loading data"))
    contigs_csv_df = pd.read_csv(contigs_fp)
    contigs_csv_df["contig_id"] = contigs_csv_df["contig_id"].str.replace(" ", "~")
    contigs_csv_df.index.name = "pos"
    contigs_csv_df.reset_index(inplace=True)
    #
    profiles_df = pd.read_csv(profiles_fp)
    profiles = profiles_df.copy()
    #
    # ClusterONE can't handle spaces
    profiles["contig_id"] = profiles["contig_id"].str.replace(" ", "~")
    # Filtering the PC profiles that appears only once
    before_filter = len(profiles)
    cont_by_pc = profiles.groupby("pc_id")["contig_id"].count().reset_index()
    #
    # get the number of contigs for each pcs and add it to the dataframe
    cont_by_pc.columns = ["pc_id", "nb_proteins"]
    pcs_csv_df = pd.merge(
        pd.read_csv(pcs_fp), cont_by_pc, left_on="pc_id", right_on="pc_id", how="left"
    ).fillna({"nb_proteins": 0})
    #
    # Drop the pcs that <= 1 contig from the profiles.
    pcs_csv_df = pcs_csv_df[pcs_csv_df["nb_proteins"] > 1].reset_index(drop=True)
    pcs_csv_df.index.name = "pos"
    pcs_csv_df = pcs_csv_df.reset_index()

    at_least_a_cont = cont_by_pc[
        cont_by_pc["nb_proteins"] > 1
    ]  # cont_by_pc.query("nb_contigs>1")
    # profiles = profiles.query("pc_id in at_least_a_cont.pc_id")
    profiles = profiles[profiles["pc_id"].isin(at_least_a_cont["pc_id"])]
    #
    logger.debug("Read {} entries from {}".format(len(contigs_csv_df), contigs_fp))
    logger.info(
        "Read {} entries (dropped {} singletons) from {}".format(
            len(profiles), (before_filter - len(profiles)), profiles_fp
        )
    )
    #
    profiles_matrix_singletons = vcontact2.pcprofiles.build_pc_matrices(
        profiles.copy(), contigs_csv_df, pcs_csv_df
    )
    return (contigs_csv_df, pcs_csv_df, profiles_df), profiles_matrix_singletons


def get_dfs(
    output_dir, cluster_one_fp, args=None
) -> tuple[
    tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame],
    tuple[
        vcontact2.pcprofiles.sparse.coo_matrix, vcontact2.pcprofiles.sparse.csr_matrix
    ],
]:
    fps = {
        f"{name}_fp": os.path.join(output_dir, "vConTACT_{}.csv".format(name))
        for name in ("contigs", "pcs", "profiles", "proteins")
    }
    if args is None or not args.force_overwrite:
        for fp in fps.values():
            if not os.path.exists(fp):
                assert args is not None
                break
        else:
            # No longer need to re-load dfs
            logger.info(f"Files {fps} exists and will be used. Use -f to overwrite.")
            return read_dfs(**fps)

    if args.db != "None":
        merged_aa_fp = os.path.join(output_dir, "merged.faa")

        if not os.path.exists(merged_aa_fp):
            logger.info("Merging {} to user sequences...".format(args.db))
            proteins_aa_fp = vcontact2.protein_clusters.merge_aa(
                args.raw_proteins, ref_dbs[args.db], merged_aa_fp
            )
        else:
            logger.info(
                f"Identified existing 'merged.faa' in output path: re-using {args.db}"
            )
            proteins_aa_fp = merged_aa_fp

    else:
        proteins_aa_fp = args.raw_proteins

    similarity_fp = ""
    if args.raw_proteins:
        if args.rel_mode == "BLASTP":
            blastp_out_fn = "{}.self-blastp.tab".format(
                os.path.basename(proteins_aa_fp).rsplit(".", 1)[0]
            )
            blastp_out_fp = os.path.join(output_dir, blastp_out_fn)

            if not os.path.exists(blastp_out_fn):
                logger.info("Creating BLAST database and running BLASTP...")
                db_fp = vcontact2.protein_clusters.make_blast_db(proteins_aa_fp)
                similarity_fp = vcontact2.protein_clusters.run_blastp(
                    proteins_aa_fp, db_fp, args.evalue, args.threads, blastp_out_fp
                )
            else:
                logger.info("Re-using existing BLASTP file...")
                similarity_fp = blastp_out_fp

        elif args.rel_mode == "Diamond":    # need CPU lock for better performance
            diamond_out_fn = "{}.self-diamond.tab".format(
                os.path.basename(proteins_aa_fp).rsplit(".", 1)[0]
            )
            diamond_out_fp = os.path.join(output_dir, diamond_out_fn)

            if not os.path.exists(diamond_out_fp):
                C_lock = acquire_lock()
                logger.info("Creating Diamond database and running Diamond...")
                db_fp = vcontact2.protein_clusters.make_diamond_db(
                    proteins_aa_fp, output_dir, args.threads
                )
                similarity_fp = vcontact2.protein_clusters.run_diamond(
                    proteins_aa_fp,
                    db_fp,
                    args.threads,
                    args.evalue,
                    args.pc_alignments,
                    diamond_out_fp,
                )
                release_lock(C_lock)
            else:
                logger.info("Re-using existing Diamond file...")
                similarity_fp = diamond_out_fp

        else:
            logger.error(
                "Unable to identify which method to use for generating protein-protein similarities."
            )

    # If they don't have raw proteins but did provide a BLASTP file
    elif (not args.raw_proteins) and args.blast_fp:
        logger.info(f"User provided BLASTp file {args.db}")
        similarity_fp = args.blast_fp
    else:
        logger.error(
            "User did not provide a proteins file (amino acid fasta file) OR a BLASTp file. One must be "
            "provided"
        )

    # Did somehow things get to this point?
    if not similarity_fp:
        logger.error("No similarity file identified?")

    print("\n\n" + "{:-^80}".format("Protein clustering"))

    # Load the proteins <-> contigs associations...
    # Should there be a "joining databases" block or just separate? (as in this case)
    logger.info("Loading proteins...")

    usr_gene2genome_fp = args.proteins_fp
    if os.path.basename(usr_gene2genome_fp).endswith(".csv"):
        usr_gene2genome_df = pd.read_csv(
            usr_gene2genome_fp, sep=",", header=0
        )  # Don't rely on pandas to identify
    elif os.path.basename(usr_gene2genome_fp).endswith(".tsv"):
        usr_gene2genome_df = pd.read_csv(usr_gene2genome_fp, sep="\t", header=0)
    else:
        logger.error(
            "Unable to identify the file format for the gene-to-genome file {}. Does it end with the "
            "extension *.tsv or *.csv?".format(os.path.basename(usr_gene2genome_fp))
        )
    if args.db != "None":
        logger.info("Merging {} to user gene-to-genome mapping...".format(args.db))

        ref_proteins_df = pd.read_csv(ref_g2c[args.db], sep=",", header=0)
        gene2genome_df = pd.concat([usr_gene2genome_df, ref_proteins_df])
    else:
        gene2genome_df = usr_gene2genome_df

    logger.debug(
        "Read {} proteins from {}.".format(len(gene2genome_df), args.proteins_fp)
    )

    pcs_mode = args.pcs_mode
    # pc_overlap, pc_penalty, pc_haircut, pc_inflation = 0.8, 2.0, 0.1, 2.0
    pc_overlap, pc_penalty, pc_haircut, pc_inflation = (
        args.pc_overlap,
        args.pc_penalty,
        args.pc_haircut,
        args.pc_inflation,
    )

    similarity_fn = os.path.basename(similarity_fp)
    if pcs_mode == "ClusterONE":
        pcs_fn = "{}_c1_{}_{}_{}.clusters".format(
            similarity_fn, pc_overlap, pc_penalty, pc_haircut
        )
    elif pcs_mode == "MCL":
        pcs_fn = "{}_mcl{}.clusters".format(similarity_fn, int(pc_inflation * 10))
    else:
        logger.error("A mode must be selected. Use ClusterONE or MCL to generate PCs.")

    pcs_fp = os.path.join(output_dir, pcs_fn)
    # Run clustering tool on the blast results...
    if not os.path.exists(pcs_fp) or args.force_overwrite:
        if pcs_mode == "MCL":
            pcs_fp = vcontact2.protein_clusters.make_protein_clusters_mcl(
                # similarity_fp, output_dir, args.pc_inflation, threads=args.threads
                similarity_fp, output_dir, args.pc_inflation
            )
        elif pcs_mode == "ClusterONE":
            pcs_fp = vcontact2.protein_clusters.make_protein_clusters_one(
                similarity_fp,
                cluster_one_fp,
                output_dir,
                args.pc_overlap,
                args.pc_penalty,
                args.pc_haircut,
            )
    else:
        logger.debug(
            "File {} exists and will be used. Use -f to overwrite.".format(pcs_fn)
        )

    # Load the clusters...
    logger.info(
        "Building the cluster and profiles (this may take some time...)\n"
        "If it fails, try re-running using --blast-fp flag "
        "and specifying merged.self-diamond.tab (or merged.self-blastp.tab)"
    )
    (
        protein_df,
        clusters_df,
        profiles_df,
        contigs_df,
    ) = vcontact2.protein_clusters.build_clusters(pcs_fp, gene2genome_df, pcs_mode)
    # protein_df = protein_id, contig_id, keywords, cluster
    # clusters_df = pc_id, size (#), annotated (#), keys
    # profiles_df = contig_id, pc_id
    # contigs_df = contig_id, proteins (#)

    # Export csv files...
    logger.info("Saving intermediate files...")  # Save the dataframes

    for name, df in {
        "proteins": protein_df,
        "contigs": contigs_df,
        "pcs": clusters_df,
    }.items():
        # hindsight
        df.set_index(name.strip("s") + "_id").to_csv(fps[f"{name}_fp"])

    profiles_df.to_csv(fps["profiles_fp"], index=False)
    return read_dfs(**fps)


def main(args):
    # There's almost no point in providing the ability to combine references with user data if they already have a
    # pre-run BLASTP (as is used with vConTACT-PCs. Otherwise the user would need to re-BLAST everything with the DB.

    # In essence, by allowing the opportunity to use "raw proteins" and have vConTACT handle everything, the
    # functionality provided by vConTACT-PCs is eliminated.
    init_logger(args.verbose)

    print("\n\n" + "{:-^80}".format("Pre-Analysis"))
    cluster_one_fp = check_deps(args)

    # Global
    logger.info("Identified {} CPUs".format(args.threads))
    logger.info("Using reference database: {}".format(args.db))

    if args.db == "None":
        logger.warning(
            "No reference database selected! vConTACT2 will be unable to run taxonomy-based performance "
            "metrics unless taxonomic information is included in the --proteins-fp file."
        )

    # Prepare output directory
    output_dir = args.output_dir
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
        logger.info("Output directory {} created.".format(output_dir))
    else:
        logger.info("Using existing directory {}.".format(output_dir))

    print("\n\n" + "{:-^80}".format("Reference databases"))

    # Merge potential reference databases to source

    # Test for presence of headers here, before Diamond and other seq analysis is done. Though out of order, better be
    # done before cycle investment
    if os.path.exists(args.proteins_fp):
        # Pull out into function?
        with open(args.proteins_fp) as gene2genome_fh:
            header = gene2genome_fh.readline().strip()

            if "contig_id" not in header:
                logger.error(
                    "Unable to identify header for gene-to-genome file. Please ensure that you have the"
                    ' following headers: "contig_id", "genome_id" and "keywords".'
                )
                raise KeyError(
                    "Unable to identify header for gene-to-genome file. Please ensure that you have the"
                    ' following headers: "contig_id", "genome_id" and "keywords".'
                )

    # Check if exists
    (contigs_csv_df, pcs_csv_df, profiles_df), profiles_matrix_singletons = get_dfs(
        output_dir, cluster_one_fp, args
    )  # read_dfs(**fps)

    # Loader
    merged_fp = os.path.join(output_dir, "merged_df.csv")
    # Protein Cluster Profile
    # .contigs, .pcs, .matrix, .singletons, .ntw, .ntw_modules
    name = "vConTACT2"
    pcp_fp = os.path.join(output_dir, name + ".pkle")

    if not os.path.exists(pcp_fp) or args.force_overwrite:
        taxonomies = ["Organism/Name", "origin", "order", "family", "genus"]
        extended_taxonomies = [
            "Organism/Name",
            "origin",
            "kingdom",
            "phylum",
            "order",
            "class",
            "family",
            "genus",
        ]
        if args.db != "None":
            print("\n\n" + "{:-^80}".format("Adding Taxonomy"))

            ref_df = pd.read_csv(ref_tax[args.db], header=0)
            if set(extended_taxonomies).issubset(ref_df.columns):
                ref_tax_ext = ref_df[extended_taxonomies]
            else:
                ref_tax_ext = ref_df[taxonomies]

            ref_tax_ext["Organism/Name"] = ref_tax_ext["Organism/Name"].str.replace(
                " ", "~"
            )  # Ensure c1 can process names

            merged_df = contigs_csv_df.merge(
                ref_tax_ext, how="left", left_on="contig_id", right_on="Organism/Name"
            ).drop("Organism/Name", axis="columns")

            merged_df.to_csv(merged_fp)
            # in any case, NO taxonomy assigned to user~given~contigs

            # Rid of pesky "Annotated" but not really
            merged_df.loc[
                merged_df["order"].isnull()
                & merged_df["family"].isnull()
                & merged_df["genus"].isnull(),
                "origin",
            ] = np.nan
        else:
            merged_df = contigs_csv_df
            merged_df.to_csv(merged_fp)

        print("\n\n" + "{:-^80}".format("Calculating Similarity Networks"))
        pcp = vcontact2.pcprofiles.PCProfiles(
            merged_df,
            pcs_csv_df,
            profiles_matrix_singletons,
            # args.threads,
            1,
            name,
            args.sig,  # 1.0
            args.max_sig,  # 300
            args.mod_sig,  # 1.0
            args.mod_shared_min,  # 3
        )
        if not args.force_overwrite:
            pcp.to_pickle(pcp_fp)
    else:
        logger.info(f"Re-using existing pickle {pcp_fp}...")
        pcp = vcontact2.pcprofiles.read_pickle(pcp_fp)

    # ntw = sparse matrix (contig x contig, edge weights)
    # pcs = pos, pc_id, size, annotated, keys, nb_proteins
    # contigs = pos, contig_id, proteins

    logger.debug("Network Contigs:\n {}".format(pcp.ntw.todense()))  # pcm contigs

    print("\n\n" + "{:-^80}".format("Contig Clustering & Affiliation"))
    try:
        cluster_one_args = {
            "--min-density": args.vc_density,
            "--min-size": args.vc_size,
            "--max-overlap": args.vc_overlap,
            "--penalty": args.vc_penalty,
            "--haircut": args.vc_haircut,
            "--merge-method": args.vc_method,
            "--similarity": args.vc_similarity,
            "--seed-method": args.vc_seed,
        }

        gc = vcontact2.contig_clusters.ContigCluster(
            pcp,
            output_dir,
            cluster_one_fp,
            cluster_one_args,
            # inflation=2.0, threshold=1.0, membership_simple=True, mode="ClusterONE")
            inflation=args.vc_inflation,
            threshold=args.sig,
            # permissive -> use abundance (default)
            membership_simple=not args.permissive,
            mode=args.vc_mode,
            # threads=args.threads,
            threads=1
        )

        # gc.pcs = pos, pc_id, size, annotated, keys, nb_proteins
        # gc.contigs = contig_id, index, pos, proteins, origin (if ref), order (ref), family (ref), genus (ref),
        # pos_cluster, membership, pos_cluster_mebship
        # pos.clusters = id, pos, size
        # gc.cluster_results = list of list cluster members

        # if len(gc.levels):  # Add predicted_{}
        #     quality = gc.total_affiliation()
        #     logger.info(quality)

    except Exception as e:
        logger.error("Error in contig clustering")
        raise e

    try:
        vc = vcontact2.cluster_refinements.ViralClusters(
            gc.contigs, profiles_df, optimize=args.optimize
        )
    except Exception as e:
        logger.error("Error in viral clusters")
        raise e

    print("\n\n" + "{:-^80}".format("Protein modules"))

    try:
        vm = vcontact2.modules.Modules(
            pcp,
            output_dir,
            # inflation=5.0, threshold=1.0, shared_min=3)
            inflation=args.mod_inflation,
            threshold=args.mod_sig,
            shared_min=args.mod_shared_min,
            # threads=args.threads,
            threads=1
        )

        # modules has saved module df, both pcs and and modules (2 files)
        # modules.matrix matrix
        # modules.matrix_module csc matrix

    except Exception as e:
        logger.error("Error in protein module computation")
        raise e

    # by the time the linking happens, contigs are already forced into clusters that they could be sharing....
    # VCs and PC modules are independent of each other, as each takes into account the totality of the other, yet isn't
    # influenced by either's result

    print("\n\n" + "{:-^80}".format("Link modules and clusters"))

    try:
        link = vm.link_modules_and_clusters_df(
            gc.clusters,
            gc.contigs,
            # thres=1.0, own_threshold=0.5)
            thres=args.link_sig,
            own_threshold=args.link_prop,
        )
    except Exception as e:
        link = None
        logger.error("Error in linking modules and clusters")
        raise e

    print("\n\n" + "{:-^80}".format("Exporting results files"))

    # OUTPUT
    if "csv" in args.exports:
        try:
            vcontact2.exports.csv.complete(output_dir, pcp, gc, vm, link)
        except Exception as e:
            logger.error("Error in CSV export")
            raise e

    if "profiles" in args.exports:
        try:
            vcontact2.exports.profiles.push()
        except Exception as e:
            logger.error("Error in profiles export")
            raise e

    if "cytoscape" in args.exports:
        cytoscape_file = output_dir + "cytoscape_export"

        try:
            vcontact2.exports.cytoscape.contigs(
                gc.network, gc.contigs, args.cluster_filter, cytoscape_file + "_contigs"
            )
        except Exception as e:
            logger.error("Error in Cytoscape export of contigs")
            raise e

        try:
            vcontact2.exports.cytoscape.clusters(
                gc.network,
                gc.contigs,
                gc.clusters,
                args.criterion,
                cytoscape_file + "_clusters",
            )

            vcontact2.exports.cytoscape.membership(
                cytoscape_file + "_membership",
                gc.matrix["B"],
                gc.contigs,
                gc.clusters,
                args.cluster_filter,
                args.criterion,
            )
        except Exception as e:
            logger.error("Error in Cytoscape export of clusters")
            raise e

    if "krona" in args.exports:
        krona_file = output_dir + "krona_export.txt"
        try:
            vcontact2.exports.krona.textfile(gc.contigs, krona_file)
        except Exception as e:
            logger.error("Error in Krona export")
            raise e

        try:
            subprocess.check_call("ktImportText {}".format(krona_file), shell=True)
        except Exception as e:
            logger.error(
                (
                    "ktImportText {} failed ({}). \n (Is krona executable in the PATH?). You can still do it "
                    "manually."
                ).format(krona_file, e)
            )

    # Final summary outputs

    # Collect all the genomes that failed to make it to this point in the analysis
    logger.info(
        "Identifying genomes that are not clustered (i.e. singletons, outliers and overlaps"
    )
    try:
        excluded = vcontact2.exports.summaries.find_excluded(
            merged_fp, os.path.join(output_dir, gc.name + ".ntw"), gc.df
        )
    except Exception as e:
        logger.error(
            f"Error in identifying excluded genomes (i.e. those dropped for being singletons or outliers):"
            f" {e}"
        )
        raise e

    if len(excluded) == 0:
        logger.warning(
            f"There were 0 genomes identified as singleton, outlier or overlaps. This is highly unusual "
            f"unless its a very simple dataset or only user sequences were provided. Carefully examine the "
            f"outputs."
        )

    print(
        f"There were {len(excluded)} genomes (including refs) that were singleton, outlier or overlaps."
    )

    # Final cleanup and adding back in those genomes
    logger.info("Building final summary table")
    try:
        vcontact2.exports.summaries.final_summary(
            output_dir,
            vc.contigs,
            os.path.join(output_dir, gc.name + ".ntw"),
            profiles_df,
            vc,
            excluded,
        )
    except Exception as e:
        logger.error(f"Error in exporting the final summary table: {e}")
        raise e


if __name__ == "__main__":
    options = parser.parse_args()

    # Logging config
    main(options)
