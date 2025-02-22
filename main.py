import os
import pandas
import json
import subprocess
import argparse
from Bio import SeqIO

def translate(inputFile, outputFile, g2gFile):
    # from DNA sequence to protein
    # DNAfasta = "/Data/ICTVData/Challenge.fasta"
    # DNAfasta = "/Data/ICTVData/ChallengeData/ICTVTaxoChallenge_10163.fasta"
    # proteinFasta = "/Data/ICTVData/ProteinFasta/Challenge.fasta"
    # proteinFasta = "/Data/ICTVData/ProteinFasta/ChallengeData/ICTVTaxoChallenge_10163.fasta"
    command = "prodigal-gv -i " + inputFile + " -a " + outputFile + ".raw -p meta -q"
    # print("translating DNA to protein sequences")
    # os.system(command)
    subprocess.run(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    focusedSamples = set()

    # then remove comments in the raw file
    with open(outputFile, 'wt') as fp:
        with open(outputFile + ".raw") as rawFP:
            for line in rawFP:
                if line.startswith('>'):
                    fp.write(line.strip().split(' ')[0] + "\n")
                else:
                    fp.write(line)

    lines = list()
    lines.append("protein_id,contig_id,keywords\n")
    seqs = SeqIO.parse(outputFile, 'fasta')
    for seq in seqs:
        contigName = seq.id[:seq.id.rfind('_')]
        lines.append(f"{seq.id},{contigName},na\n")
        focusedSamples.add(contigName)
    with open(g2gFile, 'wt') as fp:
        fp.writelines(lines)

    return focusedSamples


def extractOutput(outputFolder, focusedSamples):
    resultFile = f"{outputFolder}/viral_cluster_overview.csv"
    summaryFile = f"{outputFolder}/summary.json"
    results = dict()
    df = pandas.read_csv(resultFile)
    ranks = ["kingdom", "phylum", "class", "order", "family", "genus"]
    for idx, row in df.iterrows():
        samples = row.Members.split(',')
        result = "N/A"
        for rank in reversed(ranks):  # try genus first, if the result is confirmed then no need to check family
            candidates = json.loads(row[rank].replace("'", '"'))
            if ("Unassigned" in candidates):
                candidates.remove("Unassigned")
            if (len(candidates) == 1):
                result = candidates[0]
                break
        for sample in samples:
            if (sample in focusedSamples):
                if (sample not in results):
                    results[sample] = result
                elif results[sample] != result:
                    print(f"{sample} has multiple results")
    with open(summaryFile, 'wt') as fp:
        json.dump(results, fp, indent=2, sort_keys=True)

def main(args):
    # proteinFasta = "/Data/ICTVData/ProteinFasta/Challenge.fasta"
    # proteinFasta = "/Data/ICTVData/ProteinFasta/ChallengeData/ICTVTaxoChallenge_10163.fasta"
    os.makedirs(args.output, exist_ok=True)

    proteinFasta = f"{args.output}/protein.fasta"
    g2gFile = f"{args.output}/g2g.csv"
    focusedSamples = translate(args.input, proteinFasta, g2gFile)

    command = "conda run -n vContact2 python bin/v2.py --raw-proteins " + proteinFasta + \
          " --proteins-fp " + g2gFile + \
          " --db 'ProkaryoticViralRefSeq211-Merged' " + \
          " --output-dir " + args.output
        #   " --c1-bin ./cluster_one-1.0.jar " + \
        #   " --rel-mode Diamond " + \
        #   " --threads 16 " + \
    # os.system(command)

    subprocess.run(command, shell=True)

    extractOutput(args.output, focusedSamples)

# translate("test_data/test.fasta", "test_data/test_protein.fasta", "test_data/test_g2g.csv")
# main()

if (__name__ == '__main__'):
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, required=True)
    parser.add_argument('--output', type=str, required=True)
    args = parser.parse_args()
    main(args)