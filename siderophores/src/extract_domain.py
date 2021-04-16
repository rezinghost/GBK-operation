# Copyright (c) 2015 Xiaowen Lu, modified based on 'extract_A_domains.py'


# Libraries to import
import sys
import os
import string
import subprocess
import argparse
import re


# ----------------------------------
# basic functions
# -----------------------------------
def sort_XonY(x, y, r):
    """
    This is the function to order value in list x according to the ordering in list y, where y_i corresponds to x_i
    :param x: the list needed to be sorted
    :param y: a list with the same length of x, which can be ranked
    :param r: True--order x when y is ordered decreasingly; False--order x when y is ordered increasingly
    :return: a list, which is sorted x
    """
    if r == False:
        x_sorted = [x for (y, x) in sorted(zip(y, x), key=lambda pair: pair[0])]
    elif r == True:
        x_sorted = [x for (y, x) in sorted(zip(y, x), reverse=True, key=lambda pair: pair[0])]
    return x_sorted


try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO
import warnings

with warnings.catch_warnings():  # Don't display the SearchIO experimental warning.
    warnings.simplefilter("ignore")
    from Bio import SearchIO

# Global settings
global cpu
##Number of CPU cores used for multiprocessing
cpu = 4


def execute(commands, input=None):
    "Execute commands in a system-independent manner"

    if input is not None:
        stdin_redir = subprocess.PIPE
    else:
        stdin_redir = None

    proc = subprocess.Popen(commands, stdin=stdin_redir,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    try:
        out, err = proc.communicate(input=input)
        retcode = proc.returncode
        return out, err, retcode
    except OSError as e:
        print("%r %r returned %r" % (commands, input[:40], e))
        raise


def run_hmmsearch(hmm_Version, query_hmmfile, fasta_file, cut_off, hmmoutput):
    "Run hmmsearch"
    if hmm_Version == 'v3':
        print(
        "Goes for hmmer3!"
        format = 'hmmsearch3-domtab')
        if cut_off == 'tc':
            command = ["hmmsearch", "--cpu", "4", "--cut_tc", "--domtblout", hmmoutput, query_hmmfile, fasta_file]
        elif cut_off == 'none':
            command = ["hmmsearch", "--cpu", "4", "--domtblout", hmmoutput, query_hmmfile, fasta_file]
        else:
            command = ["hmmsearch", "--cpu", "4", "--domE", cut_off, "--domtblout", hmmoutput, query_hmmfile,
                       fasta_file]
    else:
        format = 'hmmer2-text'
        print
        "Goes for hmmer2!"
        if cut_off == 'tc':
            command = ["../external/hmmer-2.3.2/src/hmmsearch", "--cut_tc", query_hmmfile, fasta_file]
        elif cut_off == 'none':
            command = ["../external/hmmer-2.3.2/src/hmmsearch", query_hmmfile, fasta_file]
        else:
            command = ["../external/hmmer-2.3.2/src/hmmsearch", "--domE", cut_off, query_hmmfile, fasta_file]
    try:
        out, err, retcode = execute(command)
        if hmm_Version == 'v3':
            out = open(hmmoutput, "r").read()
            res_stream = StringIO(out)
            results = list(SearchIO.parse(res_stream, format))
        else:
            # out = open("/Users/Xiaowen/Documents/trans_AT/temp_hmm2_1.txt", 'r').read()
            out = out.split("Alignments of top-scoring domains:")[0]
            out = out.split("Parsed for domains:")[1]
            results = out
        return results
    except OSError:
        return []


def read_fasta_file(fastafile):
    "Get sequences from fasta file and store in dictionary"
    ###My idea:
    infile = open(fastafile).read()
    entries = infile.split(">")[1:]
    # print entries
    fastadict = {}  # Accession as key, sequence as value
    for entry in entries:
        accession = entry.partition("\n")[0].replace("\r", "")
        # accession += "|"
        # print accession
        sequence = entry.partition("\n")[2].partition("\n")[0]
        # print sequence
        fastadict[accession] = sequence
    return fastadict


def run_HMMer(hmmversion, hmmfile, fastafile, cut_off, hmmoutput):
    "Run HMMer to find start and end sites of A domains"
    runresults = run_hmmsearch(hmm_Version=hmmversion, query_hmmfile=hmmfile, fasta_file=fastafile, cut_off=cut_off,
                               hmmoutput=hmmoutput)
    if hmmversion == 'v3':
        results_by_id = {}
        for runresult in runresults:
            # Store results in dictionary by NRPS accession
            for hsp in runresult.hsps:
                if not results_by_id.has_key(hsp.hit_id):
                    results_by_id[hsp.hit_id] = [hsp]
                else:
                    results_by_id[hsp.hit_id].append(hsp)
    else:
        results_by_id = {}
        runresults_1 = runresults.split('\n')[3:]
        runresults_1 = [l for l in runresults_1 if l != '']
        if '\t[no hits above thresholds]' not in runresults_1:
            id_list = []
            loc_list = []
            for r in runresults_1:
                runresult = [i for i in r.split(' ') if i != '']
                id = runresult[0]
                start = runresult[2]
                end = runresult[3]
                id_list.append(id)
                loc_list.append([start, end])
            for x, y in zip(id_list, loc_list):
                results_by_id.setdefault(x, []).append(y)
        else:
            results_by_id = {}
    return results_by_id


def write_fasta(fastadict, outfile, hmmer_results, domain_ab):
    # For each HMM, print a FASTA file with all domain hits
    out_file = open(outfile, "w")
    domnr = 1
    for cds in hmmer_results.keys():
        # Get sequence from fastadicts
        ###TO BE DONE
        cds_sequence = fastadict[cds]
        ###TO BE DONE
        cds_name = re.sub(r'(\:|\'|\(|\)|\,|\?|\<|\>\;)', '', cds)

        # For each hit, write out domain name + sequence in FASTA format
        for hit in hmmer_results[cds]:
            domain_sequence = cds_sequence[hit.hit_start:hit.hit_end]
            loc = '-'.join((str(hit.hit_start), str(hit.hit_end)))
            domain_name = "%s|%s|%s%s" % (cds_name, loc, domain_ab, domnr)
            out_file.write(">%s\n%s\n" % (domain_name, domain_sequence))
            domnr += 1
        domnr = 1
    out_file.close()


def write_fasta_hmm2(fastadict, outfile, hmmer_results, domain_ab):
    # For each HMM, print a FASTA file with all domain hits
    out_file = open(outfile, "w")
    faa_keys = fastadict.keys()
    if len(hmmer_results) > 0:
        for cds in hmmer_results.keys():
            # Get sequence from fastadicts

            cds_name1 = [k for k in faa_keys if cds in k][0]
            cds_name = re.sub(r'(\:|\'|\(|\)|\,|\?|\<|\>)', '', cds_name1)
            cds_sequence = fastadict[cds_name1]
            # For each hit, write out domain name + sequence in FASTA format
            hit_start_sort = []
            domain_name_sort = []
            domain_sequence_sort = []

            for hit in hmmer_results[cds]:
                hit_start_sort.append(int(hit[0]))
                domain_sequence = cds_sequence[(int(hit[0]) - 1):int(hit[
                                                                         1])]  # since the string are index from zero, so we need to get the first aa of the domain by int(hit[0])-1
                # loc = '-'.join(hit)
                # loc = '_'
                domain_name_sort.append("%s_%s" % (cds_name, domain_ab))
                domain_sequence_sort.append(domain_sequence)
            domain_name_sorted = sort_XonY(domain_name_sort, hit_start_sort, r=False)
            domain_sequence_sorted = sort_XonY(domain_sequence_sort, hit_start_sort, r=False)
            index = range(1, len(hit_start_sort) + 1)

            for name, nr, seq in zip(domain_name_sorted, index, domain_sequence_sorted):
                out_file.write(">%s%s\n%s\n" % (name, nr, seq))
        out_file.close()
    else:
        out_file.close()


if __name__ == "__main__":
    "Run this file from here as a script"
    # Check if parameters are provided; if not, exit with explanation
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument('--hmmfile')
    parser.add_argument('--fastafile')
    parser.add_argument('--outfile')
    parser.add_argument('--domain_ab')
    parser.add_argument('--hmmversion')
    parser.add_argument('--cut_off', action='store', default="0.01")

    parser.add_argument('--hmmoutput')

    args = parser.parse_args()

    if any(t == [] for t in vars(args).values()) == False:
        print
        args
        hmmfile = args.hmmfile
        fastafile = args.fastafile
        outfile = args.outfile
        domain_ab = args.domain_ab
        hmmversion = args.hmmversion
        cut_off = args.cut_off
        hmmoutput = args.hmmoutput

        fastadict = read_fasta_file(fastafile)
        hmmer_results = run_HMMer(hmmversion, hmmfile, fastafile, cut_off, hmmoutput)

        if hmmversion == 'v3':
            write_fasta(fastadict, outfile, hmmer_results, domain_ab)
        else:
            write_fasta_hmm2(fastadict, outfile, hmmer_results, domain_ab)
    else:
        print
        """Please provide needed profiles"""
        sys.exit()