#!/usr/bin/env python
"""
Command line tool for distance 2 self calculation


"""
import sys
import argparse
import csv
import logging
import os
import pandas

import itertools as itr
from Distance2SelfBinding import Distance2Self
from DistanceMatrix import DistanceMatrix

from Fred2.Core import Allele


def read_hla_input(input, hla_header):
    """
    reads in the hla file
    header are defined as:

    :param hla_file:
    :return: list(Allele)
    """
    return map(Allele, set(pandas.DataFrame.from_csv(input, sep="\t", index_col=False)[hla_header]))


def load_blossum(blos):
    """
    loads a BLOSUM matrix

    :param str blos: Specifeis the BLOSUm matrix to lead
    :return: dict(str, dict(str, float)) - A BLOSUM1 matrix
    """
    try:
        mod = __import__('DistanceMatrices', fromlist=[blos])
        return getattr(mod, blos)
    except:
        mod = __import__('DistanceMatrices', fromlist=["BLOSUM50_distances"])
        return DistanceMatrix(getattr(mod, "BLOSUM50_distances"))

def main():
    parser = argparse.ArgumentParser(
        description="Distance to self calculation",

    )

    subparsers = parser.add_subparsers(help='Distance2Self offers two sub-command', dest="sub_command")
    parser_gen = subparsers.add_parser('generate',
                                       help='Command lets you generate an distance trie based on a provided peptide list')

    parser_gen.add_argument("-i", "--input",
                        required=True,
                        type=str,
                        help="Peptide with immunogenicity file (from epitopeprediction)",
                        )

    parser_gen.add_argument("-s", "--sequence",
                        required=False,
                        default="neopeptide",
                        type=str,
                        help="The columns name of the peptide sequences",
                        )

    parser_gen.add_argument("-o", "--output",
                        required=True,
                        type=str,
                        help="Specifies the output path. Results will be written to CSV",
                        )

    parser_gen.add_argument("-b", "--blosum",
                        required=False,
                        default="BLOSUM50",
                        type=str,
                        help="Specifies BLOSUM distance matrix (default BLOSUM50; available BLOSUM45, BLOSUM90)",
                        )

    #Prediction sub-command
    parser_pred = subparsers.add_parser('predict',
                                        help='Command calculates the distance to self for a provided list of peptides')
    parser_pred.add_argument("-t", "--trie",
                             required=False,
                             default=None,
                             type=str,
                             help="Specifies a custom distance trie to use"
                             )

    parser_pred.add_argument("-s", "--sequence",
                        required=False,
                        default="neopeptide",
                        type=str,
                        help="The columns name of the peptide sequences",
                        )

    parser_pred.add_argument("-k", "--k",
                             required=False,
                             default=1,
                             type=int,
                             help="Specifies the number of closest self-peptides to find"
                             )

    parser_pred.add_argument("-b", "--blosum",
                        required=False,
                        default="BLOSUM50",
                        type=str,
                        help="Specifies BLOSUm distance matrix (default BLOSUM50; available BLOSUM45, BLOSUM90)",
                        )

    parser_pred.add_argument("-a", "--alleles",
                             required=False,
                             default="HLA",
                             type=str,
                             help="Specifies the HLA allele column header of the peptide input file",
                             )

    parser_pred.add_argument("-i", "--input",
                        required=True,
                        type=str,
                        help="Peptide with immunogenicity file (from epitopeprediction)",
                        )

    parser_pred.add_argument("-o", "--output",
                        required=True,
                        type=str,
                        help="Specifies the output path. Results will be written to CSV",
                        )

    args = parser.parse_args()
    blos = load_blossum("{blos}_distance".format(blos=args.blosum.strip().upper()))
    dist2self = Distance2Self(blos,saveTrieFile=True)
    df = pandas.DataFrame.from_csv(args.input, sep="\t", index_col=False)
    peps = list(set(df[args.sequence]))

    if args.sub_command == "generate":
        peps.sort(key=len)
        for plength, peps in itr.groupby(peps, key=len):
            dist2self.generate_trie(peps, peptideLength=plength, outfile="{path}_l{peplength}.trie".format(
                                                                  path=os.path.splitext(args.output)[0],
                                                                  peplength=plength))

    else:
        peps.sort(key=len)
        for plength, peps in itr.groupby(peps, key=len):
            alleles = read_hla_input(args.input, args.alleles)
            pathToTrie = args.trie if args.trie is not None and os.path.isfile(args.trie) else None
            res = dist2self.calculate_distances(peps, alleles=alleles, hla_header=args.alleles, pep_header=args.sequence,
                                                pathToTrie=pathToTrie, n=args.k)
            merged = pandas.merge(df, res, how="outer",on=[args.sequence,args.alleles])
            merged.to_csv(args.output, sep="\t",index=False)

if __name__ == "__main__":
    sys.exit(main())

