#!/usr/bin/env python3

# Copyright (c) 2025 Outbreak Preparedness and Response Group at BNITM
# This file is part of ViMOP and is licensed under the MIT License.
# See the LICENSE file in the root of this repository for full license details.


import argparse
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_args():
    parser = argparse.ArgumentParser(
        description='Build a consensus FASTA from a Medaka gVCF, masking low-quality and low-coverage sites.'
    )
    parser.add_argument(
        '--ref',
        required=True,
        help='Reference fasta file.',
    )
    parser.add_argument(
        '--gvcf',
        required=True,
        help='Input gVCF file (Medaka output, annotated with DP).',
    )
    parser.add_argument(
        '--min_depth',
        type=int,
        default=20,
        help='Minimum read depth (DP) to keep a site unmasked. Default: 20.',
    )
    parser.add_argument(
        '--min_qual',
        type=int,
        default=20,
        help='Minimum variant call quality (QUAL) to keep a site unmasked. Default: 20.',
    )
    parser.add_argument(
        '--model',
        default='Unknown',
        help='Model used for medaka (only for output file info)',
    )
    parser.add_argument(
        '--sample',
        default='Unknown',
        help='Samle name (only for output file info)',
    )
    parser.add_argument(
        '--out',
        default='consensus.fasta',
        help='Output consensus FASTA path. Default: consensus.fa',
    )
    return parser.parse_args()


def read_ref(fname):
    return SeqIO.read(fname, 'fasta')


def read_vcf(fname):
    return list(pysam.VariantFile(fname))


def filter_variants(variants, min_qual, min_depth):
    return [
        v
        for v in variants
        if v.info.get("DP", 0) >= min_depth and v.qual is not None and v.qual >= min_qual
    ]


def build_consensus(variants, reflen):
    ref_called = {
        v.pos - 1: v.ref
        for v in variants
        if v.alts is None
    }
    var = {
        v.pos - 1: (len(v.ref), v.alts[0])
        for v in variants
        if v.alts is not None
    }
    cons = []
    i = 0
    while i < reflen:
        if i in var:
            varlen, alt = var[i]
            cons.append(alt)
            i += varlen
        elif i in ref_called:
            cons.append(ref_called[i])
            i += 1
        else:
            cons.append('N')
            i += 1
    return ''.join(cons)


def write_consensus(fname, consensus, ref_id, model, sampleid):
    description = f'method=medaka reference={ref_id} sample={sampleid} medaka_model={model}'
    record = SeqRecord(
        id='consensus',
        name='consensus',
        seq=Seq(consensus),
        description=description
    )
    SeqIO.write(record, fname, 'fasta')


def main():

    args = parse_args()
    ref = read_ref(args.ref)
    variants = read_vcf(args.gvcf)
    reflen = len(ref.seq)
    filtered = filter_variants(variants, args.min_qual, args.min_depth)
    consensus = build_consensus(filtered, reflen)
    write_consensus(args.out, consensus, ref.id, args.model, args.sample)


if __name__ == '__main__':
    main()
