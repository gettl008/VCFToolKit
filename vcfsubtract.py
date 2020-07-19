#!/usr/bin/env python2

import argparse
import os
import sys
import re

###############################################


def subtractVcfs(vcf1, vcf2, strict):
    outlines = [line.rstrip() for line in vcf1 if line.startswith('##')]
    header = [s for s in vcf1 if '#CHROM' in s][0].split('\t')[:9]
    samples1 = [s for s in vcf1 if '#CHROM' in s][0].split('\t')[9:]
    samples1[-1] = samples1[-1].rstrip()
    header = header + samples1
    outlines.append('\t'.join(header))
    for line1 in vcf1:
        if not line1.startswith('#'):
            line1 = line1.split('\t')
            line1[-1] = line1[-1].rstrip()
            match = '\t'.join(line1[0:2]) + '\t'
            line2 = [l.rstrip() for l in vcf2 if l.startswith(match)]
            if len(line2) > 0:
                line2 = line2[0].split('\t')
                line2[-1] = line2[-1].rstrip()
                if line1[4] != line2[4]:
                    alt1 = line1[4].split(',')
                    alt2 = line2[4].split(',')
                    overlap = False
                    nonoverlap = False
                    for a in alt1:
                        if a in alt2:
                            overlap = True
                        else:
                            nonoverlap = True
                    if overlap and nonoverlap and not strict:
                        outlines.append('\t'.join(line1))
                    elif not overlap:
                        outlines.append('\t'.join(line1))
            else:
                outlines.append('\t'.join(line1))
    return outlines


def main():
    parser = argparse.ArgumentParser(description='Returns all values unique to vcfA (not in vcfB)')
    parser.add_argument('vcfA', metavar='vcfA', type=str, nargs=1,
                        help='VCF file')
    parser.add_argument('vcfB', metavar='vcfB', type=str, nargs='+',
                        help='VCF file')
    parser.add_argument('-o', '--output_file',
                        help='Output VCF file. File will be overwritten if it already exists.')
    parser.add_argument('-s', '--strict',
                        help='Remove entries even if only partially overlap \
                             (i.e. share one or more of multiple alternative alleles).', action="store_true")
    args = parser.parse_args()

    try:
        args.output_file
        if os.path.exists(args.output_file):
            os.remove(args.output_file)
        try:
            out = open(args.output_file, "w+")
        except:
            print "Could not open file %s" % args.output_file
            sys.exit(1)
    except:
        pass

    try:
        vcf1 = open(args.vcfA[0], "r").readlines()
    except:
        print "Could not open vcf %s" % args.vcfA[0]
        sys.exit(1)

    while len(args.vcfB) >= 1:
        try:
            vcf2 = open(args.vcfB[0], "r").readlines()
        except:
            print "Could not open vcf %s" % args.vcfB[0]
            sys.exit(1)
        del args.vcfB[0]
        vcf1 = subtractVcfs(vcf1, vcf2, args.strict)
    try:
        out.write('\n'.join(vcf1))
    except:
        print '\n'.join(vcf1)


main()