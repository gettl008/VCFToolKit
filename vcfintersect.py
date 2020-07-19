#!/usr/bin/env python2

import argparse
import os
import sys
import re

###############################################


def findIntersections(vcf1, vcf2):
    outlines = [line.rstrip() for line in vcf1 if line.startswith('##')]
    header = [s for s in vcf1 if '#CHROM' in s][0].split('\t')[:9]
    samples1 = [s for s in vcf1 if '#CHROM' in s][0].split('\t')[9:]
    samples1[-1] = samples1[-1].rstrip()
    header = header + samples1
    samples2 = [s for s in vcf2 if '#CHROM' in s][0].split('\t')
    samples2[-1] = samples2[-1].rstrip()
    new_samples2 = []
    for s in samples2[9:]:
        if s not in samples1:
            header.append(s)
            new_samples2.append(samples2.index(s))
    outlines.append('\t'.join(header))
    for line1 in vcf1:
        if not line1.startswith('#'):
            line1 = line1.split('\t')
            line1[-1] = line1[-1].rstrip()
            match = '\t'.join(line1[0:5]) + '\t'
            line2 = [l.rstrip() for l in vcf2 if l.startswith(match)]
            if len(line2) > 0:
                line2 = line2[0].split('\t')
                line2[-1] = line2[-1].rstrip()
                if line1[4] != line2[4]:
                    alt1 = line1[4].split(',')
                    alt2 = line2[4].split(',')
                    alt_all = alt1[:]
                    alt2_dict = {0: 0}
                    for i, a in enumerate(alt2):
                        if a not in alt1:
                            alt_all.append(a)
                        alt2_dict[i + 1] = alt_all.index(a) + 1
                    if len(alt_all) < len(alt1 + alt2):
                        outline = line1[:4]
                        outline.append(','.join(alt_all))
                        outline = outline + line1[5:]
                        for s in new_samples2:
                            samp_data = line2[s].split(':')
                            genotype = re.split("[/|]", samp_data[0])
                            for i, haplotype in enumerate(genotype):
                                if haplotype == ".":
                                    haplotype = 0
                                genotype[i] = str(alt2_dict[int(haplotype)])
                            genotype.sort()
                            genotype = '/'.join(genotype)
                            samp_data = [genotype] + samp_data[1:]
                            outline.append(':'.join(samp_data))
                        outlines.append('\t'.join(outline))
                else:
                    outline = line1
                    for s in new_samples2:
                        outline.append(line2[s])
                    outlines.append('\t'.join(outline))
    return outlines


def main():
    parser = argparse.ArgumentParser(description='Generates a VCF of overlapping variants from a set of VCF files.')
    parser.add_argument('vcfs', metavar='vcf', type=str, nargs='+',
                        help='VCF file')
    parser.add_argument('-o', '--output_file',
                        help='Output VCF file. File will be overwritten if it already exists.')
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
        vcf1 = open(args.vcfs[0], "r").readlines()
    except:
        print "Could not open vcf %s" % args.vcfs[0]
        sys.exit(1)
    del args.vcfs[0]
    if len(args.vcfs) == 0:
        outlines = [line.rstrip() for line in vcf1]
        vcf1 = outlines
    while len(args.vcfs) >= 1:
        try:
            vcf2 = open(args.vcfs[0], "r").readlines()
        except:
            print "Could not open vcf %s" % args.vcfs[0]
            sys.exit(1)
        del args.vcfs[0]
        vcf1 = findIntersections(vcf1, vcf2)

    try:
        out.write('\n'.join(vcf1))
    except:
        print '\n'.join(vcf1)


main()