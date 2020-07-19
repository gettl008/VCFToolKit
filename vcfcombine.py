#!/usr/bin/env python2

import sys
import argparse
import re
import os
import sys
import numpy
if sys.version_info[0] < 3:
    from StringIO import StringIO
else:
    from io import StringIO
import pandas as pd


###############################################

def findUnion(vcf1, vcf2, combine):
    outlines = [line.rstrip() for line in vcf1 if line.startswith('##')]
    samples1 = [s for s in vcf1 if '#CHROM' in s][0].split('\t')
    samples1[-1] = samples1[-1].rstrip()
    header = samples1[:]
    samples2 = [s for s in vcf2 if '#CHROM' in s][0].split('\t')
    samples2[-1] = samples2[-1].rstrip()
    new_samples2 = []
    old_samples2 = {}
    for s in samples2[9:]:
        if combine:
            if s not in samples1:
                header.append(s)
                new_samples2.append(samples2.index(s))
            else:
                index1 = samples1.index(s)
                index2 = samples2.index(s)
                old_samples2[index2] = index1
        else:
            new_samples2.append(samples2.index(s))
            if s in samples1[9:]:
                s = "%s.X" % s
            header.append(s)
    outlines.append('\t'.join(header))
    for line1 in vcf1:
        if not line1.startswith('#'):
            # Cycle lines of first VCF
            line1 = line1.split('\t')
            line1[-1] = line1[-1].rstrip()
            match = '\t'.join(line1[0:2]) + '\t'
            # If matching position in second VCF try to reconcile
            line2 = [l.rstrip() for l in vcf2 if l.startswith(match)]
            if len(line2) > 0:
                line2 = line2[0].split('\t')
                line2[-1] = line2[-1].rstrip()
                # Reconcile differences in alt alleles
                if line1[4] != line2[4]:
                    alt1 = line1[4].split(',')
                    alt2 = line2[4].split(',')
                    alt_all = alt1[:]
                    alt2_dict = {0: 0}
                    for i, a in enumerate(alt2):
                        if a not in alt1:
                            alt_all.append(a)
                        alt2_dict[i + 1] = alt_all.index(a) + 1
                    outline = line1[:4]
                    outline.append(','.join(alt_all))
                    outline = outline + line1[5:]
                    for s in new_samples2:
                        samp_data = line2[s].split(':')
                        genotype = re.split("[/|]", samp_data[0])
                        for i, haplotype in enumerate(genotype):
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
            else:
                outline = line1
                field_info_len = len(line1[8].split(':')) - 1
                new_sampleinfo = ':'.join(['.'] * field_info_len)
                new_sampleinfo = "0/0:" + new_sampleinfo
                for s in new_samples2:
                    outline.append(new_sampleinfo)
                outlines.append('\t'.join(outline))
    count = 0
    for line1 in vcf2:
        if not line1.startswith('#'):
            line1 = line1.split('\t')
            line1[-1] = line1[-1].rstrip()
            match = '\t'.join(line1[0:2])
            line2 = [l.rstrip() for l in vcf1 if l.startswith(match)]
            if len(line2) == 0:
                count = count + 1
                outline = line1[:9]
                field_info_len = len(line1[8].split(':')) - 1
                new_sampleinfo = ':'.join(['.'] * field_info_len)
                new_sampleinfo = "0/0:" + new_sampleinfo
                for s in samples1[9:]:
                    outline.append(new_sampleinfo)
                for s in new_samples2:
                    outline.append(line1[s])
                for s in old_samples2.keys():
                    outline[old_samples2[s]] = line1[s]
                outlines.append('\t'.join(outline))
    return outlines


def sortVcf(vcf):
    outlines = [line.rstrip() for line in vcf if line.startswith('#')]
    contig_lines = [line.rstrip() for line in vcf if line.startswith('##contig=<ID=')]
    if len(contig_lines) == 0:
        print "Contig order not provided in header. Returning unsorted VCF."
        return vcf
    else:
        contigs = []
        for line in contig_lines:
            contigs.append(re.split('[=,]', line)[2])
        for contig in contigs:
            match = contig + '\t'
            contigmatches = [line.rstrip() for line in vcf if line.startswith(match)]
            if len(contigmatches) > 0:
                contigdata = StringIO('\n'.join(contigmatches))
                pddata = pd.read_table(contigdata, header=None)
                pddata.sort_values(by=1)
                for line in pddata.values:
                    linelist = [str(i) for i in line.tolist()]
                    outlines.append('\t'.join(linelist))
    return outlines


def main():
    parser = argparse.ArgumentParser(description='Returns a VCF of all unique variant in all input VCFs.')
    parser.add_argument('vcfs', metavar='vcf', type=str, nargs='+',
                        help='VCF file')
    parser.add_argument('-o', '--output_file',
                        help='Output VCF file. File will be overwritten if it already exists.')
    parser.add_argument('-c', '--combine',
                        help='Combine samples with same name into one column.', action="store_true")
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

    while len(args.vcfs) >= 1:
        try:
            vcf2 = open(args.vcfs[0], "r").readlines()
        except:
            print "Could not open vcf %s" % args.vcfs[0]
            sys.exit(1)
        del args.vcfs[0]
        vcf1 = findUnion(vcf1, vcf2, args.combine)
    vcf1 = sortVcf(vcf1)
    try:
        out.write('\n'.join(vcf1))
    except:
        print '\n'.join(vcf1)


main()