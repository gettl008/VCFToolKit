#!/usr/bin/env python2

import argparse
import os
import sys



def filterByDepth(vcf, depth_stats, args):
    outlines = []
    outlines = outlines + [line.rstrip() for line in vcf if line.startswith('#')]
    statsdict = {}
    for line in depth_stats:
        if not line.startswith('#'):
            line = line.rstrip()
            stats_split = line.split('\t')
            statsdict[stats_split[0]] = [(float(stats_split[1]) - (float(stats_split[2]) * float(args.sd_lims))),
                                         (float(stats_split[1]) + (float(stats_split[2]) * float(args.sd_lims)))]
    for line in vcf:
        if line.startswith('#CHROM'):
            line = line.rstrip()
            samples = line.split('\t')
        elif not line.startswith('#'):
            line = line.rstrip()
            line = line.split('\t')
            format = line[8].split(':')
            gt = format.index('GT')
            ad = format.index('AD')
            dp = format.index('DP')
            return_line = 0
            for i, samp in enumerate(line[9:]):
                samp = samp.split(':')
                while len(samp) < len(format):
                    samp.append(".")
                if samp[ad] == ".":
                    samp[ad] = "0,0"
                if samp[dp] == ".":
                    samp[dp] = "0"
                sampstats = statsdict[samples[i + 9]]
                samp[gt] = samp[gt].replace('|', '/')
                split_gt = samp[gt].split('/')
                split_ad = samp[ad].split(',')
                if samp[gt] != "0/0" and samp[gt] != "./." and sampstats[0] < float(samp[dp]) < sampstats[1]:
                    # print "P1\t" + samp[gt] + "\t" + samp[dp]
                    if split_gt[0] == split_gt[1] and \
                                    (float(split_ad[int(split_gt[1])])/float(samp[dp])) > args.hom_altfreq:
                        # print samp
                        return_line = 1
                    elif split_gt[0] != split_gt[1] and \
                                    (float(split_ad[int(split_gt[1])])/float(samp[dp])) > args.het_altfreq:
                        # print samp
                        return_line = 1
                elif samp[gt] != "0/0" and (float(samp[dp]) < sampstats[0] or float(samp[dp]) > sampstats[1] or samp[gt] == "./."):
                    if samp[gt] == "./." or float(samp[dp]) == 0:
                        samp[gt] = "0/0"
                    elif split_gt[0] == split_gt[1] and \
                                    (float(split_ad[int(split_gt[1])])/float(samp[dp])) < args.hom_altfreq:
                        samp[gt] = "0/0"
                    elif split_gt[0] != split_gt[1] and \
                                    (float(split_ad[int(split_gt[1])])/float(samp[dp])) < args.het_altfreq:
                        samp[gt] = "0/0"
                line[i+9] = ":".join(samp)
            if return_line == 1:
                outlines.append('\t'.join(line))
            # else:
            #     print '\t'.join(line)
    return outlines




def main():
    parser = argparse.ArgumentParser(description='Returns CSV containing allele information and frequency')
    parser.add_argument('-v', '--vcf',
                        help='VCF file')
    parser.add_argument('-o', '--output',
                        help='Output vcf. File will be overwritten if it already exists.')
    parser.add_argument('-d', '--depth_stats',
                        help='Three column tab-delimited text file containing <SAMPLE> <MEAN> <STANDARD_DEVIATION>.'
                             'If present, header line should begin with "#"')
    parser.add_argument('-s', '--sd_lims', type=int,
                        help='Number of standard deviations from value to set min and max depth cutoffs. Default = 1.',
                        default=1)
    parser.add_argument('-f', '--het_altfreq', type=float,
                        help='Minimum frequency of alternative allele if heterozygous. Default = 0',
                        default=0.0)
    parser.add_argument('-F', '--hom_altfreq',  type=float,
                        help='Minimum frequency of alternative allele if homozygous. Default = 0',
                        default=0.0)
    args = parser.parse_args()
    try:
        args.output
        if os.path.exists(args.output):
            os.remove(args.output)
        try:
            out = open(args.output, "w+")
        except:
            print "Could not open file %s" % args.output
            sys.exit(1)
    except:
        pass

    try:
        vcf = open(args.vcf, "r").readlines()
    except:
        print "Could not open vcf %s" % args.vcf
        sys.exit(1)
    try:
        dpstats = open(args.depth_stats, "r").readlines()
    except:
        print "Could not open vcf %s" % args.depth_stats
        sys.exit(1)
    outlines = filterByDepth(vcf, dpstats, args)
    try:
        out.write('\n'.join(outlines))
    except:
        print '\n'.join(outlines)


main()