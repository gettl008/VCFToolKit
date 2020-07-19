#!/usr/bin/env python2

import argparse
import os
import sys
import numpy


def alleleFreqs(vcf, args):
    outfile = []
    print "allelefreqs"
    samples = [s for s in vcf if '#CHROM' in s][0].split('\t')[9:]
    x = 0
    for line in vcf:
        if not line.startswith('#'):
            # Break down VCF components into array
            print line
            line = line.split('\t')
            # Get the the variant position
            outline_pos = line[0:2]
            gt_info = line[8].split(':')
            dp_index = gt_info.index('DP')
            ad_index = gt_info.index('AD')
            alleles = [line[3]] + line[4].split(',')
            ref_allele = alleles[0]
            for i, s in enumerate(line[9:]):
                sample = samples[i - 1]
                outline = [sample] + outline_pos[:] + [ref_allele, alleles[1]]
                s = s.split(':')
                dps = []
                if len(s) < len(gt_info):
                    allele_freq = 'NA'
                    dps.append(0.0)
                else:
                    dp = s[dp_index]
                    ad = s[ad_index].split(',')
                    while len(ad) < len(alleles):
                        ad.append('0')
                    if dp == '.' or int(dp) == 0:
                        allele_freq = str(0)
                        dp = 0
                    else:
                        if args.dp_out:
                            allele_freq = str(dp)
                        elif args.AF:
                            af_index = gt_info.index('AF')
                            af = s[af_index].split(',')[0]
                            print af
                            if af != '.':
                                allele_freq = str(float(af))
                            else:
                                allele_freq = 0
                            print allele_freq
                        else:
                            allele_freq = str(float(ad[i + 1]) / float(dp))
                dps.append(float(dp))
                outline.append(allele_freq)
                outline = '\t'.join(outline)
                outfile.append(outline)
        x += 1
    print x
    return outfile


def main():
    parser = argparse.ArgumentParser(description='Returns CSV containing allele information and frequency')
    parser.add_argument('vcf', metavar='vcf', type=str, nargs=1,
                        help='VCF file')
    parser.add_argument('-o', '--output_base',
                        help='Output basename. Files will be overwritten if it already exists.')
    parser.add_argument('-s', '--samples',
                        help='Quoted string of samples in the order in which to sort.')
    parser.add_argument('-d', '--diff',
                        help='Minimum difference between max and min frequencies to pass filter.')
    parser.add_argument('-a', '--autocorrelation',
                        help='Minimum autocorrelation with lag=1.')
    parser.add_argument('-m', '--min_depth',
                        help='Minimum average depth for a position.')
    parser.add_argument('-f', '--min_freq',
                        help='Minimum average alternate allele frequency for a position.')
    parser.add_argument('-F', '--max_freq',
                        help='Maximum average alternate allele frequency for a position.', default = 1)
    parser.add_argument('-A', '--output_annotations',
                        help='Annotations to output (separated by commas) e.g MQ,ReadPosRankSum).', default = "")
    parser.add_argument('-D', '--dp_out',
                        help='Only output site read depth instead of freq.', action='store_true')
    parser.add_argument('-l', '--AF',
                        help='Use allele frequency information provided by the caller', action='store_true')
    args = parser.parse_args()

    try:
        args.output_base
        outvcfpath = args.output_base + ".vcf"
        outfreqpath = args.output_base + ".txt"
        if os.path.exists(outvcfpath):
            os.remove(outvcfpath)
        if os.path.exists(outfreqpath):
            os.remove(outfreqpath)
        try:
            outVCFpath = open(outvcfpath, "w+")
        except:
            print "Could not open file %s" % outvcfpath
            sys.exit(1)
        try:
            outFREQpath = open(outfreqpath, "w+")
        except:
            print "Could not open file %s" % outfreqpath
            sys.exit(1)
    except:
        pass

    try:
        vcf = open(args.vcf[0], "r").readlines()
    except:
        print "Could not open vcf %s" % args.vcfs[0]
        sys.exit(1)
    if args.samples:
        vcf = vcfSortBySample(vcf, args.samples)
    outlines = alleleFreqs(vcf, args)
    # print outfreqs[0]
    try:
        outFREQpath.write('\n'.join(outlines))
    except:
        print '\n'.join(outlines)
    try:
        outVCFpath.write(''.join(outlines))
    except:
        print ''.join(outlines)



main()