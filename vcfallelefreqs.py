#!/usr/bin/env python2

import argparse
import os
import sys
import numpy


def autocorr(x, t=1):
    a = numpy.array([x[0:len(x)-t], x[t:len(x)]]).astype(numpy.float)
    return numpy.corrcoef(a)[0][1]


def vcfSortBySample(vcf, samples):
    print "vcfsortbysample"
    outlines = [line.rstrip() for line in vcf if line.startswith('##')]
    header = [s for s in vcf if '#CHROM' in s][0].split('\t')
    header[-1] = header[-1].rstrip()
    samples = samples.split(',')
    sample_dict = {}
    for i, sample in enumerate(samples):
        sample_dict[i] = header.index(sample)
    outheader = header[0:9] + samples
    outlines.append('\t'.join(outheader))
    for line in vcf:
        if not line.startswith('#'):
            line = line.split('\t')
            line[-1] = line[-1].rstrip()
            outline = line[:9]
            for i in sample_dict.keys():
                outline.append(line[sample_dict[i]])
            outlines.append('\t'.join(outline))
    print "vcfsortbysample finish"
    return outlines


def alleleFreqs(vcf, annots, args):
    outdict = {}
    print "allelefreqs"
    outdict['header'] = '\n'.join([line.rstrip() for line in vcf if line.startswith('#')])
    for line in vcf:
        if not line.startswith('#'):
            # Break down VCF components into array
            line = line.split('\t')
            # Get the the variant position
            outline_pos = line[0:2]
            gt_info = line[8].split(':')
            dp_index = gt_info.index('DP')
            ad_index = gt_info.index('AD')
            alleles = [line[3]] + line[4].split(',')
            ref_allele = alleles[0]
            annot_info = line[7].replace('=',';').split(';')
            out_annots = []
            for annot in annots:
                if annot in annot_info:
                    annot_index = annot_info.index(annot) + 1
                    out_annots.append(annot_info[annot_index])
                else:
                    out_annots.append("NA")
            for i, allele in enumerate(alleles[1:]):
                outline = outline_pos[:]
                outline = outline + [ref_allele, allele] + out_annots
                dps = []
                for s in line[9:]:
                    s = s.split(':')
                    # Make sure that genotype fields match description
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
                                af = s[af_index].split(',')
                                if i < len(af):
                                    allele_freq = str(float(af[i]))
                                else:
                                    allele_freq = 0
                            else:
                                allele_freq = str(float(ad[i + 1]) / float(dp))
                        dps.append(float(dp))
                    outline.append(allele_freq)
                outline.append(str(sum(dps) / len(dps)))
                if '\t'.join(line) not in outdict:
                    outdict['\t'.join(line)] = ['\t'.join(outline)]
                else:
                    outdict['\t'.join(line)].append('\t'.join(outline))
    print "Finished allele freqs"
    return outdict


def filterVcfFreqs(vcfdict, vcf, args):
    print "filtervcffreqs"
    outdict = {'header': vcfdict['header']}
    outvcf = []
    # Q = 0
    for vcfline in vcf:
        # Q = Q + 1
        if not vcfline.startswith('#'):
            for f_line in vcfdict[vcfline]:
                freqline = f_line.split('\t')
                avg_dp = float(freqline.pop())
                annot_extra = len(args.output_annotations.split(','))
                freqs = [float(i) for i in freqline[4 + annot_extra:] if i != 'NA']
                if (args.min_depth and avg_dp >= float(args.min_depth)) or not args.min_depth:
                    fdiff = abs(max(freqs) - min(freqs))
                    avg_freq = sum(freqs) / len(freqs)
                    if (args.min_freq and float(args.min_freq) <= float(avg_freq) <= float(args.max_freq)) or not args.min_freq:
                            if args.autocorrelation and len(freqs) > 1:
                                # print(Q)
                                a = autocorr(freqs)
                            else:
                                a = 0
                            if (args.autocorrelation and abs(float(a)) >= float(args.autocorrelation)) or not args.autocorrelation:
                                outfreq = freqline[:4 + annot_extra] + [str(avg_dp), str(fdiff), str(a)] + freqline[4 + annot_extra:]
                                if vcfline not in outdict:
                                    outdict[vcfline] = ['\t'.join(outfreq)]
                                    outvcf.append(vcfline)
                                else:
                                    outdict[vcfline].append('\t'.join(outfreq))
        else:
            outvcf.append(vcfline)
    print "filtervcffreqs finish"
    return outdict, outvcf


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
                        help='Annotations to output (separated by commas) e.g MQ,ReadPosRankSum).', default="MQ")
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
    annots = args.output_annotations.split(',')
    outdict = alleleFreqs(vcf, annots, args)
    if args.diff or args.autocorrelation or args.min_depth or args.min_freq or args.max_freq:
        outdict, outvcf = filterVcfFreqs(outdict, vcf, args)
    #args.diff, args.autocorrelation, args.min_depth, args.min_freq,
    #                                 args.max_freq)
    samples = [s for s in outvcf if '#CHROM' in s][0].split('\t')[9:]
    outfreqs = ["CHROM\tPOS\tREF\tALT\t" + '\t'.join(annots) + "\tDP\tDIFF\tRk\t" + '\t'.join(samples).rstrip()]
    # print outfreqs[0]
    for line in outvcf:
        if not line.startswith('#'):
            outfreqs = outfreqs + outdict[line]
    try:
        outFREQpath.write('\n'.join(outfreqs))
    except:
        print '\n'.join(outfreqs)
    try:
        outVCFpath.write(''.join(outvcf))
    except:
        print ''.join(outvcf)



main()



# Issues:
# Chromosomes need to be in correct order in  header
# BK006938.2