#!/usr/bin/env python2

import argparse
import os
import sys
import re

###############################################

def lc_dictionary(lc):
    lc_dict = {}
    for line in lc:
        if not line.startswith('#'):
            line = line.split('\t')
            if line[0] not in lc_dict.keys():
                lc_dict[line[0]] = {}
            lc_dict[line[0]][line[3]] = line[4]
    return(lc_dict)

def remove_lc(vcf, lc_dict):
    filtered_lines = [line.rstrip() for line in vcf if line.startswith('#')]
    lc_lines = [line.rstrip() for line in vcf if line.startswith('#')]
    for line in vcf:
        if not line.startswith('#'):
            splitline = line.split('\t')
            chr = splitline[0]
            pos = int(splitline[1])
            lc_keys = map(int, lc_dict[chr].keys())
            lc_keys.sort()
            # lc_keys = int(lc_dict[chr].keys())
            # print lc_keys
            i = 0
            print_line = True
            while(i < len(lc_dict[chr].keys())):
                s = lc_keys[i]
                e = int(lc_dict[chr][str(s)])
                if s <= pos <= e:
                    # print str(s) + " " + str(pos) + " " + str(e)
                    lc_lines.append(line.rstrip())
                    print_line = False
                    break
                i += 1
            if print_line:
                filtered_lines.append(line.rstrip())
    return(filtered_lines, lc_lines)





def main():
    parser = argparse.ArgumentParser(description='Generates a VCF of overlapping variants from a set of VCF files.')
    parser.add_argument('vcfs', metavar='vcf', type=str, nargs=1,
                        help='VCF file')
    parser.add_argument('-o', '--output_base',
                        help='Output VCF files basename. Files will be overwritten if they already exists.',
                        default=os.getcwd() + "/out")
    parser.add_argument('-R', '--repeat_regions',
                        help='Repeat/low complexity regions in GFF format. Output from RepeatMasker. Required.')
    args = parser.parse_args()
    try:
        vcf = open(args.vcfs[0], "r").readlines()
    except:
        print "Could not open vcf %s" % args.vcfs[0]
        sys.exit(1)
    try:
        repeats = open(args.repeat_regions, "r").readlines()
    except:
        print "Could not open repeat region gff %s" % args.repeat_regions
        sys.exit(1)
    try:
        args.output_base
        filtered_out = args.output_base + ".filtered.vcf"
        lc_out = args.output_base + ".lc.vcf"
        if os.path.exists(filtered_out):
            os.remove(filtered_out)
        try:
            out = open(filtered_out, "w+")
        except:
            print "Could not open file %s" % filtered_out
            sys.exit(1)
        if os.path.exists(lc_out):
            os.remove(lc_out)
        try:
            lcout = open(lc_out, "w+")
        except:
            print "Could not open file %s" % lc_out
            sys.exit(1)
    except:
        pass

    lc = lc_dictionary(repeats)
    filtered_vcf, lc_vcf = remove_lc(vcf, lc)
    try:
        out.write('\n'.join(filtered_vcf))
    except:
        print '\n'.join(filtered_vcf)
    try:
        lcout.write('\n'.join(lc_vcf))
    except:
        print "No LC output"

main()