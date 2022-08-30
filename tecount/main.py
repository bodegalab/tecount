#!/usr/bin/env python

import argparse, gzip, os, sys, subprocess
from shutil import rmtree
from datetime import datetime
from functools import partial
import multiprocessing
from copy import deepcopy
from pysam import AlignmentFile, index

__version__ = '0.6.0'

def parseArguments():
    parser = argparse.ArgumentParser(prog='TEcount',
                                     usage='TEcount -b <file.bam> -r <rmsk.bed.gz> [OPTIONS]',
                                     description='Count reads mapping on Transposable Elements subfamilies, families and classes.',
                                     epilog='Documentation and issue tracker: https://github.com/bodegalab/tecount')
    parser.add_argument('-b','--bam', required=True, help='scRNA-seq reads aligned to a reference genome')
    parser.add_argument('-r', '--rmsk', required=True, help='Genomic TE coordinates in bed format, with subfamily, family and class on columns 7, 8 and 9. Plain text or gzip-compressed.')
    parser.add_argument('-o', '--overlap', type=int, default=1, help='Minimum bp overlap between read and feature (default: 1).')
    parser.add_argument('-s', '--strandness', type=str, default='unstranded', help='Strandness of the library. One of: "unstranded" (default), "forward", "reverse".')
    parser.add_argument('--prefix', type=str, default='', help='Prefix for output file names (default: no prefix).')
    parser.add_argument('--outdir', type=str, default='./', help='Output directory (default: current directory).')
    parser.add_argument('--tmpdir', type=str, default='./tmp', help='Temporary files directory (default: ./tmp).')
    parser.add_argument('--keeptmp', default=False, action='store_true', help='Keep temporary files (default: False).')
    parser.add_argument('-p','--threads', type=int, default=1, help='Number of cpus to use (default: 1)')
    parser.add_argument('-V','--version', action='version', version='%(prog)s {}'.format(__version__), help='Print software version and quit.')
    return parser

# Wrapper function to execute subprocesses
def run_shell_cmd(cmd):
    p = subprocess.Popen(
        ['/bin/bash', '-o', 'pipefail'],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        preexec_fn=os.setsid,
        text=True
    )
    pid = p.pid
    pgid = os.getpgid(pid)
    stdout, stdin = p.communicate(cmd)
    return stdout.strip('\n')

# Small function to write out the current date and time
def time():
    return datetime.now().strftime("%m/%d/%Y - %H:%M:%S")

# Small function to write a message to stderr with timestamp
def writerr(msg, send=True):
    if send:
        message = f'[{time()}] '
        if not msg[-1]=='\n':
            msg += '\n'
        message += msg
        sys.stderr.write(message)

# Check if bam file is indexed
def checkIndex(bamFile):
    with AlignmentFile(bamFile) as bam:
        if not bam.has_index():
            writerr('BAM index not found. Attempting to index the BAM...')
            try:
                index(bamFile)
            except:
                sys.exit(f'Couldn\'t index the BAM file. Please do so manually with `samtools index {bamFile}`.')
            else:
                writerr('BAM indexing done.')
        else:
            writerr(f'Found index for BAM file {bamFile}')

# Function to retrieve reference names from bam file
def getRefs(bamFile, tmpdir):
    os.makedirs(tmpdir, exist_ok=True)
    chrFile = tmpdir + '/' + 'chrFile.txt'
    chrNames = list()
    cmd = f'samtools idxstat {bamFile} > {chrFile}'
    run_shell_cmd(cmd)
    with open(chrFile, 'r') as f:
        for line in f:
            l = line.strip().split('\t')
            if int(l[2])>0 and not l[0]=='*':
                chrNames.append(l[0])
    return chrNames

# Check if file is gzip
def testGz(input_file):
    '''Test if file is gzip'''
    with gzip.open(input_file, 'rb') as f:
        try:
            f.read(1)
            return True
        except gzip.BadGzipFile:
            return False

# returns 3 dictionaries of empty sets of features in a tuple
def getFeaturesDicts(bedFile, tmpdir):
    os.makedirs(tmpdir, exist_ok=True)
    if testGz(bedFile):
        openfile = lambda x: gzip.open(x, 'rb')
        rl = lambda x: x.decode('utf-8').strip().split('\t')
    else:
        openfile = lambda x: open(x, 'r')
        rl = lambda x: x.strip().split('\t')
    with openfile(bedFile) as f:
        subfams = set()
        fams = set()
        classes = set()
        for line in f:
            l = rl(line)
            sfam, fam, clas = l[6:9]
            subfams.add(sfam)
            fams.add(fam)
            classes.add(clas)
    subfams = { x: set() for x in subfams }
    fams = { x: set() for x in fams }
    classes = { x: set() for x in classes }
    return subfams, fams, classes

# intersect reads with features using bedtools
def isec(bamFile, bedFile, overlap, tmpdir, chr):
    os.makedirs(tmpdir, exist_ok=True)
    refFile = tmpdir + '/' + chr + '.bed.gz'
    isecFile = tmpdir + '/' + chr + '.isec.gz'

    # extract chromosome
    if testGz(bedFile):
        cmd0 = f'zcat {bedFile} | awk \'$1=="{chr}"\' | gzip > {refFile}'
    else:
        cmd0 = f'awk \'$1=="{chr}"\' {bedFile} > {refFile}'

    # intersect
    cmd = f'samtools view -u {bamFile} {chr} | '
    cmd += ' bedtools bamtobed -i stdin -splitD | LC_ALL=C sort -k1,1 -k2,2n | '
    #cmd += ' bedtools bamtobed -i stdin -bed12 -split -splitD | '
    cmd += f' bedtools intersect -a stdin -b {refFile} -wo -sorted | '
    cmd += f' awk \' $NF>={overlap} '
    cmd += ' { l=split($4,readid,"/"); readname=readid[1]; '
    # if single-end, set mate=1 even if no mate is recorded in read name
    cmd += ' if(length(l)==1){mate=1}else{mate=readid[2]} '
    cmd += ' ori=($6==$12 ? "same" : "opposite"); '
    #cmd += ' ori=($6==$18 ? "same" : "opposite"); '
    cmd += ' OFS="\\t"; print readname,mate,$13,$14,$15,ori; } \' | '
    cmd += f' gzip > {isecFile} '

    writerr(f'Extracting {chr} reference')
    run_shell_cmd(cmd0)
    writerr(f'Intersecting alignments with {chr} reference')
    run_shell_cmd(cmd)

    return isecFile

# count reads mapping on features
def count(filesList, featuresTuple, strandness, prefix, outdir):
    os.makedirs(outdir, exist_ok=True)
    # create one dictionary for each output file
    sfcounts, fmcounts, clcounts = deepcopy(featuresTuple)
    sfcountsS, fmcountsS, clcountsS = deepcopy(featuresTuple)
    sfcountsAS, fmcountsAS, clcountsAS = deepcopy(featuresTuple)
    for file in filesList:
        with gzip.open(file, 'rb') as f:
            for line in f:
                l = line.decode('utf-8').strip().split('\t')
                readid, mate, sfam, fam, clas, ori = l
                mate = int(mate)
                # count unstranded
                sfcounts[sfam].add(readid)
                fmcounts[fam].add(readid)
                clcounts[clas].add(readid)
                if not strandness == 'unstranded':
                    if (
                        strandness == 'forward' and (
                            (
                                mate == 1 and ori == 'same'
                            ) or (
                                mate == 2 and ori == 'opposite'
                            )
                        )
                    ) or (
                        strandness == 'reverse' and (
                            (
                                mate == 1 and ori == 'opposite'
                            ) or (
                                mate == 2 and ori == 'same'
                            )
                        )
                    ):
                        # count sense
                        sfcountsS[sfam].add(readid)
                        fmcountsS[fam].add(readid)
                        clcountsS[clas].add(readid)
                    else:
                        # count antisense
                        sfcountsAS[sfam].add(readid)
                        fmcountsAS[fam].add(readid)
                        clcountsAS[clas].add(readid)
    
    
    # write output files
    if len(prefix) > 0 and prefix[-1] != '_':
        prefix += '_'

    sfunsfile = f'{prefix}subfamily_unstranded.count.txt'
    sfunsfile = os.path.join(outdir, sfunsfile)
    fmunsfile = f'{prefix}family_unstranded.count.txt'
    fmunsfile = os.path.join(outdir, fmunsfile)
    clunsfile = f'{prefix}class_unstranded.count.txt'
    clunsfile = os.path.join(outdir, clunsfile)

    with open(sfunsfile, 'w') as f:
        f.writelines([ f'{x}\t{str(len(y))}\n' for x, y in sfcounts.items() ])
    with open(fmunsfile, 'w') as f:
        f.writelines([ f'{x}\t{str(len(y))}\n' for x, y in fmcounts.items() ])
    with open(clunsfile, 'w') as f:
        f.writelines([ f'{x}\t{str(len(y))}\n' for x, y in clcounts.items() ])

    if not strandness == 'unstranded':
        sffwdfile = f'{prefix}subfamily_sense.count.txt'
        sffwdfile = os.path.join(outdir, sffwdfile)
        fmfwdfile = f'{prefix}family_sense.count.txt'
        fmfwdfile = os.path.join(outdir, fmfwdfile)
        clfwdfile = f'{prefix}class_sense.count.txt'
        clfwdfile = os.path.join(outdir, clfwdfile)

        sfrevfile = f'{prefix}subfamily_antisense.count.txt'
        sfrevfile = os.path.join(outdir, sfrevfile)
        fmrevfile = f'{prefix}family_antisense.count.txt'
        fmrevfile = os.path.join(outdir, fmrevfile)
        clrevfile = f'{prefix}class_antisense.count.txt'
        clrevfile = os.path.join(outdir, clrevfile)

        with open(sffwdfile, 'w') as f:
            f.writelines([ f'{x}\t{str(len(y))}\n' for x, y in sfcountsS.items() ])
        with open(fmfwdfile, 'w') as f:
            f.writelines([ f'{x}\t{str(len(y))}\n' for x, y in fmcountsS.items() ])
        with open(clfwdfile, 'w') as f:
            f.writelines([ f'{x}\t{str(len(y))}\n' for x, y in clcountsS.items() ])

        with open(sfrevfile, 'w') as f:
            f.writelines([ f'{x}\t{str(len(y))}\n' for x, y in sfcountsAS.items() ])
        with open(fmrevfile, 'w') as f:
            f.writelines([ f'{x}\t{str(len(y))}\n' for x, y in fmcountsAS.items() ])
        with open(clrevfile, 'w') as f:
            f.writelines([ f'{x}\t{str(len(y))}\n' for x, y in clcountsAS.items() ])


def main():

    parser = parseArguments()
    args = parser.parse_args()

    checkIndex(args.bam)

    chrNames = getRefs(args.bam, args.tmpdir)

    pool = multiprocessing.Pool(args.threads)
    func = partial(isec, args.bam, args.rmsk, args.overlap, args.tmpdir)
    isecFiles = pool.map(func, chrNames)

    featuresDicts = getFeaturesDicts(args.rmsk, args.tmpdir)

    count(isecFiles, featuresDicts, args.strandness, args.prefix, args.outdir)

    if not args.keeptmp:
        writerr('Cleaning up temporary files...')
        rmtree(args.tmpdir)

    writerr('Done.')

if __name__ == "__main__":
    main()
