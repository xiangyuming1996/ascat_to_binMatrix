# ascat to bin matrix python script
# for CNV heatmap

import sys
import pandas as pd

def make_bins(filename):
    """
    return all genomic bins
    as list of [[chr, start, end],...]
    """
    res = []
    for line in open(filename):
        chromosome, length = line.strip('\n').split('\t')[:2]
        chromosome, length = map(int, [chromosome, length])

        # only autosome, chr1-22
        if chromosome >= 23:
            continue

        start = 1
        while start + 1000000 < length:
            end = start + 1000000 -1
            res.append([chromosome, start, end])
            start += 1000000
        res.append([chromosome, start, length])

    print('make bins succed, in total: %s\n' % len(res))
    return res

def overlap_section(a, b):
    """
    return the overlap section of a and b
    return False if a, b is not overlapped
    """
    left, right = max(a[0], b[0]), min(a[1], b[1])
    if left < right:
        return [left, right]
    else:
        return False

def in_section(a,b):
    """
    return True if section a is in section b
    """
    return a[0] > b[0] and a[1] < b[1]

def readascat(filename):
    """read ASCAT format file
    return a dictionary with sample name as key,
    list of [chromosome, start, end, copynumber] as value
    24 for chromosome x
    25 for y
    """
    d = {}
    segcount = 0
    for line in open(filename):
        l = line.strip('\n').split('\t')
        [sample, chromosome, start, end, nMajor, nMinor] = map(lambda x: int(x) if x.isdigit() else x, l)

        # skip table head
        if sample == 'sample':
            continue

        # only autosome
        if not chromosome <= 22:
            continue

        # nMinor should always be 0
        cn = sum([nMajor, nMinor])

        if sample not in d:
            d[sample] = []
        #print [sample, chromosome, start, end, cn]
        d[sample].append([chromosome, start, end, cn])
        segcount += 1

    print('read ascat file succed\nsamples count: %s\nsegment count: %s\n' % (len(d), segcount))

    return d

def align_ascat_to_bins(ascat, bins):
    """
    heavy lift jobs code
    read each sample's segment data, align and call real copy number for
    each bins
    """
    res = {}

    print('start align segments in ascat to each window(bin)')
    print('%s bins in total' % len(bins))
    print('%s samples in total' % len(ascat))

    for sample, segments in ascat.items():
        print(sample)
        res[sample] = []
        for win in bins:
            contain = []
            win_chr, win_start, win_end = win
            for seg in segments:
                seg_chr, seg_start, seg_end, seg_cn = seg
                if win_chr != seg_chr:
                    continue
                overlap = overlap_section([win_start, win_end], [seg_start, seg_end])
                if overlap:
                    ratio = (overlap[1] - overlap[0])/(1000000 - 1.0)
                    relcn += ratio*seg_cn
                    contain.append([ratio, seg_cn])
                    #print win, seg, overlap
            if not contain:
                relcn = 2
            else:
                contain.append([1 - sum([i[0] for i in contain]), 2])
                primary_seg = sorted(contain, key=lambda x:x[0], reverse=True)
                relcn = primary_seg[0][1]
            res[sample].append(relcn)
            #print sample, win, contain, relcn
            #print sample, win, contain, relcn
    return res

def machine(filename):
    ascat = readascat(filename)
    bins = make_bins('chr')

    res = align_ascat_to_bins(ascat, bins)
    loc_index = pd.Index(['.'.join([str(i) for i in b]) for b in bins])
    df = pd.DataFrame(res, index=loc_index).T

    col = pd.Index(df.columns, name='sample')
    df.columns = col
    outfilename = filename + '.binMat.txt'

    df.to_csv(outfilename, sep='\t', header=True, index=True)
    print df
    print 'file saved as', outfilename
    print 'done'
    return 0

def usage():
    print '\nUsage:'
    print '\tpython ascat_to_binMatrix.py file_of_ascat_cnv_segment.txt'
    print '\noutput will write to file'
    return 0

def main():
    if len(sys.argv) == 1:
        usage()
    elif sys.argv[-1] == '-h':
        usage()
    else:
        f = sys.argv[1]
        machine(f)

if __name__ == "__main__":
    main()
