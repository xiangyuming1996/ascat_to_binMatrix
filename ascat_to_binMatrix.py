
# ascat to bin matrix python script
# for CNV heatmap

import sys
import pandas as pd
from collections import Counter

def make_bins(filename='hg19.chrom.size'):
    """
    return all genomic bins
    as list of [[chr, start, end],...]
    """

    print('make bins from file: %s' % filename)

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

    print('make bins successfully, in total: %s\n' % len(res))
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
    """
    d = {}
    segcount = 0
    print('\nreading ascat file: %s' % filename)

    for line in open(filename):
        l = line.strip('\n').split('\t')
        l = list(map(lambda x: int(x) if x.isdigit() else x, l))
        [sample, chromosome, start, end, nMajor, nMinor] = l

        # skip table head
        if sample == 'sample':
            continue

        # only autosome
        if isinstance(chromosome, str):
            continue
        if chromosome >= 23:
            continue

        # nMinor should always be 0
        cn = sum([nMajor, nMinor])

        if sample not in d:
            d[sample] = []
        #print [sample, chromosome, start, end, cn]
        d[sample].append([chromosome, start, end, cn])
        segcount += 1

    print()
    print('read ascat file successfully!')
    print('samples count: %s, %s' % (len(d), ' '.join(d.keys())))
    print('segment count: %s' % segcount)
    print()

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

            # for each window, record segment ratio and copy number
            contain = []

            win_chr, win_start, win_end = win
            for seg in segments:
                seg_chr, seg_start, seg_end, seg_cn = seg

                if win_chr != seg_chr:
                    continue

                overlap = overlap_section([win_start, win_end], [seg_start, seg_end])
                if overlap:
                    ratio = (overlap[1] - overlap[0])/(1000000 - 1.0)

                    # I don't know what's this, it's really strange code.
                    # maybe mis-typing ?
                    #relcn += ratio*seg_cn

                    contain.append([ratio, seg_cn])
                    #print win, seg, overlap
            if not contain:
                relcn = 2
            else:
                # all other space will be set to copy number of 2, by default.
                contain.append([1 - sum([i[0] for i in contain]), 2])

                # find longest segment, use corresponding cn for this window.
                primary_seg = sorted(contain, key=lambda x:x[0], reverse=True)
                relcn = primary_seg[0][1]

            res[sample].append(relcn)
            #print sample, win, contain, relcn
    return res

def machine(filename):

    ascat = readascat(filename)
    bins = make_bins()

    res = align_ascat_to_bins(ascat, bins)
    loc_index = pd.Index(['.'.join([str(i) for i in b]) for b in bins])
    df = pd.DataFrame(res, index=loc_index).T

    #col = pd.Index(df.columns, name='sample')
    df.index.name = 'sample'
    outfilename = filename + '.binMat.txt'

    df.to_csv(outfilename, sep='\t', header=True, index=True)
    print(df)
    print('file saved as', outfilename)
    print()

    make_chr_colorband(bins)
    make_chr_nameprobe(bins)

    return 0

def make_chr_colorband(bins):
    f = 'chromosome.colorband.txt'
    with open(f, 'w') as handle:
        for i in bins:
            chrn = i[0]
            if chrn % 2 == 0:
                handle.write('even\t')
            else:
                handle.write('odd\t')
        handle.write('\n')
    print('make chromosome colorband file successfully!, file saved to: "%s"' % f)
    print()
    return 0

def make_chr_nameprobe(bins):

    res = []
    chrs = [i[0] for i in bins]
    chr_count = Counter(chrs).most_common(22)
    chr_count = sorted(chr_count, key=lambda x: x[0])
    res.append(chr_count[0][1] / 2)

    for i in range(1,22):
        res.append(sum([cc[1] for cc in chr_count[:i]]) + chr_count[i][1] / 2 )

    res = list(map(round, res))
    f = 'chromosome.nameprobe.txt'
    open(f, 'w').write('\t'.join(map(str, res)) + '\n')
    print('make chr name probe index successfully, file saved to: "%s"' % f)
    print()
    return 0

def usage():
    print('\nUsage:')
    print('\tpython ascat_to_binMatrix.py file_of_ascat_cnv_segment.txt')
    print('\noutput will write to file')

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
    print('done')
