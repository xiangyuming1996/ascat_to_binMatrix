# ascat to bin matrix python script

import sys

def make_bin(filename):
    """
    return all genomic bins
    """
    res = []
    for line in open(filename):
        chromosome, length = line.strip('\n').split('\t')[:2]
        chromosome, length = map(int, [chromosome, length])

        start = 1
        while start + 1000000 < length:
            end = start + 1000000 -1
            res.append([chromosome, start, end])
            start += 1000000
        res.append([chromosome, start, length])

    return res

def overlap_section(a, b):
    """
    return the overlap section of a and b
    """
    left, right = max(a[0], b[0]), min(a[1], b[1])
    if left < right:
        return [left, right]
    else:
        return False

def in_section(a,b):
    """
    return True if section a in section b
    """
    return a[0] > b[0] and a[1] < b[1]

def readascat(filename):
    """read ASCAT format file
    return a dictionary with sample name as key
    24 for chromosome x
    25 for y
    """
    d = {}

    for line in open(filename):
        l = line.strip('\n').split('\t')
        [sample, chromosome, start, end, nMajor, nMinor] = map(lambda x: int(x) if x.isdigit() else x, l)
        if sample == 'sample':
            continue
        cn = sum([nMajor, nMinor])
        if sample not in d:
            d[sample] = []
        #print [sample, chromosome, start, end, cn]
        d[sample].append([chromosome, start, end, cn])

    return d

def machine(filename):
    ascat = readascat(filename)
    bins = make_bin('chr')
    for sample, segments in ascat.items():
        #print sample
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

            #print sample, win, contain, relcn
            print sample, win, contain, relcn


def main():
    f = sys.argv[1]
    machine(f)

if __name__ == "__main__":
    main()
