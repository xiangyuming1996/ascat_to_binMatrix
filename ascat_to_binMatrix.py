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
    return max(a[0], b[0]) < min(a[1], b[1])

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
        d[sample].append([chromosome, start, end,cn])

    return d

def machine(filename):
    data = readascat(filename)
    for key, value in data.items():
        print key, value

def main():
    f = sys.argv[1]
    machine(f)

if __name__ == "__main__":
    main()
