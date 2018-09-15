import sys

def make_bin(filename):
    res = []
    for line in open(filename):
        start = 0
        while


def machine(filename):
    """read ASCAT format file
    """
    d = {}
    for line in open(filename):
        l = line.strip('\n').split('\t')
        if len(l) != 6:
            raise Exception('column number not 6')
        [sample, chromosome, start, end, nMajor, nMnior] = l

def main():
    f = sys.argv[1]
    machine(f)

if __name__ == "__main__":
    main()
