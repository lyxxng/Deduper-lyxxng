import argparse
import pathlib

def get_args():
    """
    Create argparse object for script.

    Args:
        None
    
    Returns:
        ArgumentParser: argparse object specifying -f, -o, and -u arguments.
    """
    parser = argparse.ArgumentParser(
        prog="Deduper",
        description="Reference Based PCR Duplicate Removal tool",
    )
    parser.add_argument("-f", "--file", type=pathlib.Path, required=True,
                        help="designates absolute file path to sorted sam file")
    parser.add_argument("-o", "--outfile", type=pathlib.Path, required=True,
                        help="designates absolute file path to deduplicated sam file")
    parser.add_argument("-u", "--umi", type=pathlib.Path, required=False, default="STL96.txt",
                        help="designates file containing the list of UMIs (default: %(default)s)")

    return parser.parse_args()

def get_all_umis(umi_file):
    """
    Get set of UMIs.

    Args:
        umi_file (string): File path of list of UMIs.
    
    Returns:
        set: All UMIs from file.
    """
    umi_set = set()

    with open(umi_file, 'r') as f:
        for line in f:
            umi_set.add(line.strip('\n'))

    return umi_set

def parse_cigar(cigar_str):
    """
    Parse CIGAR string.

    Args:
        cigar_str (string): CIGAR string from SAM file.
    
    Returns:
        list: Parsed CIGAR separating numbers and letters.
    """
    temp = ""
    cigar_lst = []

    for char in cigar_str:
        if char.isalpha():
            cigar_lst.extend([temp, char])
            temp = ""
        else:
            temp = temp + char

    return cigar_lst

def get_read_start(reverse, pos, cigar):
    """
    Get 5' start of read.

    Args:
        reverse (bool): Whether sequence is the reverse strand.
        pos (int): 1-based left-most mapping position.
        cigar (list): Parsed CIGAR string.
    
    Returns:
        int: 5' start of read given strand, position, and CIGAR.
    """
    if not reverse and cigar[1] == 'S':
        pos -= int(cigar[0])
    elif reverse:
        for i in range(0, len(cigar), 2):
            pos += int(cigar[i])
    return pos

def get_umi(sam_line):
    """
    Get UMI from line in SAM file.

    Args:
        sam_line (string): Single line from SAM file.
    
    Returns:
        string: UMI.
    """
    qname = sam_line[0]

    return qname.split(':')[-1]

def dedup(infile, outfile, umi_file):
    """
    Removes all PCR duplicates and writes deduplicated reads to new SAM file.

    Args:
        infile (string): File path of sorted SAM file.
        outfile (string): File path of deduplicated SAM file.
        umi_file (string): File path of list of UMIs.
    
    Returns:
        None
    """
    umi_set = get_all_umis(umi_file)
    curr_chrom = 0
    curr_set = set()

    with open(infile, 'r') as in_sam, open(outfile, 'w') as out_sam:
        for line in in_sam:
            stripped = line.strip('\n')
            # TODO: Check if I should write these lines to out file
            if stripped[0] == '@':
                continue
            sam_line = stripped.split('\t')

            umi = get_umi(sam_line)
            if umi not in umi_set:
                continue

            chrom = int(sam_line[2])
            if chrom != curr_chrom:
                curr_chrom = chrom
                curr_set.clear()

            pos = int(sam_line[3])
            reverse = int(sam_line[1]) & 16 == 16
            cigar = parse_cigar(sam_line[5])

            read_start = get_read_start(reverse, pos, cigar)

            seq = (umi, chrom, reverse, read_start)
            if seq not in curr_set:
                curr_set.add(seq)
                out_sam.write(line)
        return

def main():
    args = get_args()

    infile = args.file
    outfile = args.outfile
    umi_file = args.umi

    dedup(infile, outfile, umi_file)

if __name__ == "__main__":
    main()