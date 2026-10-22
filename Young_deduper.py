import argparse
import pathlib

EPILOG="""\033[91mInput SAM file must be sorted by chromosome!\033[0m
        e.g., samtools sort -o <in.sorted.sam> <in.sam>"""

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
        epilog=EPILOG
    )
    parser.add_argument("-f", "--file", type=pathlib.Path, required=True,
                        help="designates absolute file path to sorted sam file")
    parser.add_argument("-o", "--outfile", type=pathlib.Path, required=True,
                        help="designates absolute file path to deduplicated sam file")
    parser.add_argument("-u", "--umi", type=pathlib.Path, required=False, default="STL96.txt",
                        help="designates file containing the list of UMIs \033[4m(default: %(default)s)\033[0m")

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
    if not reverse:
        if cigar[1] == 'S':
            pos -= int(cigar[0])
        return pos
    
    for i in range(0, len(cigar), 2):
        length, op = int(cigar[i]), cigar[i + 1]
        if i == 0 and op == 'S':
            continue
        if op in ('M', 'D', 'N', 'S'):
            pos += length
    return pos - 1  

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
        list: Number of header lines, unique reads, incorrect
        UMIs, duplicates removed, and chromosome information.
    """
    umi_set = get_all_umis(umi_file)
    curr_chrom = "0"
    curr_set = set()

    header = 0
    unique = 0
    bad_umi = 0
    dup_rem = 0
    chrom_dict = {}

    with open(infile, 'r') as in_sam, open(outfile, 'w') as out_sam:
        for line in in_sam:
            stripped = line.strip('\n')

            if stripped[0] == '@':
                out_sam.write(line)
                header += 1
                continue

            sam_line = stripped.split('\t')

            umi = get_umi(sam_line)
            if umi not in umi_set:
                bad_umi += 1
                continue

            chrom = sam_line[2]
            if chrom != curr_chrom:
                curr_chrom = chrom
                curr_set.clear()
                chrom_dict[chrom] = 0

            pos = int(sam_line[3])
            reverse = int(sam_line[1]) & 16 == 16
            cigar = parse_cigar(sam_line[5])

            read_start = get_read_start(reverse, pos, cigar)

            seq = (umi, reverse, read_start)
            if seq not in curr_set:
                curr_set.add(seq)
                out_sam.write(line)
                unique += 1
                chrom_dict[chrom] += 1
            else:
                dup_rem += 1

        return [header, unique, bad_umi, dup_rem, chrom_dict]

def print_results(stats):
    """
    Prints statistics on PCR duplicate removal.

    Args:
        stats (list): Number of header lines, unique reads, incorrect
        UMIs, duplicates removed, and chromosome information.
    
    Returns:
        None.
    """
    header, unique, bad_umi, dup_rem, chrom_dict = stats

    print(f"Number of header lines: {header}")
    print(f"Number of unique reads: {unique}")
    print(f"Number of wrong UMIs: {bad_umi}")
    print(f"Number of duplicates removed: {dup_rem}")

    print("\nNumber of reads per chromosome")
    for k, v in chrom_dict.items():
        print(f"{k}\t{v}")

    return

def main():
    args = get_args()

    infile = args.file
    outfile = args.outfile
    umi_file = args.umi

    stats = dedup(infile, outfile, umi_file)
    print_results(stats)

if __name__ == "__main__":
    main()