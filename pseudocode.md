### Define the problem
Strands are duplicated during PCR. These PCR duplicates can give a false sense of sequencing depth, since not every piece of DNA will amplify equally. Essentially, it will affect anything where we are quantifying how often a sequence appears in a sample. To remove PCR duplicates, we need to take into account the chromosome, 5' start of the read, strand, and UMI.

### Pseudocode
```
def get_args():
    ''' Function that creates and returns an argparse object that
    defines the file, outfile, umi, and help options. '''
    # Input: None
    # Output: argparse object

    define argparse object
    add file, outfile, umi, help

    return argparse object

def get_umis(umi_file):
    ''' Function that takes UMI file path as input and returns
    a set containing all 96 unique UMIs'''
    # Input: UMI file
    # Output: Set with all UMIs from file

    initialize empty umi set

    open umi file to read:
        read line, remove '\n'
        add line to umi set
    
    return set

def parse_cigar(cigar_str):
    ''' Function that takes the CIGAR string and converts
    it to a list that splits up the numbers and characters'''
    # Input: CIGAR string (e.g., "6S64M")
    # Output: Parsed CIGAR string as list (e.g., ["6", "S", "64", "M"])

    initialize empty cigar str
    initialize empty cigar list

    for each char in cigar:
        if char is letter:
            append str to lst
            str = char
            append str to list
            str = ""
        else:
            append number to str
    
    return cigar list

def get_read_start(reverse, pos, cigar):
    '''Function that takes the strand, position, and CIGAR
    string, and calculates the 5' start of read'''
    # Inputs:
    # Boolean indicating if seq is reverse strand (e.g., 0 for forward)
    # Parsed CIGAR string (e.g., ["6", "S", "64", "M"])
    # 1-based leftmost matching position (e.g., 20)

    # Output: Read start (e.g., 14 for the inputs specified above)

    if (not reverse and S at cigar[1]):
        pos = pos - cigar[0]
    elif reverse:
        for 0, 2, 4, ... in cigar:
            pos += cigar
    retun pos

def dedup(infile, outfile, umi_set):
    ''' Main deduper logic, parses each line from input SAM file
    and writes unique sequences to output SAM file as they appear '''
    # Inputs: Input and output SAM file paths, UMI set
    # Output: None (unique SEQs written to output file)

    intialize curr_chrom to 0
    initialize curr_set to empty set

    open input sam file to read, open output sam file to write:
        read line, remove '\n'
        if line starts with '@':
            continue
        sam_line list = tab delimited line

        grab umi from sam_line
        if umi not in umi set:
            continue

        grab chrom from sam_line
        if chrom not curr_chrom:
            curr_chrom = chrom
            curr_set is empty
        
        grab pos from sam_line
        reverse = sam_line flag & 16 == 16

        grab cigar from sam_line
        parse cigar
        get read start from pos and cigar

        curr_seq = list of umi, chrom, reverse, read start
        if curr_seq not in curr_set:
            add to curr_set
            write line to output file
        
    return

def main():
    ''' Get argparse arguments, UMI set, and run dedup function '''
    # Input and Output: None
    
    get_args
    get_umis
    dedup(args, umi set)
```