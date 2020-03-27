#!/usr/bin/env python
import math
from simplesam import Writer, Sam


def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>&#64;': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '&#64;+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break


def construct_suffix_array(suffix):
    suffix_array = sorted([(suffix[i:], i) for i in range(len(suffix))])
    temp_arr = []
    for i in range(0, len(suffix_array)):
        temp_arr.append(suffix_array[i][1])
    return temp_arr


def burrows_wheeler(string):
    input_array = list(string)
    array = []
    for i in range(len(input_array)):
        temp_string = string[1: len(string)] + string[0]
        string = ''.join(temp_string)
        array.append(string)
    sorted_array = sorted(array)
    burrows_wheeler_string = ""
    for i in range(len(sorted_array)):
        burrows_wheeler_string += sorted_array[i][-1]
    return burrows_wheeler_string


def occ(string):
    occ_arr = [[0 for i in range(4)] for j in range(len(string))]
    for i in range(0, len(string)):
        if string[i] == 'A':
            if i == 0:
                occ_arr[i][0] = 1
            else:
                occ_arr[i][0] = occ_arr[i - 1][0] + 1
                occ_arr[i][1] = occ_arr[i - 1][1]
                occ_arr[i][2] = occ_arr[i - 1][2]
                occ_arr[i][3] = occ_arr[i - 1][3]
        elif string[i] == 'C':
            if i == 0:
                occ_arr[i][1] = 1
            else:
                occ_arr[i][1] = occ_arr[i - 1][1] + 1
                occ_arr[i][0] = occ_arr[i - 1][0]
                occ_arr[i][2] = occ_arr[i - 1][2]
                occ_arr[i][3] = occ_arr[i - 1][3]
        elif string[i] == 'G':
            if i == 0:
                occ_arr[i][2] = 1
            else:
                occ_arr[i][2] = occ_arr[i - 1][2] + 1
                occ_arr[i][1] = occ_arr[i - 1][1]
                occ_arr[i][0] = occ_arr[i - 1][0]
                occ_arr[i][3] = occ_arr[i - 1][3]
        elif string[i] == 'T':
            if i == 0:
                occ_arr[i][3] = 1
            else:
                occ_arr[i][3] = occ_arr[i - 1][3] + 1
                occ_arr[i][1] = occ_arr[i - 1][1]
                occ_arr[i][2] = occ_arr[i - 1][2]
                occ_arr[i][0] = occ_arr[i - 1][0]
        else:
            occ_arr[i][1] = occ_arr[i - 1][1]
            occ_arr[i][2] = occ_arr[i - 1][2]
            occ_arr[i][0] = occ_arr[i - 1][0]
            occ_arr[i][3] = occ_arr[i - 1][3]
    return occ_arr


def index(reference_file, reference_idx):
    f = open(reference_file)
    read = f.read()
    arr = read.split('\n')
    f.close()
    genome = ""
    for i in range(1, len(arr) - 1):
        genome += arr[i]

    output_file = open(reference_idx, "a+")
    output_file.write(str(arr[0]) + " " + str(len(genome)) + "\n")
    output_file.write(genome + "\n")
    last_col = list(burrows_wheeler(genome))
    output_file.write(str(last_col) + "\n")
    first_col = sorted(last_col)
    output_file.write(str(first_col) + "\n")
    suff_arr = construct_suffix_array(genome)
    output_file.write(str(suff_arr) + "\n")
    occ_arr = occ(last_col)
    output_file.write(str(occ_arr))
    output_file.close()


def ref_positions(suff_arr, interval, genome, seed_start, seed_end):
    padded_list = []
    pos_list = []
    if len(interval) != 2:
        return [], [0]
    for i in range(interval[0], interval[1] + 1):
        start = 0
        end = len(genome) - 1
        seed_st = suff_arr[i] - (seed_start + 5)
        seed_fi = suff_arr[i] + (125 - seed_end)
        if seed_st > 0:
            start = seed_st
        else:
            end = 110
        if seed_fi < end:
            end = seed_fi
        else:
            start = end - 110
        pos_list.append(start)
        padded_list.append(genome[start: end])
    return padded_list, pos_list


def get_interval(arr, genome, read):
    i = 0
    interval = []
    flag = False
    while genome[arr[i]] != read[0]:
        i += 1
        if i == len(genome) - 1:
            return []
    while genome[arr[i]] == read[0]:
        if arr[i] + len(read) > len(genome) and read[0:len(genome) - arr[i]] == genome[arr[i]:len(genome)]:
            if len(interval) == 0:
                flag = True
                interval = [i]
        elif arr[i] + len(read) <= len(genome) and read == genome[arr[i]: arr[i] + len(read)]:
            if len(interval) == 0:
                flag = True
                interval = [i]
        elif flag:
            break
        if i == len(genome) - 1:
            break
        i += 1
    interval.append(i - 1)
    return interval


def find_alignment(dp, read, genome, gap, back_pointer):
    backward_str = ""
    y = len(read)
    start_pos = 0
    curr_max = float("-inf")
    count = 0
    for i in range(0, len(genome)):
        if dp[i][len(read)] >= curr_max:
            curr_max = dp[i][len(read)]
            start_pos = i
    x = start_pos

    cont = True
    while cont:
        curr_x, curr_y = back_pointer[x][y]
        if curr_x == x and curr_y == y - 1:
            backward_str = '-' + backward_str
            y -= 1
        elif curr_x == x - 1 and curr_y == y - 1:
            backward_str = read[y - 1] + backward_str
            x -= 1
            y -= 1
        else:
            backward_str = '+' + backward_str
            x -= 1
        count += 1
        if y == 0 or x == 0 or count == 100:
            cont = False

    return backward_str


def fitting_alignment(read, reference, gap):
    dp = [[0 for i in range(0, len(read) + 1)] for j in range(0, len(reference) + 1)]
    backwards_pointer = [[(0, 0) for i in range(0, len(read) + 1)] for j in range(0, len(reference) + 1)]

    for i in range(0, len(read) + 1):
        if i == 0:
            backwards_pointer[0][i] = (0, 0)
        else:
            backwards_pointer[0][i] = (i - 1, 0)
        dp[0][i] = i * (-1 * gap)

    for i in range(1, len(reference) + 1):
        for j in range(1, len(read) + 1):
            cost = 0
            if reference[i - 1] != read[j - 1]:
                cost = -2
            curr_max = max(cost + dp[i - 1][j - 1], dp[i - 1][j] - gap, dp[i][j - 1] - gap)
            dp[i][j] = curr_max
            if curr_max == cost + dp[i - 1][j - 1]:
                backwards_pointer[i][j] = (i - 1, j - 1)
            elif curr_max == dp[i - 1][j] - gap:
                backwards_pointer[i][j] = (i - 1, j)
            else:
                backwards_pointer[i][j] = (i, j - 1)

    return find_alignment(dp, read, reference, gap, backwards_pointer), dp[len(reference)][len(read)], reference


def create_cigar(seq):
    ins_count = 0
    del_count = 0
    match_count = 0
    count = 0
    cigar = ""

    for i in range(0, len(seq)):
        count += 1
        if seq[i] == '+':
            ins_count += 1
            if del_count != 0:
                cigar += str(del_count) + "D"
                del_count = 0
            elif match_count != 0:
                cigar += str(match_count) + "M"
                match_count = 0
        elif seq[i] == '-':
            del_count += 1
            if ins_count != 0:
                cigar += str(ins_count) + "I"
                ins_count = 0
            elif match_count != 0:
                cigar += str(match_count) + "M"
                match_count = 0
        else:
            match_count += 1
            if del_count != 0:
                cigar += str(del_count) + "D"
                del_count = 0
            elif ins_count != 0:
                cigar += str(ins_count) + "I"
                ins_count = 0

    if del_count != 0:
        cigar += str(del_count) + "D"
    if ins_count != 0:
        cigar += str(ins_count) + "I"
    if match_count != 0:
        cigar += str(match_count) + "M"

    return cigar


def align(reference_idx, reads, output_file):
    f = open(reference_idx)
    read_file = f.read()
    idx_file = read_file.split('\n')
    f.close()

    genome_name = (idx_file[0]).split(" ")
    genome_reference_name = genome_name[0][1:len(genome_name[0])]
    genome_reference_len = genome_name[len(genome_name) - 1]
    genome = idx_file[1]
    last = idx_file[2]
    first = idx_file[3]
    suff_arr = idx_file[4]
    suff_arr = [int(i) for i in suff_arr[1: len(suff_arr) - 1].split(', ')]
    occ = idx_file[5]

    ninf = float("-inf")
    seed_skip = lambda l: math.floor(l / 5.0)
    gap = 2
    out_file = open(output_file, "a+")
    writer = Writer(out_file)
    out_file.write("@SQ\tSN:" + genome_reference_name + "\tLN:" + genome_reference_len + "\n")

    with open(reads) as f:
        flag = True
        for name, seq, qual in readfq(f):
            alignments = []
            read_len = len(seq)
            best_score = ninf
            seed_pos = 0
            skip = seed_skip(read_len)
            for seed_start in range(0, read_len, skip):
                seed_end = min(read_len, seed_start + skip)
                # get an interval from the bwt transform that matches the slice
                interval = get_interval(suff_arr, genome, seq[seed_start:seed_end])

                # list of slices from genome with padding where we may have exact match
                reference, sequence = ref_positions(suff_arr, interval, genome, seed_start, seed_end)
                for i in range(0, len(reference)):

                    alignment, score, position = fitting_alignment(seq, reference[i], gap)
                    a = (alignment, score, position, sequence[i])
                    if score > best_score:
                        best_score = score
                        alignments = [a]
                    elif score == best_score:
                        if a not in alignments:
                            alignments.append(a)

            for a in alignments:
                print(name)
                sam_record = Sam()
                if flag:
                    sam_record.flag = 0
                    flag = False
                else:
                    sam_record.flag = 256
                ali, score, sequ, start_pos = a
                sam_record.tlen = 0
                sam_record.pos = start_pos
                sam_record.rnext = '*'
                sam_record.pnext = 0
                sam_record.mapq = 255
                sam_record.qual = '*'
                sam_record.qname = name
                sam_record.rname = genome_reference_name
                sam_record.cigar = create_cigar(ali)
                sam_record.seq = seq
                writer.write(sam_record)
    out_file.close()


if __name__ == "__main__":
    import sys
    func = sys.argv[1]
    if func == "index":
        index(sys.argv[2], sys.argv[3])
    elif func == "align":
        align(sys.argv[2], sys.argv[3], sys.argv[4])
    else:
        print("Misformed input")


#print(construct_suffix_array("AAATAAGTAGAAT$"))
#print(get_interval(construct_suffix_array("AAATAAGTAGAAT$"), "AAATAAGTAGAAT$", "AA"))
#index('./data/2019-nCoV.fa', './data/test_output.txt')
#align('./data/test_output.txt', './data/reads.fa', './data/temp.sam')
# print("Done")
