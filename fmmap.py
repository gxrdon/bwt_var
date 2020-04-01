#!/usr/bin/env python
import math
from simplesam import Writer, Sam


def readfq(fp): # this is a generator function
    """

    readfq from https://github.com/lh3/readfq/blob/master/readfq.py

    """
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
    genome_with_stop = genome + "$"
    output_file = open(reference_idx, "a+")
    output_file.write(str(arr[0]) + " " + str(len(genome)) + "\n")
    output_file.write(genome + "\n")
    last_col = burrows_wheeler(genome_with_stop)
    output_file.write(last_col + "\n")
    first_col = sorted(list(last_col))
    first = ""
    for i in range(0, len(first_col)):
        first += first_col[i]
    output_file.write(str(first) + "\n")
    suff_arr = construct_suffix_array(genome_with_stop)
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


def get_first_int(first_char, first, occur):
    if first_char == "A":
        return 1, 1 + int(occur[len(first) - 2][0])
    elif first_char == "C":
        return 1 + int(occur[len(first) - 2][0]), 1 + int(occur[len(first) - 2][0]) + int(occur[len(first) - 2][1])
    elif first_char == "G":
        return 1 + int(occur[len(first) - 2][0]) + int(occur[len(first) - 2][1]), 1 + int(occur[len(first) - 2][0]) + int(occur[len(first) - 2][1]) + int(occur[len(first) - 1][2])
    else:
        return 1 + int(occur[len(first) - 2][0]) + int(occur[len(first) - 2][1]) + int(occur[len(first) - 2][2]), 1 + int(occur[len(first) - 2][0]) + int(occur[len(first) - 2][1]) + int(occur[len(first) - 2][2]) + int(occur[len(first) - 2][3])


def get_curr_int(first_char, first, occur, start_offset, end_offset):
    if first_char == "A":
        return 1 + start_offset, 1 + end_offset
    elif first_char == "C":
        return 1 + occur[len(first) - 2][0] + start_offset, 1 + occur[len(first) - 2][0] + end_offset
    elif first_char == "G":
        return 1 + occur[len(first) - 2][0] + occur[len(first) - 2][1] + start_offset, 1 + occur[len(first) - 2][0] + occur[len(first) - 2][1] + end_offset
    else:
        return 1 + occur[len(first) - 2][0] + occur[len(first) - 2][1] + occur[len(first) - 2][2] + start_offset, 1 + occur[len(first) - 2][0] + occur[len(first) - 2][1] + occur[len(first) - 2][2] + end_offset


def find_occ_pos(s):
    if s == 'A':
        return 0
    elif s == 'C':
        return 1
    elif s == 'G':
        return 2
    else:
        return 3


def get_interval_fm(read, first, last, occur):
    int_start, int_end = get_first_int(read[0], first, occur)
    for i in range(1, len(read)):
        pos = find_occ_pos(read[i])
        p = int(int_start) - 1
        start_offset = occur[p][pos]
        end_offset = occur[int_end][pos]
        int_start, int_end = get_curr_int(read[i], first, occur, start_offset, end_offset)

    return [int_start, int_end]


def find_alignment(dp, read, genome, gap, back_pointer):
    backward_str = ""
    y = len(read)
    start_pos = 0
    curr_max = float("-inf")
    count = 0
    for i in range(0, len(genome)):
        if dp[i][y] >= curr_max:
            curr_max = dp[i][y]
            start_pos = i
    x = start_pos
    num_ins_or_match = 0
    count_match = 0
    cont = True
    while cont:
        curr_x, curr_y = back_pointer[x][y]
        if curr_x == x and curr_y == y - 1:
            backward_str = '+' + backward_str
            y -= 1
            num_ins_or_match += 1
        elif curr_x == x - 1 and curr_y == y - 1:
            backward_str = read[y - 1] + backward_str
            x -= 1
            y -= 1
            num_ins_or_match += 1
            count_match += 1
        else:
            backward_str = '-' + backward_str
            x -= 1
        count += 1
        if y == 0 or num_ins_or_match == 100:
            cont = False

    # May help us speed up the program slightly
    perfect_match = False
    if count_match == 100:
        perfect_match = True

    return backward_str, perfect_match, x


def fitting_alignment(read, reference, gap):
    dp = [[0 for i in range(0, len(read) + 1)] for j in range(0, len(reference) + 1)]
    backwards_pointer = [[(0, 0) for i in range(0, len(read) + 1)] for j in range(0, len(reference) + 1)]
    for i in range(0, len(read) + 1):
        if i == 0:
            backwards_pointer[0][i] = (0, 0)
        else:
            backwards_pointer[0][i] = (0, i - 1)
        dp[0][i] = i * (-1 * gap)

    for i in range(1, len(reference) + 1):
        backwards_pointer[i][0] = (i - 1, 0)

    max_return = float("-inf")
    for i in range(1, len(reference) + 1):
        for j in range(1, len(read) + 1):
            cost = 0
            if reference[i - 1] != read[j - 1]:
                cost = 2
            curr_max = max(dp[i - 1][j - 1] - cost, dp[i - 1][j] - gap, dp[i][j - 1] - gap)
            dp[i][j] = curr_max
            if curr_max == dp[i - 1][j - 1] - cost:
                backwards_pointer[i][j] = (i - 1, j - 1)
            elif curr_max == dp[i - 1][j] - gap:
                backwards_pointer[i][j] = (i - 1, j)
            else:
                backwards_pointer[i][j] = (i, j - 1)
            if j == len(read) and dp[i][j] >= max_return:
                max_return = dp[i][j]

    return find_alignment(dp, read, reference, gap, backwards_pointer), max_return


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
    last = list(idx_file[2])
    first = list(idx_file[3])
    suff_arr = idx_file[4]
    suff_arr = [int(i) for i in suff_arr[1: len(suff_arr) - 1].split(', ')]
    occur = idx_file[5]
    occur = occur[2:len(occur) - 2]
    occurence_list = occur.split("], [")
    temp_list = [0 for i in range(0, 29822)]
    for i in range(0, len(occurence_list)):
        temp = (str(occurence_list[i])[0: len(occurence_list[i])]).split(", ")
        for j in range(0, len(temp)):
            if j == 0:
                temp_list[i] = [int(temp[j])]
            else:
                temp_list[i].append(int(temp[j]))
    occur = temp_list
    genome_len = len(genome)
    ninf = float("-inf")
    seed_skip = lambda l: math.floor(l / 5.0)
    gap = 2
    out_file = open(output_file, "a+")
    writer = Writer(out_file)
    out_file.write("@SQ\tSN:" + genome_reference_name + "\tLN:" + genome_reference_len + "\n")

    with open(reads) as f:
        flag = True
        for name, seq, qual in readfq(f):
            if "N" in seq:
                continue
            alignments = []
            read_len = len(seq)
            best_score = ninf
            seed_pos = 0
            skip = 20
            visited_dict = {}
            perfect_match = False
            for seed_start in range(0, read_len, skip):
                seed_end = min(read_len, seed_start + skip)

                # get an interval from the bwt transform that matches the slice
                curr_read = seq[seed_start:seed_end]
                interval = get_interval_fm(curr_read[len(curr_read)::-1], first, last, occur)

                # list of slices from genome with padding where we may have exact match
                reference, pos_list = ref_positions(suff_arr, interval, genome, seed_start, seed_end)
                for i in range(0, len(pos_list)):
                    # Eliminates function calls to same reference
                    if pos_list[i] in visited_dict:
                        continue
                    else:
                        visited_dict[pos_list[i]] = True
                    match_ret, score = fitting_alignment(seq, reference[i], gap)
                    alignment, perf_match, position = match_ret
                    a = (alignment, score, position, pos_list[i])
                    if score > best_score:
                        best_score = score
                        alignments = [a]
                    elif score == best_score:
                        if a not in alignments:
                            alignments.append(a)

                    if perf_match:
                        alignments = [a]
                        perfect_match = True
                if perfect_match:
                    break

            # Linear operation: shouldn't be a time sink
            for a in alignments:
                print(name)
                sam_record = Sam()
                ali, score, position, start_pos = a
                fix_index = 1
                if start_pos + 105 < genome_len:
                    fix_index += position
                sam_record.flag = 0
                sam_record.tlen = 0
                sam_record.pos = start_pos + fix_index
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
