import gzip
import numpy
from scipy.stats import norm


def parse_file(file_input):
    transcripts = {}
    transcript_list = []
    alignments = []
    alignment_bool = False
    num_transcripts = 0

    with gzip.open(file_input, 'rb') as f:
        count = 0
        align_idx = -1
        for line in f:
            temp = str(line).split('\'')
            if count == 0:
                num_transcripts = int(temp[1][0: len(temp[1]) - 2])
                count += 1
            elif alignment_bool:
                if len(str(line).split("\\t")) == 1:
                    alignments.append([])
                    align_idx += 1
                else:
                    curr_align = temp[1].split("\\t")
                    align_name = str(curr_align[0])
                    align_ori = str(curr_align[1])
                    align_pos = int(curr_align[2])
                    align_prob = curr_align[3][0: len(curr_align[3]) - 2].split("e")

                    if len(align_prob) > 1:
                        num = float(align_prob[0])
                        power = float(align_prob[1])
                        align_prob = num ** power
                    else:
                        align_prob = float(align_prob[0])
                    align_obj = {"index": transcripts[align_name], "ori": align_ori, "pos": align_pos, "prob": align_prob}
                    alignments[align_idx].append(align_obj)
            else:
                count += 1
                if count > num_transcripts:
                    print("done transcript parsing")
                    alignment_bool = True
                n = temp[1].split("\\t")
                key = n[0]
                val = int(n[1][0: len(n[1]) - 2])
                transcripts[key] = count - 2
                trans_obj = {"length": val, "name": key}
                transcript_list.append(trans_obj)
    print("done alignment parsing")
    return num_transcripts, transcript_list, alignments


def get_effective_length(norm_len):
    mu = 200
    sd = 25

    d = norm(mu, sd)
    p = numpy.array([d.pdf(i) for i in range(norm_len)])
    p /= p.sum()
    cond_mean = numpy.sum([i * p[i] for i in range(len(p))])
    eff_len = norm_len - cond_mean
    return eff_len


def full_model_em(align_file, output_file, num_transcripts, transcripts, alignments):
#def full_model_em(align_file, output_file):
#    num_transcripts, transcripts, alignments = parse_file(align_file)

    n_arr = [1/num_transcripts] * num_transcripts
    not_converged = True

    while not_converged:
        n_update = [0] * num_transcripts
        not_converged = False
        print("Full model not converged")
        for align in alignments:
            # Short circuit if only one transcript instead of running EM on it
            if len(align) == 1:
                n_update[align[0]["index"]] += 1
                continue
            for trans in align:
                p_1 = n_arr[trans["index"]]
                p_2 = 1 / get_effective_length(transcripts[trans["index"]]["length"])
                p_3 = norm.cdf(transcripts[trans["index"]]["length"] - trans["pos"])
                p_4 = trans["prob"]
                n_update[trans["index"]] += (p_1 * p_2 * p_3 * p_4)
                print("Name: " + str(transcripts[trans["index"]]["name"]) + " Score: " + str(n_update[trans["index"]]))

        # If we find that one value hasn't converged, continue process
        for i in range(0, num_transcripts):
            if int(n_arr[i] * 1000) != int(n_update[i] * 1000):
                not_converged = True
        n_arr = n_update
    print("Writing to full EM file")
    f = open(output_file, "w")
    for i in range(0, num_transcripts):
        f.write(str(transcripts[i]["name"]) + "\t" + str(get_effective_length(transcripts[i]["length"])) + "\t" + str(n_arr[i]))
    f.close()


def equivalence_class_em(align_file, output_file, num_transcripts, transcripts, alignments):
#def equivalence_class_em(align_file, output_file):
#    num_transcripts, transcripts, alignments = parse_file(align_file)

    eq_classes = {}

    for align in alignments:
        tset = tuple(sorted(a["index"] for a in align))
        if tset in eq_classes:
            eq_classes[tset] += 1
        else:
            eq_classes[tset] = 1

    n_arr = [1 / num_transcripts] * num_transcripts
    not_converged = True

    while not_converged:
        n_update = [0] * num_transcripts
        not_converged = False
        print("Equiv not converged")
        count = 0
        for eq_class in eq_classes:
            print("Iteration: " + str(count))
            # Short circuit if only one transcript instead of running EM on it
            if len(eq_class) == 1:
                n_update[eq_class[0]] += 1
                continue
            summation_value = 0
            for i in range(0, len(eq_class)):
                summation_value += n_arr[eq_class[i]] * (1 / get_effective_length(transcripts[eq_class[i]]["length"]))
            for trans in eq_class:
                p_1 = n_arr[trans]
                p_2 = 1 / get_effective_length(transcripts[trans]["length"])
                n_update[trans] += (eq_classes[eq_class] * p_1 * p_2) / summation_value

        # If we find that one value hasn't converged, continue process
        for i in range(0, len(eq_classes)):
            if int(n_arr[i] * 1000) != int(n_update[i] * 1000):
                not_converged = True
        n_arr = n_update

    print("Writing to equiv EM file")
    f = open(output_file, "w")
    for i in range(0, num_transcripts):
        f.write(str(transcripts[i]["name"]) + "\t" + str(get_effective_length(transcripts[i]["length"])) + "\t" + str(n_arr[i]))
        print()
    f.close()


def call_both_em_functions(align_file, equiv_output, full_output):
    """
    Function created to run both the functions at once so we don't have to parse the file twice.
    """
    num_transcripts, transcripts, alignments = parse_file(align_file)
    equivalence_class_em(align_file, equiv_output, num_transcripts, transcripts, alignments)
    full_model_em(align_file, full_output, num_transcripts, transcripts, alignments)


#full_model_em("./data/alignments.cs423.gz", "./data/output.txt")
call_both_em_functions("./data/alignments.cs423.gz", "./data/equiv_output.txt", "./data/full_output.txt")
"""if __name__ == "__main__":
    import sys
    version = sys.argv[1]
    if version == "--eqc":
        equivalence_class_em(sys.argv[2], sys.argv[3])
    else:
        full_model_em(sys.argv[1], sys.argv[2])
"""