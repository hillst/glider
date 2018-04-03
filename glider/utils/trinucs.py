def build_trinuc_lookup():
    '''
    Builds a trinuc-to-index and index-to-trinuc dictionaries for decoding indicators values and trinucleotide strings.
    :return: trinuc: index, index: trinuc
    '''
    values = []
    index = 0
    for c in "ACGT":
        for d in "ACGT":
            for ref in "CT":
                for alt in 'ACGT':
                    if ref == alt:
                        continue
                    value = "{}[{}>{}]{}".format(c, ref, alt, d)
                    values.append(value)
                    index += 1
    return {v: i for i, v in enumerate(values)}, {i: v for i, v in enumerate(values)}


def trinuc_order():
    '''
    Static trinucleotide order used by COSMIC.
    :return:
    '''
    from functools import cmp_to_key
    return sorted(build_trinuc_lookup()[0].keys(), key=trinuc_cmp)


def trinuc_colors():
    lookup = {"CA": "lightblue", "CG": "black", "CT": "red", "TA": "lightgrey", "TC": "forestgreen", "TG": "magenta"}
    return [lookup[trinuc[2] + trinuc[4]] for trinuc in trinuc_order()]


def invert_phred(qscore):
    '''
    Invers the PHRED valued q-score to a probability of sequencing error.

    :param qscore:
    :return: float [0-1]
    '''
    p = 10 ** ((-1 * float(qscore)) / 10)
    return p


def vectorize_row(mrbq, vbq, mapq, pir, frag, trinuc_index, normalize=False, reads="", mutations="", callset="",
                  biotype=""):
    '''
    Vectorizes the combination of QC sample statistics and creates a categorical vector for the trinucleotides.

    :param mrbq:
    :param vbq:
    :param mapq:
    :param pir:
    :param frag:
    :param trinuc_index:
    :param normalize: Scales all the values to be in: [0,1]. Useful for modelling.
    :param reads:
    :param mutations:
    :param callset:
    :param biotype:
    :return:
    '''
    if normalize:
        return [invert_phred(mrbq), invert_phred(vbq), float(mapq) / 60., abs(75 - int(pir)) / 151.,
                abs(float(frag)) / 300.] + [1 if i == int(trinuc_index) else 0 for i in range(96)] + [reads, mutations,
                                                                                                      callset, biotype]
    else:
        return [float(mrbq), float(vbq), float(mapq), float(pir), abs(float(frag))] + [
            1 if i == int(trinuc_index) else 0 for i in range(96)] + [reads, mutations, callset, biotype]


def trinuc_cmp(A):
    '''
    Comparison function for comparing two tri-nucleotides of the format A[C>T]A
    :param A:
    :return:
    '''
    return "".join((A[2], A[4], A[0], A[6]))
