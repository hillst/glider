def load_data(filename, MAX=-1, biotype="luad", batch=None):
    """
    Loads data into vectorized form.
    """
    from glider.utils.trinucs import vectorize_row
    import os
    sample_name = os.path.basename(filename).strip(".tsv")
    read_source, mutation_source = sample_name.split("--")
    mutation_sample, mutation_callset = mutation_source.split(".")

    data = []
    for i, line in enumerate(open(filename, 'r')):
        if MAX > 0 and i > MAX:
            break
        line_arr = line.strip().split("\t")
        chrm, pos, ref, alt, read_id, mrbq, vbq, mapq, pir, frag, trinuc, trinuc_index = line_arr
        data.append(vectorize_row(mrbq, vbq, mapq, pir, frag, trinuc_index, normalize=False, reads=read_source,
                                  mutations=mutation_sample, callset=mutation_callset, biotype=biotype))

    return data


def load_all_data():
    """
    lol this is supposed to be static so ?

    :return:
    """
    pass
