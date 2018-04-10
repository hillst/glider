'''
cfdna glider performs the QC and interpretation of cfDNA datas.

This should be run on the cluster lol.
'''


def load_all_data(datatype="union"):
    from glider.data.annotations import load_data
    data = []
    data += load_data("BB053/trinuc_tsvs/BB053--BB1104." + datatype + ".tsv", biotype="control")
    data += load_data("BB053/trinuc_tsvs/BB053--BB1116." + datatype + ".tsv", biotype="control")
    data += load_data("BB053/trinuc_tsvs/BB053--BB1122." + datatype + ".tsv", biotype="control")
    data += load_data("BB053/trinuc_tsvs/BB053--BB1125." + datatype + ".tsv", biotype="control")
    data += load_data("BB053/trinuc_tsvs/BB053--BB1146." + datatype + ".tsv", biotype="control")
    data += load_data("BB342/trinuc_tsvs/BB342--BB1104." + datatype + ".tsv", biotype="control")
    data += load_data("BB342/trinuc_tsvs/BB342--BB1116." + datatype + ".tsv", biotype="control")
    data += load_data("BB342/trinuc_tsvs/BB342--BB1122." + datatype + ".tsv", biotype="control")
    data += load_data("BB342/trinuc_tsvs/BB342--BB1125." + datatype + ".tsv", biotype="control")
    data += load_data("BB342/trinuc_tsvs/BB342--BB1146." + datatype + ".tsv", biotype="control")
    data += load_data("BB427/trinuc_tsvs/BB427--BB1104." + datatype + ".tsv", biotype="control")
    data += load_data("BB427/trinuc_tsvs/BB427--BB1116." + datatype + ".tsv", biotype="control")
    data += load_data("BB427/trinuc_tsvs/BB427--BB1122." + datatype + ".tsv", biotype="control")
    data += load_data("BB427/trinuc_tsvs/BB427--BB1125." + datatype + ".tsv", biotype="control")
    data += load_data("BB427/trinuc_tsvs/BB427--BB1146." + datatype + ".tsv", biotype="control")
    data += load_data("BB600-cfDNA/trinuc_tsvs/BB600-cfDNA--BB1104." + datatype + ".tsv", biotype="control")
    data += load_data("BB600-cfDNA/trinuc_tsvs/BB600-cfDNA--BB1116." + datatype + ".tsv", biotype="control")
    data += load_data("BB600-cfDNA/trinuc_tsvs/BB600-cfDNA--BB1122." + datatype + ".tsv", biotype="control")
    data += load_data("BB600-cfDNA/trinuc_tsvs/BB600-cfDNA--BB1125." + datatype + ".tsv", biotype="control")
    data += load_data("BB600-cfDNA/trinuc_tsvs/BB600-cfDNA--BB1146." + datatype + ".tsv", biotype="control")
    data += load_data("BB601-cfDNA/trinuc_tsvs/BB601-cfDNA--BB1104." + datatype + ".tsv", biotype="control")
    data += load_data("BB601-cfDNA/trinuc_tsvs/BB601-cfDNA--BB1116." + datatype + ".tsv", biotype="control")
    data += load_data("BB601-cfDNA/trinuc_tsvs/BB601-cfDNA--BB1122." + datatype + ".tsv", biotype="control")
    data += load_data("BB601-cfDNA/trinuc_tsvs/BB601-cfDNA--BB1125." + datatype + ".tsv", biotype="control")
    data += load_data("BB601-cfDNA/trinuc_tsvs/BB601-cfDNA--BB1146." + datatype + ".tsv", biotype="control")
    data += load_data("BB672/trinuc_tsvs/BB672--BB1104." + datatype + ".tsv", biotype="control")
    data += load_data("BB672/trinuc_tsvs/BB672--BB1116." + datatype + ".tsv", biotype="control")
    data += load_data("BB672/trinuc_tsvs/BB672--BB1122." + datatype + ".tsv", biotype="control")
    data += load_data("BB672/trinuc_tsvs/BB672--BB1125." + datatype + ".tsv", biotype="control")
    data += load_data("BB672/trinuc_tsvs/BB672--BB1146." + datatype + ".tsv", biotype="control")
    data += load_data("Pre-BB1104/trinuc_tsvs/Pre-BB1104--BB1104." + datatype + ".tsv", biotype="luad")
    data += load_data("Pre-BB1116/trinuc_tsvs/Pre-BB1116--BB1116." + datatype + ".tsv", biotype="luad")
    data += load_data("Pre-BB1122/trinuc_tsvs/Pre-BB1122--BB1122." + datatype + ".tsv", biotype="luad")
    data += load_data("Pre-BB1125/trinuc_tsvs/Pre-BB1125--BB1125." + datatype + ".tsv", biotype="luad")
    data += load_data("Pre-BB1146/trinuc_tsvs/Pre-BB1146--BB1146." + datatype + ".tsv", biotype="luad")
    return data


def load_test():
    from glider.data.annotations import load_data, as_dataframe
    data = load_data("test_data/Pre-BB1104--BB1104.union.tsv")
    print len(data)
    print as_dataframe(data).describe()


print len(load_test())
