filters = [
    ('all', lambda x: x),
    ('coding', lambda x: x[ x['location'].str.match('coding')]),
    ('splicesite', lambda x: x[ x['location'].str.match('splicesite')]),
    ('5primeUTR', lambda x: x[ x['location'].str.match('5primeUTR')]),
    ('3primeUTR', lambda x: x[ x['location'].str.match('3primeUTR')]),
    ('deepintronic', lambda x: x[ x['location'].str.match('deepintronic')]),
    ('noncodingRNA', lambda x: x[ x['location'].str.match('noncodingRNA')]),
    ('mithocondrial', lambda x: x[x['location'].str.match('mithocondrial')]),
    ('unknown', lambda x: x[ x['location'].str.match('unknown')])
]

filters_var_type = [
    ('all_types', lambda x: x),
    ('snps', lambda x: x[x['type'].str.match('snp')]),
    ('indels', lambda x: x.query('type == "indel" & subtype == "ins" | type == "indel" & subtype == "del"')),
    ('insertions', lambda x:  x.query('type == "indel" & subtype == "ins"')),
    ('deletions', lambda x: x.query('type == "indel" & subtype == "del"')),
    ('mnps', lambda x: x.query('type == "indel" & subtype == "mnp"'))
]

filter_intronic_bins = [
    ('all_intronic', lambda x: x[~x['intron_bin'].isnull()]),
    ('all_except_0-10', lambda x: x[~x['intron_bin'].str.match('0-10')]),
    ('0-10', lambda x: x[x['intron_bin'].str.match('0-10')]),
    ('10-30', lambda x: x[x['intron_bin'].str.match('10-30')]),
    # ('30-60', lambda x: x[x['intron_bin'].str.match('30-60')]),
    # ('60-100', lambda x: x[x['intron_bin'].str.match('60-100')]),
    ('30-100', lambda x: x[x['intron_bin'].str.match('30-100')]),
    ('100-200', lambda x: x[x['intron_bin'].str.match('100-200')]),
    ('200-500', lambda x: x[x['intron_bin'].str.match('200-500')]),
    ('500+', lambda x: x[x['intron_bin'].str.match('500\+')])
]


def subset_variants_by_type(filters_var_type, types_to_analyse):

    if not types_to_analyse:
        return filters_var_type
    else:
        return [var_types for var_types in filters_var_type if var_types[0] in types_to_analyse]
