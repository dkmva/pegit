from itertools import zip_longest
import json
import mmap
import os.path
import subprocess

import django.conf
import regex


DEFAULT_CODON_TABLE = {
    'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAT': 'N', 'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AGA': 'R', 'AGC': 'S', 'AGG': 'R', 'AGT': 'S', 'ATA': 'I', 'ATC': 'I', 'ATG': 'M', 'ATT': 'I',
    'CAA': 'Q', 'CAC': 'H', 'CAG': 'Q', 'CAT': 'H', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R', 'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'GAA': 'E', 'GAC': 'D', 'GAG': 'E', 'GAT': 'D', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G', 'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'TAA': '$', 'TAC': 'Y', 'TAG': '$', 'TAT': 'Y', 'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TGA': '$', 'TGC': 'C', 'TGG': 'W', 'TGT': 'C', 'TTA': 'L', 'TTC': 'F', 'TTG': 'L', 'TTT': 'F'
}


DEGENERATE_CODON_TABLE = {
    'A': 'GCN', 'C': 'TGY', 'E': 'GAR', '$': 'TRR', 'G': 'GGN', 'F': 'TTY', 'I': 'ATH', 'H': 'CAY',
    'K': 'AAR', 'M': 'ATG', 'L': 'YTN', 'N': 'AAY', 'Q': 'CAR', 'P': 'CCN', 'S': 'WSN', 'R': 'MGN',
    'T': 'ACN', 'W': 'TGG', 'V': 'GTN', 'Y': 'TAY', 'D': 'GAY'
}

degenerate_to_nucleotides = {
    'A': {'A'}, 'C': {'C'}, 'G': {'G'}, 'T': {'T'},
    'M': {'A', 'C'}, 'R': {'A', 'G'}, 'W': {'A', 'T'}, 'S': {'C', 'G'}, 'Y': {'C', 'T'}, 'K': {'G', 'T'},
    'V': {'A', 'C', 'G'}, 'H': {'A', 'C', 'T'}, 'D': {'A', 'G', 'T'}, 'B': {'C', 'G', 'T'},
    'N': {'A', 'C', 'G', 'T'}
}

# HUMAN CODON TALBE https://www.genscript.com/tools/codon-frequency-table
DEFAULT_CODON_USAGE_TABLE = {
    'TTT': ['F', 0.45], 'TTC': ['F', 0.55], 'TTA': ['L', 0.07], 'TTG': ['L', 0.13],
    'TCT': ['S', 0.18], 'TCC': ['S', 0.22], 'TCA': ['S', 0.15], 'TCG': ['S', 0.06],
    'TAT': ['Y', 0.43], 'TAC': ['Y', 0.57], 'TAA': ['$', 0.28], 'TAG': ['$', 0.20],
    'TGT': ['C', 0.45], 'TGC': ['C', 0.55], 'TGA': ['$', 0.52], 'TGG': ['W', 1.00],

    'CTT': ['L', 0.13], 'CTC': ['L', 0.20], 'CTA': ['L', 0.07], 'CTG': ['L', 0.41],
    'CCT': ['P', 0.28], 'CCC': ['P', 0.33], 'CCA': ['P', 0.27], 'CCG': ['P', 0.11],
    'CAT': ['H', 0.41], 'CAC': ['H', 0.59], 'CAA': ['Q', 0.25], 'CAG': ['Q', 0.75],
    'CGT': ['R', 0.08], 'CGC': ['R', 0.19], 'CGA': ['R', 0.11], 'CGG': ['R', 0.21],

    'ATT': ['I', 0.36], 'ATC': ['I', 0.48], 'ATA': ['I', 0.16], 'ATG': ['M', 1.00],
    'ACT': ['T', 0.24], 'ACC': ['T', 0.36], 'ACA': ['T', 0.28], 'ACG': ['T', 0.12],
    'AAT': ['N', 0.46], 'AAC': ['N', 0.54], 'AAA': ['K', 0.42], 'AAG': ['K', 0.58],
    'AGT': ['S', 0.15], 'AGC': ['S', 0.24], 'AGA': ['R', 0.2], 'AGG': ['R', 0.2],

    'GTT': ['V', 0.18], 'GTC': ['V', 0.24], 'GTA': ['V', 0.11], 'GTG': ['V', 0.47],
    'GCT': ['A', 0.26], 'GCC': ['A', 0.40], 'GCA': ['A', 0.23], 'GCG': ['A', 0.11],
    'GAT': ['D', 0.46], 'GAC': ['D', 0.54], 'GAA': ['E', 0.42], 'GAG': ['E', 0.58],
    'GGT': ['G', 0.16], 'GGC': ['G', 0.34], 'GGA': ['G', 0.25], 'GGG': ['G', 0.25],
}


def load_codon_table(name):
    """Load a codon usage table from a json file."""
    try:
        with open(os.path.join(django.conf.settings.DESIGN_CODON_USAGE_FOLDER, f'{name}.json'), 'r') as f:
            return json.load(f)
    except FileNotFoundError:
        return DEFAULT_CODON_USAGE_TABLE


def dump_codon_table(codon_table, name):
    """Write a codon usage table to a json file."""
    with open(os.path.join(django.conf.settings.DESIGN_CODON_USAGE_FOLDER, f'{name}.json'), 'w') as f:
        json.dump(codon_table, f)


def is_valid_dna(seq):
    """Checks if sequence is a valid DNA sequence."""
    valid = 'actgACTG'
    return all(c in valid for c in seq)


def is_valid_degenerate(seq):
    """Checks if sequence is a valid degenerate DNA sequence."""
    valid = 'abcdghkmnrstvwyABCDGHKMNRSTVWY'
    return all(c in valid for c in seq)


COMPLEMENT = str.maketrans('ATGCRYMKBVDHatgcrymkbvdh',
                           'TACGYRKMVBHDtacgyrkmvbhd')


def complement(seq):
    """Returns the complement of a sequence."""
    try:
        return seq.complement()
    except AttributeError:
        return seq.translate(COMPLEMENT)


def reverse_complement(seq):
    """Returns the reverse complement of a sequence."""
    try:
        return seq.reverse_complement()
    except AttributeError:
        return complement(seq)[::-1]


def translate(sequence, codon_table=None):
    """Translates a DNA sequence into a protein sequence."""
    if codon_table is None:
        codon_table = DEFAULT_CODON_TABLE
    peptide = list()
    for codon in zip_longest(*[iter(sequence.upper())] * 3, fillvalue=None):
        try:
            amino_acid = codon_table[''.join(codon)]
            peptide.append(amino_acid)
        except TypeError:
            break
    return ''.join(peptide)


def gc(seq):
    """Calculate GC content of sequence."""
    if not len(seq):
        return 0
    return (seq.count('C') + seq.count('G')) / len(seq)


def dgn_to_regex(dgn):
    """Convert a degenerate sequence to regex."""
    ret = []
    for n in dgn:
        n = ''.join(degenerate_to_nucleotides[n.upper()])
        ret.append(f'[{n}]')

    return ''.join(ret)


def get_num_lines(file_path):
    fp = open(file_path, "r+")
    buf = mmap.mmap(fp.fileno(), 0)
    lines = 0
    while buf.readline():
        lines += 1
    return lines


def parse_gff_file(file_path, gff_record_fields=None, disable=False):
    import tqdm

    if gff_record_fields is None:
        gff_record_fields = ["chromosome", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]

    with open(file_path, 'r') as f:
        for line in tqdm.tqdm(f, total=get_num_lines(file_path), disable=disable):
            if line.startswith('#'):
                continue
            line = line.strip('\n')
            columns = [None if col == '.' else col for col in line.split('\t')]
            assert len(columns) == len(gff_record_fields)
            record = dict(zip(gff_record_fields, columns))
            record['start'] = int(record['start'])
            record['end'] = int(record['end'])
            if 'attributes' in record.keys():
                record['attributes'] = dict() if record['attributes'] == '.' else dict(
                    item.split('=') for item in record['attributes'].split(';'))
            yield record


def extract_sequence(assembly, chromosome, strand, start=None, end=None):
    """Get DNA sequence from a 2bit file"""
    if start is not None:
        start = f' -start={start}'
    else:
        start = ''
    if end is not None:
        end = f' -end={end}'
    else:
        end = ''
    result = subprocess.run(f'{django.conf.settings.DESIGN_TWO_BIT_TO_FA_PATH} -seq={chromosome}{start}{end} {django.conf.settings.DESIGN_ASSEMBLIES_FOLDER}/{assembly}.2bit stdout',
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, encoding='ascii')

    m = regex.search(r'>= seqSize \((\d+)\)', result.stderr)
    if m:
        result = subprocess.run(
            f'{django.conf.settings.DESIGN_TWO_BIT_TO_FA_PATH} -seq={chromosome}{start} -end={m.groups()[0]} {django.conf.settings.DESIGN_ASSEMBLIES_FOLDER}/{assembly}.2bit stdout',
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, encoding='ascii')
    output = result.stdout.split('\n')
    header, sequence = output[0], ''.join(output[1:])
    if strand == '-' and sequence:
        sequence = reverse_complement(sequence)
    return sequence, result.stderr


def extract_sequences(seqlist, assembly, outfile):
    """Get multiple DNA sequences from a 2bit file"""
    result = subprocess.run(f"{django.conf.settings.DESIGN_TWO_BIT_TO_FA_PATH} -seqList={seqlist} {django.conf.settings.DESIGN_ASSEMBLIES_FOLDER}/{assembly}.2bit {outfile}",
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, encoding='ascii')
    m = regex.search(r'twoBitReadSeqFrag in \S+ end \((\d+)\) >= seqSize', result.stderr)

    while m:
        with open(seqlist, 'r') as inf, open(f'{seqlist}-2', 'w') as outf:
            for line in inf:
                if not line.endswith(f'{m.groups()[0]}\n'):
                    outf.write(line)
        os.rename(f'{seqlist}-2', seqlist)
        result = subprocess.run(
            f"{django.conf.settings.DESIGN_TWO_BIT_TO_FA_PATH} -seqList={seqlist} {django.conf.settings.DESIGN_ASSEMBLIES_FOLDER}/{assembly}.2bit {outfile}",
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, encoding='ascii')
        m = regex.search(r'twoBitReadSeqFrag in \S+ end \((\d+)\) >= seqSize', result.stderr)

    return outfile


def run_bowtie(query, assembly):
    """Align a single sequence with bowtie."""
    cmd = f'{django.conf.settings.DESIGN_BOWTIE_PATH} -v 3 -a --best --sam-nohead -x {django.conf.settings.DESIGN_BOWTIE_GENOMES_FOLDER}/{assembly}/{assembly} --suppress 1,5,6,7 -c {query} -y --threads {django.conf.settings.DESIGN_BOWTIE_THREADS}'
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
    for line in iter(p.stdout.readline, ''):
        yield line
    p.stdout.close()
    #return subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)


def run_bowtie_multi(queries, assembly):
    """Align a list of sequences with bowtie."""
    cmd = f"{django.conf.settings.DESIGN_BOWTIE_PATH} -v 3 -a --best --sam-nohead -x {django.conf.settings.DESIGN_BOWTIE_GENOMES_FOLDER}/{assembly}/{assembly} --suppress 5,6,7 -r {queries} -y --threads {django.conf.settings.DESIGN_BOWTIE_THREADS}"
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
    for line in iter(p.stdout.readline, ''):
        yield line
    #return subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)


def run_bowtie_pairs(pair1, pair2, assembly, product_min_size, product_max_size, **options):
    """Pairwise alignment of two lists of sequences using bowtie."""
    cmd = f"{django.conf.settings.DESIGN_BOWTIE_PATH} -v 3 -k 50 --best --sam-nohead -x {django.conf.settings.DESIGN_BOWTIE_GENOMES_FOLDER}/{assembly}/{assembly} --suppress 2,6,7 -r -1 {pair1} -2 {pair2} -I {int(product_min_size/3)} -X {int(product_max_size*3)} -y --threads {django.conf.settings.DESIGN_BOWTIE_THREADS}"

    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)

    for line in iter(p.stdout.readline, ''):
        yield line

    #return subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)


def run_primer3(sequence, start, length, product_min_size, product_max_size, primer_min_length, primer_max_length, primer_opt_length, primer_min_tm, primer_max_tm, primer_opt_tm, **options):
    """Design primers using primer3."""
    primer3_template = f"""PRIMER_THERMODYNAMIC_PARAMETERS_PATH={django.conf.settings.DESIGN_PRIMER3_CONFIG_PATH}
SEQUENCE_TEMPLATE={sequence}
SEQUENCE_TARGET={start},{length}
PRIMER_PRODUCT_SIZE_RANGE={product_min_size}-{product_max_size}
PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE=3
PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE=3
PRIMER_SECONDARY_STRUCTURE_ALIGNMENT=1
PRIMER_NUM_RETURN=10
PRIMER_MAX_HAIRPIN_TH=47.00
PRIMER_INTERNAL_MAX_HAIRPIN_TH=47.00
PRIMER_MAX_END_STABILITY=9.0
PRIMER_EXPLAIN_FLAG=1
PRIMER_LIBERAL_BASE=1
PRIMER_FIRST_BASE_INDEX=1
PRIMER_MIN_SIZE={primer_min_length}
PRIMER_OPT_SIZE={primer_opt_length}
PRIMER_MAX_SIZE={primer_max_length}
PRIMER_MIN_TM={primer_min_tm}
PRIMER_OPT_TM={primer_opt_tm}
PRIMER_MAX_TM={primer_max_tm}
="""

    p = subprocess.Popen(django.conf.settings.DESIGN_PRIMER3_PATH, stdin=subprocess.PIPE, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         universal_newlines=True)
    for line in primer3_template.split('\n'):
        p.stdin.write(line + '\n')
    p.stdin.close()

    return p.stdout.read()


def check_primer_specificity(primers, assembly, write_folder, **options):
    """Align primer pairs to find potential off-targets"""

    mm_pattern = regex.compile(r'(?P<position>\d+):(?P<to>\w)>(?P<from>\w)')

    fwd_file = os.path.join(write_folder, 'fwd')
    rev_file = os.path.join(write_folder, 'rev')
    with open(fwd_file, 'w') as fwd, open(rev_file, 'w') as rev:

        for pair in primers:
            fwd.write(f"{pair['primers']['LEFT']['SEQUENCE']}\n")
            fwd.write(f"{pair['primers']['LEFT']['SEQUENCE']}\n")
            fwd.write(f"{pair['primers']['RIGHT']['SEQUENCE']}\n")
            rev.write(f"{pair['primers']['LEFT']['SEQUENCE']}\n")
            rev.write(f"{pair['primers']['RIGHT']['SEQUENCE']}\n")
            rev.write(f"{pair['primers']['RIGHT']['SEQUENCE']}\n")

    f = run_bowtie_pairs(fwd_file, rev_file, assembly, **options)
    for plus in f:
        minus = next(f)
        plus_primer, plus_chr, plus_start, plus_seq, plus_mm = plus.split('\t')
        plus_start = int(plus_start)
        plus_seq = list(plus_seq)
        for mm in regex.findall(mm_pattern, plus_mm):
            plus_seq[int(mm[0])] = mm[1].lower()

        plus_seq = ''.join(plus_seq)
        minus_primer, minus_chr, minus_start, minus_seq, minus_mm = minus.split('\t')

        minus_seq = list(minus_seq)
        for mm in regex.findall(mm_pattern, minus_mm):
            minus_seq[int(mm[0])] = mm[1].lower()
        minus_seq = ''.join(minus_seq)
        minus_start = int(minus_start)
        pair = int(plus_primer.split('/')[0]) // 3
        product = (plus_chr, plus_seq, minus_seq, plus_start,
                   minus_start + len(minus_seq),
                   minus_start + len(minus_seq) - plus_start)

        yield pair, product
