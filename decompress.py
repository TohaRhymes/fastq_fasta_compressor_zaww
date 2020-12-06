import os.path
import sys

MIN_PHRED = 33
START_PHRED = ord('@')
MIN_QUAL = 20

import _brotli

# The compression mode.
MODE_GENERIC = _brotli.MODE_GENERIC
MODE_TEXT = _brotli.MODE_TEXT
MODE_FONT = _brotli.MODE_FONT

# The Compressor object.
Compressor = _brotli.Compressor

# The Decompressor object.
Decompressor = _brotli.Decompressor


# Deompress a byte string.
def decompress_zaww():
    with open(out_names + '.zaww', 'rb') as f1, \
            open(out_seq + '.zaww', 'rb') as f2, \
            open(out_qual + '.zaww', 'rb') as f3, \
            open(out_dels + '.zaww', 'rb') as f4:
        for name_file, file in zip([out_names, out_seq, out_qual, out_dels], [f1, f2, f3, f4]):
            sequence = file.read()
            compress_sequence = decompress(sequence)
            with open(name_file, 'wb') as file_decompress:
                file_decompress.write(compress_sequence)


def tgaps_to_qvals(tgaps):
    tgaps = [ord(q.encode('ascii')) for q in tgaps]
    q_vals = [x // 2 if x % 2 == 0 else -(x - 1) // 2 for x in tgaps]
    q_vals[0] += 2 * MIN_PHRED
    for i in range(len(tgaps) - 1):
        q_vals[i + 1] += q_vals[i]
    q_vals = bytes(bytearray(q_vals)).decode('ascii')
    return q_vals


def to_fasta_fastq(filename, names_file, seq_file, dels_file, qual_file=None, fastq=True):
    """
    Function reads 4 special files and makes one fasta/fastq-file.

    :param filename: fastq/fasta-file for output.
    :param names_file: input file with names of reads.
    :param seq_file: input file with merged sequences (one string).
    :param dels_file: input file with length of reads.
    :param qual_file: input file with merged qualities (one string).
    :param fastq: boolean, is output file FastQ (true), or Fasta (false)
    :return: Nothing:)
    """
    if fastq:
        o_qual_file = open(qual_file, 'r')

    flag = "@" if fastq else ">"

    with open(names_file, 'r') as o_names_file, \
            open(seq_file, 'r') as o_seq_file, \
            open(dels_file, 'r') as o_dels_file, \
            open(filename, 'w') as out_file:

        try:
            while 1:
                length = int(o_dels_file.readline())
                name = o_names_file.readline()[:-1]
                # name = inverse_bwt(o_names_file.readline()[:-1])
                out_file.write('%s%s\n' % (flag, name))
                out_file.write('%s\n' % o_seq_file.read(length))
                # out_file.write('%s\n' % inverse_bwt(o_seq_file.read(length+1)))
                if fastq:
                    out_file.write('+\n%s\n' % tgaps_to_qvals(o_qual_file.read(length)))
        except ValueError:
            pass

try:
    # Decompress a compressed byte string.
    decompress = _brotli.decompress

    # Raised if compression or decompression fails.
    error = _brotli.error

    # 4 output files from unsqueze
    dir_out = './script_out'


    dir_fastq = './'

    out_names = os.path.join(dir_out, 'names_file')
    out_seq = os.path.join(dir_out, 'seq_file')
    out_qual = os.path.join(dir_out, 'qual_file')
    out_dels = os.path.join(dir_out, 'dels_file')

    decompress_zaww()

    # TO ONE FILE
    try:
        out_file = sys.argv[1]

        to_fasta_fastq(os.path.join(dir_fastq, out_file),
                       out_names,
                       out_seq,
                       out_dels,
                       out_qual
                       )

        for file in [out_names,
                     out_seq,
                     out_dels,
                     out_qual]:
            os.remove(file)
    except IndexError:
        print('Usage:\n python3 decompress.py <output.fastq>')
except FileNotFoundError:
    print('There is no package script_out with 4 files')