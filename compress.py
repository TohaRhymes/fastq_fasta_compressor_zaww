import _brotli
import numpy as np
from Bio import SeqIO
import os.path
import sys

MIN_PHRED = 33
START_PHRED = ord('@')
MIN_QUAL = 20

# The compression mode.
MODE_GENERIC = _brotli.MODE_GENERIC
MODE_TEXT = _brotli.MODE_TEXT
MODE_FONT = _brotli.MODE_FONT

# The Compressor object.
Compressor = _brotli.Compressor

# The Decompressor object.
Decompressor = _brotli.Decompressor


# do_batch = False
#
# BATCH = ''
# if do_batch:
#     BATCH = '\n'


def dir_creation(path):
    try:
        os.mkdir(path)
    except OSError:
        pass


def loss_t_gaps(q_scores):
    tgaps = []
    tgaps.append(q_scores[0] - MIN_PHRED)
    for i in range(len(q_scores) - 1):
        tgaps.append(q_scores[i + 1] - q_scores[i])
    tgaps = [(2 * x) if x > 0 else (-2 * x + 1) for x in tgaps]
    tgaps = [2 ** int(np.log2(x)) if (x % 2 == 0 or x == 1) else 2 ** int(np.log2(x)) + 1 for x in tgaps]
    return tgaps


def qvals_to_tgaps(q_scores):
    # q_scores = [ord(q.compress.py('ascii')) for q in q_scores]
    loss_tgaps = loss_t_gaps(q_scores)
    # loss_tgaps = [2 if (x > 127 and x % 2 == 0) else 1 if (x > 127 and x % 2 == 1) else x for x in loss_tgaps]
    ans = bytes(bytearray(loss_tgaps)).decode('ascii')
    return ans


def from_fasta_fastq(filename, names_file, seq_file, dels_file, qual_file=None, fastq=True):
    """
    Function reads fastq-file and write 4

    :param filename: input fastq-file.
    :param names_file: output file with names of reads.
    :param seq_file: output file with merged sequences (one string).
    :param dels_file: output file with length of reads
    :param qual_file: output file with merged qualities (one string).
    :param fastq: boolean, is input file FastQ (true), or Fasta (false)
    :return: Nothing:)
    """
    if fastq:
        o_qual_file = open(qual_file, 'wb')

    flag = "fastq" if fastq else "fasta"

    with open(names_file, 'w') as o_names_file, \
            open(seq_file, 'w') as o_seq_file, \
            open(dels_file, 'w') as o_dels_file:
        for record in SeqIO.parse(filename, flag):
            if not fastq:
                o_names_file.write('%s\n' % record.name)
                # o_names_file.write('%s\n' % bwt(record.name))
                o_seq_file.write('%s' % seq)
                # o_seq_file.write('%s' % bwt(record.seq))
                o_dels_file.write('%d\n' % len(record))
            if fastq:
                qual = ''
                counter = 0
                start = True
                counter_start = 0
                for q in record.letter_annotations['phred_quality']:
                    qual += chr(q + MIN_PHRED)
                    if q + MIN_PHRED - START_PHRED > MIN_QUAL:
                        if start:
                            counter_start = counter
                        counter = 0
                        start = False
                    else:
                        counter += 1
                seq = record[counter_start:len(record) - counter].seq
                qual = record.letter_annotations['phred_quality'][counter_start:len(record) - counter]
                if len(seq) > 0:
                    o_names_file.write('%s\n' % record.name)
                    # o_names_file.write('%s\n' % bwt(record.name))
                    o_dels_file.write('%d\n' % len(seq))
                    o_seq_file.write('%s' % seq)
                    # o_seq_file.write('%s' % bwt(seq))
                    o_qual_file.write(qvals_to_tgaps(qual).encode('ascii'))
                # print(qual, seq)
    if fastq:
        o_qual_file.close()


# Compress a byte string.
def compress(string, mode=MODE_TEXT, quality=11, lgwin=22, lgblock=0):
    """Compress a byte string.
    Args:
      string (bytes): The input data.
      mode (int, optional): The compression mode can be MODE_GENERIC (default),
        MODE_TEXT (for UTF-8 format text input) or MODE_FONT (for WOFF 2.0).
      quality (int, optional): Controls the compression-speed vs compression-
        density tradeoff. The higher the quality, the slower the compression.
        Range is 0 to 11. Defaults to 11.
      lgwin (int, optional): Base 2 logarithm of the sliding window size. Range
        is 10 to 24. Defaults to 22.
      lgblock (int, optional): Base 2 logarithm of the maximum input block size.
        Range is 16 to 24. If set to 0, the value will be set based on the
        quality. Defaults to 0.
    Returns:
      The compressed byte string.
    Raises:
      brotli.error: If arguments are invalid, or compressor fails.
    """
    compressor = Compressor(mode=mode, quality=quality, lgwin=lgwin,
                            lgblock=lgblock)
    return compressor.process(string) + compressor.finish()

try:
    # Decompress a compressed byte string.
    decompress = _brotli.decompress

    # Raised if compression or decompression fails.
    error = _brotli.error

    # from where to read
    dir_in = './'

    input_file = 'example/read.fastq'

    # 4 output files from unsqueze
    dir_out = './script_out'
    dir_creation(dir_out)


    out_names = 'names_file'
    out_seq = 'seq_file'
    out_qual = 'qual_file'
    out_dels = 'dels_file'

    out_names = os.path.join(dir_out, out_names)
    out_seq = os.path.join(dir_out, out_seq)
    out_qual = os.path.join(dir_out, out_qual)
    out_dels = os.path.join(dir_out, out_dels)

    from_fasta_fastq(os.path.join(dir_in, input_file),
                     out_names,
                     out_seq,
                     out_dels,
                     out_qual)

    with open(out_names, 'r') as f1, open(out_seq, 'r') as f2, open(out_qual, 'r') as f3, open(out_dels, 'r') as f4:
        for name_file, file in zip([out_names, out_seq, out_qual, out_dels], [f1, f2, f3, f4]):
            sequance = file.read().encode('ascii')
            compress_sequance = compress(sequance)
            with open(name_file + '.zaww', 'wb') as file_zaww:
                file_zaww.write(compress_sequance)
    for file in [out_names,
                 out_seq,
                 out_dels,
                 out_qual]:
        os.remove(file)
except IndexError:
        print('Usage:\n python3 decompress.py <intput.fastq>')
