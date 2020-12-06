### Our new compression tools for fastq/fasta-files:

#### At the moment, compressing of fastq-files only  is available

__USAGE for compressing__:
```
python3 compress.py input_read.fastq
```

After this stage 4 compressed files will be available in the `script_out` directory.


__USAGE for compressing__:
```
python3 decompress.py output_read.fastq
```

Requirments:
`BioPython`, `numpy`, `brotli`
