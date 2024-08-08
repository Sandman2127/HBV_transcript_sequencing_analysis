### Introduction:
The analyze_read_alignments.py program bins aligned FLNC reads by parsing each reads alignment start and CIGAR string from SAM format file(s) to determine alignment locations and any breaks. Using the coordinates discovered the program searches for an appropriate transcript bin. To control for differences in sequencing library depth, the raw read counts are normalized via the reads per million (RPM) normalization technique. Finally, plots and tab separated value (.tsv) summary files of the transcript binning process are output.

If tests or use of the parse_cigar.py methods are desired a separate main can be directly called.

### Required Software Dependencies:
```
python = 3.11
pandas = 1.5.3
bokeh = 3.4.1
seaborn = 0.12.2
matplotlib = 3.7.1
```
### Help Menus:
```
python analyze_read_alignments.py --help
usage: A program for automated charecterization of a list of sam files based on an input transcript_bed file
       [-h] --samfiles SAMFILES [SAMFILES ...] --transcript_bed TRANSCRIPT_BED
       [--min_splice_len MIN_SPLICE_LEN] [--debug] [--min_qual MIN_QUAL]
       [--sample_weights SAMPLE_WEIGHTS] [--combined_bam COMBINED_BAM]

options:
  -h, --help            show this help message and exit
  --samfiles SAMFILES [SAMFILES ...]
                        list of sam files
  --transcript_bed TRANSCRIPT_BED
                        bed format of transcriptome, see script for details of format
  --min_splice_len MIN_SPLICE_LEN
                        minimum length of a splice to be considered, default: 300
  --debug               default: False
  --min_qual MIN_QUAL   minimum quality score to consider a read, default: 60
  --sample_weights SAMPLE_WEIGHTS
                        add sample weights in the order of input files calculated by calculation
                        script as csv formated string example: 1,1.4,0.79
  --combined_bam COMBINED_BAM
                        combined bam file of all samples


python parse_cigar.py --help
usage: parse_cigar.py [-h] [--start_pos START_POS] [--cigar CIGAR] [--sam SAM]

A method for parsing CIGAR strings with splices and determining the total sequence length

options:
  -h, --help            show this help message and exit
  --start_pos START_POS
                        The starting position of the read on the reference
  --cigar CIGAR         The CIGAR string to be parsed
  --sam SAM             The SAM file to be parsed

```

### Standard Usage:
```
python3 analyze_read_alignments.py --transcript_bed ./HBV_transcript_data/20240715_Human_HBV_transcripts_w_splices_2700-2200.bed --sample_weights 0.2120,0.2235,0.4103 --samfiles input.sam input2.sam input3.sam --combined_bam bams.merged
```

