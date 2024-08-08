#!/bin/bash

PY3="/opt/anaconda/envs/py3/bin/python"

ANALYZE_HBV_ALIGNMENTS="$HOME/progs/Utilities/HBV_seq_analysis/analyze_read_alignments.py"
TRANSCRIPT_BED="$HOME/progs/Utilities/HBV_seq_analysis/HBV_transcript_data/HBV_transcript.bed"
TRANSCRIPT_ORIG_BED="$HOME/progs/Utilities/HBV_seq_analysis/HBV_transcript_data/HBV_transcript_original.bed"
#TRANSCRIPT_SPLICED_BED="$HOME/progs/Utilities/HBV_seq_analysis/HBV_transcript_data/HBV_transcripts_w_splices.bed"
TRANSCRIPT_SPLICED_BED="$HOME/progs/Utilities/HBV_seq_analysis/HBV_transcript_data/HBV_transcripts_w_splices_110828.bed"
SAMP_WEIGHTS="0.2120,0.2536,0.2490,0.2434,0.1966,0.2190,0.2117,0.2183,0.3094,0.2233,0.2885,0.2508,0.2514"

SAMPLES="89A008_HC.ccs_polish-flnc.fa.algn.srt.sam
89A008_Day_141.ccs_polish-flnc.fa.algn.srt.sam
89A008_Day_379.ccs_polish-flnc.fa.algn.srt.sam
A2A004_HC.ccs_polish-flnc.fa.algn.srt.sam
A2A004_d85.ccs_polish-flnc.fa.algn.srt.sam
A2A004_d323.ccs_polish-flnc.fa.algn.srt.sam
A2A004_Day_351.ccs_polish-flnc.fa.algn.srt.sam
A3A006_HC.ccs_polish-flnc.fa.algn.srt.sam
A3A006_d141.ccs_polish-flnc.fa.algn.srt.sam
A3A006_Day_351.ccs_polish-flnc.fa.algn.srt.sam
88A010_57.ccs_polish-flnc.fa.algn.srt.sam
95A010_d57.ccs_polish-flnc.fa.algn.srt.sam
A4A014_HC.ccs_polish-flnc.fa.algn.srt.sam"

BAMSMERGED="merged.bam"

#$PY3 $ANALYZE_HBV_ALIGNMENTS --transcript_bed $TRANSCRIPT_BED --min_edge_buffer 50 --sample_weights $SAMP_WEIGHTS --samfiles $SAMPLES 

$PY3 $ANALYZE_HBV_ALIGNMENTS --transcript_bed ${TRANSCRIPT_SPLICED_BED} --sample_weights $SAMP_WEIGHTS --samfiles $SAMPLES --combined_bam $BAMSMERGED --debug
