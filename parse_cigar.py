import re    
import sys
import argparse

# primary regex for parsing CIGAR strings:
end_clipping_regex = re.compile(r'[0-9]+S')
splicing_regex = re.compile(r'(?P<before_splice>([0-9]+[=DIMXSH]+)+)+(?P<splice_len>[0-9]+)N(?P<after_splice>.*)')
splice_len_regex = re.compile(r'(?P<remainder>.*[=XDIS])(?P<splice_len>[0-9]+)')
# checking for softclips:
softclip_start_regex = re.compile(r'^([0-9]+)S')
softclip_end_regex = re.compile(r'.*=([0-9]+)S')

def calculate_transcript_length(CIGAR:str)->int:
    """
    A method for decoding total length of a CIGAR string
    Methods:
    1. Strip all soft clipping from the CIGAR string
    2. Split on X and = and sum all free #'s
    3. Split on I and D and sum all free #'s
    4. Add the sum of all free #'s to the abs(deletion_cnt) to get accurate positioning on reference in jbrowse, insertions are ignored for binning purposes
    
    Method change:
    11/03/23 deletions are now included in the final length of the transcript to provide accurate positioning on 
    """
    len_sum = 0
    # removes soft clipping from both ends of the CIGAR string
    CIGAR_clips_rmd = re.sub(end_clipping_regex,'',CIGAR)
    # split on X and = to get the number of positions that are added to get mapping position
    add_vals = re.split('X|=',CIGAR_clips_rmd)
    #print(add_vals)
    insertion_del_list = []
    # sum everything thats not an insertion or deletion
    for val in add_vals:
        # controls for insertions and deletions
        if 'I' not in val and 'D' not in val and 'N' not in val and val != '':
            len_sum += int(val)
        else:
            insertion_del_list.append(val)
    #print("Len sum",len_sum,"\nInsertion_deletion_list:",insertion_del_list)
    del_cnt = 0
    for indel in insertion_del_list:
        #print(indel)
        # make sure its not empty or has soft clipping or other splicing in it
        if indel != '' and 'S' not in indel:
            indel_split = re.split('D|I',indel)
            indel_split_len = len(indel_split)
            #print(indel_split,indel_split_len)
            del_cnt = insertion_deletion_running_total(indel,del_cnt)
            if indel_split_len == 3:
                if str(indel_split[2]) != '':
                    len_sum += int(indel_split[2]) 
            elif indel_split_len == 2:
                if str(indel_split[1]) != '':
                    len_sum += int(indel_split[1])
            else:
                pass
        else:
            pass
    len_sum += abs(del_cnt)
    return len_sum

def insertion_deletion_running_total(CIGAR:str,del_cnt:int)->int:
    if 'I' in CIGAR:
        if re.split('I',CIGAR)[0] != '':
            #insert_del_cnt += int(re.split('I',CIGAR)[0])
            pass
    elif 'D' in CIGAR:
        if re.split('D',CIGAR)[0] != '':
            del_cnt -= int(re.split('D',CIGAR)[0])
    else:
        return del_cnt
    #print("deletion count",del_cnt)
    return del_cnt

def decode_CIGAR(CIGAR:str,read_start:int)->list[tuple[int,int]]:
    """
    A method designed to parse advanced CIGAR strings with splices and determine the total sequence length

    returns a list of tuples with the start and end positions of each exon on the reference

    Methods:
    
    CIGAR "1621S144=1D84=865N307=1X51N1068=1D998="  yields --> ['144=1D84=865', '307=1X51', '1068=1D998=']
    
    On reference meaning we ignore insertions and deletions because they are not present on the reference and we want to categorize these relative to our reference.
    
    By using up front filtering we should be able to remove alignments that are very weak...
    """
    # removes soft clipping from both ends of the CIGAR string
    CIGAR_clips_rmd = re.sub(end_clipping_regex,'',CIGAR)
    # split on N's and get multiple fragments
    splice_frags = re.split('N',CIGAR_clips_rmd)
    # splitting on N's here: "1621S144=1D84=865N307=1X51N1068=1D998="  yields --> ['144=1D84=865', '307=1X51', '1068=1D998=']
    splice_coords = []
    exon_cnt = 1
    exon_start = 0
    splice_frag_len = len(splice_frags)
    for fragment in splice_frags:
        # the last fragment ex:'1068=1D998=' will not have an N so we simply decode its length directly...
        #print("[STDOUT]: fragment:",fragment,"\n[STDOUT]: frag length:",splice_frag_len,"\n[STDOUT]: total exon cnt:",exon_cnt)
        if splice_frag_len == exon_cnt:
            remainder_len = calculate_transcript_length(fragment)
            # handles errors with single exon transcripts
            if splice_frag_len == 1:
                next_exon_start = read_start
            exon_start = next_exon_start
            exon_end = exon_start + int(remainder_len)
            splice_coords.append((exon_start,exon_end-1))
            break
        else:
            match = splice_len_regex.match(fragment)
            if match:
                splice_len = match.group('splice_len')
                remainder = match.group('remainder')
                remainder_len = calculate_transcript_length(remainder)
                #print(remainder,remainder_len)
                if exon_cnt == 1:
                    exon_start = read_start
                    exon_end = exon_start + int(remainder_len)
                else:
                    exon_start = next_exon_start
                    exon_end = exon_start + int(remainder_len)
                # assign next exon start position based on previous splice
                next_exon_start = exon_end + int(splice_len)
                # store your found splice coords:
                splice_coords.append((exon_start,exon_end-1))
        exon_cnt += 1
    return splice_coords

def get_softclips(cigar:str,dna_seq:str):
    """
    A method for parsing the softclips from a CIGAR string and returning only the softclipped sequence
    """
    if 'S' in cigar:
        print("CIGAR:",cigar)
        # check the start
        match = softclip_start_regex.match(cigar)
        if match:
            softclip_start_len = int(match.group(1))
        else:
            softclip_start_len = 0
            
        # check the end
        match = softclip_end_regex.match(cigar)
        if match:
            softclip_end_len = int(match.group(1))
        else:
            softclip_end_len = 0
    print("Softclip start:",softclip_start_len,dna_seq[:softclip_start_len],"\nSoftclip end:",softclip_end_len,dna_seq[-softclip_end_len:],"\n")
    return(dna_seq[:softclip_start_len],dna_seq[-softclip_end_len:])

def write_to_fasta(seq:str,fasta_out:object,read_name:str):
    fasta_out.write(">"+read_name+"\n")
    fasta_out.write(seq+"\n")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="A method for parsing CIGAR strings with splices and determining the total sequence length")
    parser.add_argument('--start_pos',type=int,default=-1,help="The starting position of the read on the reference")
    parser.add_argument('--cigar',type=str,default='',help="The CIGAR string to be parsed")
    parser.add_argument('--sam',type=str,default='',help="The SAM file to be parsed")
    args = parser.parse_args()
    clip_len_threshold = 50
    # options for cigar parsing
    if args.start_pos != -1 and sys.argv[2] != '':
        read_start_pos=int(sys.argv[1])
        cigar_input=str(sys.argv[2])           #'13=1X1=1X683=1223N981=1X291=1X1=1X169=1S'
        if len(re.findall('N',cigar_input)) > 1:
            print(decode_CIGAR(cigar_input,args.start_pos))
    elif args.sam != '':
        with open(args.sam,'r') as sam:
            output_fasta = open(args.sam[:-4]+"_soft_clips.fa",'w+')
            for line in sam:
                if line[0] != '@':
                    line = line.strip().split('\t')
                    read_name = line[0]
                    genome_pos = line[3]
                    cigar_str = line[5]
                    seq = line[9]
                    # retrieve softclips from either end of the sequence if they exist:
                    softclip_start_seq,softclip_end_seq = get_softclips(cigar_str,seq)
                    # remove long poly-A runs:
                    if len(softclip_start_seq) > 0:
                        softclip_start_seq = softclip_start_seq.replace('A'*10,'')
                    if len(softclip_end_seq) > 0:
                        softclip_end_seq = softclip_end_seq.replace('A'*10,'')
                    # write your outputs to fasta:
                    if len(softclip_start_seq) > clip_len_threshold:
                        write_to_fasta(softclip_start_seq,output_fasta,read_name+"_clip_start")
                    if len(softclip_end_seq) > clip_len_threshold:
                        write_to_fasta(softclip_end_seq,output_fasta,read_name+"_clip_end")

        

# some multiple splice test cases:
# python ~/progs/Utilities/HBV_seq_analysis/parse_cigar.py 188 '1S19=1D20=2X70=1D32=1D20=1D29=1D16=1I9=1I49=1I11=1I21=262N14=1D16=1I47=1I14=1D5=1223N19=1D14=1D9=1D49=1D21=1I5=1I32=1D43=1D11=2I8=1D19=1X12=1D19=1X27=1X10=1I2=1X1=1X4=1I15=1I15=1D11=1I8=1D28=1I2=1I23=2I3=1D6=1I5=1I11=1I7=1I27=1I4=1I2=1I5=1I10=2I6=3I3=3I4=1D9=2I14=1I6=1I11=1I29=1I14=1I32=1I22=1I3=1I2=2I2=1I2=1X20=1I15=1I16=1I9=1I9=1I6=1I1=1I13=1I42=1I23=1X6=1I8=1I5=1I6=1I3=1I1X20=1I30=1I4=1I3=1I2=1I15=1I6=1I3=1I10=1I10=2I8=1X21=1X22=1I2=6I14=1I1=1I14=1I2=1D15=1I7=1I3=1X10=1I3=1I2=1I6=2I5=1I3=1I19=1X1=3I11=1I9=1I17=1I18=1I41=1I23=1I22=1I17=1I1=1X1=1X24=1X4=1D20=1I47=1I16=2I2=8I3=1I4=1X5=1I18=1I7=1I4=1I2=1I9=4S'
# [(188, 490), (753, 850), (2074, 3517)]

# python ~/progs/Utilities/HBV_seq_analysis/parse_cigar.py 282  '1S72=1D1341=1X111=1D675=1D273=1X264=1X31=1X291=1X1=1X23=47N102=1S'
# [(282, 3372), (3420, 3521)]

# python ~/progs/Utilities/HBV_seq_analysis/parse_cigar.py 204 '1S287=262N98=1223N981=1X1=1X21=1X267=1X1=1X169=2S'
# [(204, 490), (753, 850), (2074, 3518)]

# advanced single exon cases:
# python ~/progs/Utilities/HBV_seq_analysis/parse_cigar.py 192 '1S7=1D21=1I131=1I31=1I4=1I16=2I39=1I37=1I36=1I23=1I110=1I53=1I24=1I41=1I13=1D34=1D23=1I14=1I25=1D9=1I11=1I45=1I68=1D5=1D16=1I39=1I169=1I25=1I20=1I10=1I18=1X16=1I6=1I4=1I11=1I5=1X26=1I15=2X40=1D4=1D8=1I5=1I5=1I10=1I1=1I25=2X17=1I110=1D31=1I23=1I6=1I3=1D1X25=1I12=1I21=1I15=1D11=1I12=1I19=2X33=1I11=1I16=1I16=1I36=1D12=1I7=1I8=1I21=1I57=1D27=1D35=1D19=1I46=1X31=1I11=1I32=1I23=1D32=2X32=1I36=1I17=1I33=1I7=1I1X6=1I127=1I94=1I22=1I79=1I3=1I9=1X43=2I9=1I2=1I29=1I5=2I8=1I19=1I8=1I4=1I12=1I12=1I18=1I10=1I59=2I8=1I5=1I27=1I2=1I5=1I7=1X9=1I29=1D15=2I59=1X18=1I32=1I127=1X1=1X14=1I15=2I8=1I31=1I4=1I3=1I28=1I6=1X6=1I17=1I5=1I7=1I2=1I7=1I12=5S'
# [(192, 3515)]

# python ~/progs/Utilities/HBV_seq_analysis/parse_cigar.py 192 '2S17=1D5=1I9=1I31=1I34=1I7=1D49=1I8=1I119=1I15=1D20=1I163=2I18=1I43=1I34=1D15=2X70=1I144=1X43=1D93=1I119=2I21=1I6=1D9=1D16=1D21=1D5=1I134=1D75=1I40=1I73=1I12=1I23=1X13=1X58=1I39=1D48=1I25=1I5=1I65=1I40=1D29=1I3=1X30=1I13=1I8=1I26=1D38=1I40=1D39=1D41=1D27=1D16=1I20=1I56=1I16=1I6=3I2=1X1=1I5=1I10=1I27=1I14=1I8=1I4=1I5=1D7=1I7=1I12=1I3=1I12=2I4=1D16=1I13=2I6=1I6=1I33=2I6=1I38=1I9=1I3=1D3=1I6=1I3=1I5=1I18=1I1X8=1I6=1I6=1I5=2X9=1I2=1D12=1I22=1I17=1I33=1I30=1I5=1I6=1I3=1I10=1I3=1I10=1I18=1I14=1I26=1I18=1I19=1I2=2I24=1I12=1I33=1I1=1X13=2I10=1I6=1I26=1I5=1I7=1I9=1D4=1X5=1I12=2I15=1I7=1I14=2I1=2I8=1I11=2I29=1I3=1I8=1I13=1X24=1I9=2I3=1I5=1I5=1I9=1I26=2I1=1X1=1X2=1I9=1I10=1D7=1D26=2I16=1D12=1I1X9=1I5=1I3=1I13=1I9=1I3=1I29=1I8=1I9=' 
# [(192, 3523)]