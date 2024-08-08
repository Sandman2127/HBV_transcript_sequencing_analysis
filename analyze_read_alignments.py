#!/opt/anaconda/envs/py3/bin/python3
import argparse
import pysam
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import random
from itertools import chain
from parse_cigar import decode_CIGAR
from bokeh_replot_analysis import replot_analysis
sns.set_theme(style="ticks")

parser = argparse.ArgumentParser(prog='A program for automated charecterization of a list of sam files based on an input transcript_bed file')
parser.add_argument("--samfiles",nargs='+',required=True,help="list of sam files")
parser.add_argument("--transcript_bed",type=str,required=True,help="bed format of transcriptome, see script for details of format")
parser.add_argument("--min_splice_len",type=int,default=300,help="minimum length of a splice to be considered, default: 300")
parser.add_argument("--debug",action='store_true',default=False,help="default: False")
parser.add_argument("--min_qual",type=int,default=60,help="minimum quality score to consider a read, default: 60")
parser.add_argument("--sample_weights",type=str,default="0",help="add sample weights in the order of input files calculated by calculation script as csv formated string example: 1,1.4,0.79")
parser.add_argument("--combined_bam",type=str,default="",help="combined bam file of all samples")

args = parser.parse_args()
debug = args.debug
min_qual = args.min_qual
min_splice_len = args.min_splice_len
samfiles = args.samfiles
use_custom_palatte = True
is_human_sample = True

class sample_reads_classified():
    def __init__(self,transcript_dict:dict[str:list]):
        # read alignment tracker
        col_names = ['chr','read_name','read_start','read_end','sequence']
        # instantiate a transcript dictionary to count instances in each sample:
        self.read_classification_dict = {}
        for transcript in transcript_dict.keys():
            self.read_classification_dict[transcript] = pd.DataFrame(columns=col_names)
        # append other after all others have been added
        self.read_classification_dict['Other'] = pd.DataFrame(columns=col_names)

    def categorize_hbv_read(self,aligned_chr:str,read_name:str,read_start:int,seq:str,cigar:str,transcript_dict:dict[str:list]):
        aligned_transcript_data = get_transcript_alignment(aligned_chr,read_name,read_start,cigar,transcript_dict)
        aligned_transcript,read_end = aligned_transcript_data[0],aligned_transcript_data[1]
        append_df = pd.DataFrame({'chr':aligned_chr,'read_name':read_name,'read_start':read_start,'read_end':read_end,'sequence':seq},index=[0])
        # auto-add classified read to unique database
        self.read_classification_dict[aligned_transcript] = pd.concat((self.read_classification_dict[aligned_transcript],append_df),ignore_index=True)

def _print_alignment_type(read_name,aligned_transcript,read_start,read_end,tss_upstream,tss_downstream,tes_upstream,tes_downstream):
    print("[STDOUT]: read {} aligns to {} region in start: {} > {} > {} and end: {} < {} < {}".format(read_name,aligned_transcript,tss_upstream,read_start,tss_downstream,tes_upstream,read_end,tes_downstream))

def _print_alignment_type_other(read_name,aligned_transcript,read_start,read_end):
    print("[STDOUT]: read {} aligns to {} region in start {} and end {}".format(read_name,aligned_transcript,read_start,read_end))

def is_special_2x_case(transcript_data:list[tuple[int,int]])->bool:
    """
    Method to check if a transcript is a special case of 2x splicing

    This allows us to detect HBV_2x transcrpts that make full HBV S antigen with any amount of gaps after the primary transcript

    special case: 1-620 start - orf - 1321 and any amount of gaps up to intervals 2406-2700,2701-3600,3601-5586,5587-5765
    """
    transcript_start = transcript_data[0][0]
    transcript_start_end = transcript_data[0][1]
    # occasionally transcript_start_end == transcript_end as theres no splicing
    transcript_end = transcript_data[-1][1]
    # intervals for special case, I hate hardcoding variables but this is a special case for HBV
    if is_human_sample:
        HBsAg_2X_intervals = [(2340,2716),(2717,3664),(3665,5651),(5652,5931)]
        HBV_variable_region_start = 656
        HBV_variable_region_end = 1352
    else:
        HBV_variable_region_start = 625
        HBV_variable_region_end = 1321
        HBsAg_2X_intervals = [(2406,2700),(2701,3600),(3601,5586),(5587,5765)]
    # check intervals
    for interval in HBsAg_2X_intervals:
        interval_start = interval[0]
        interval_end = interval[1]
        if transcript_start <= HBV_variable_region_start and transcript_start_end >= HBV_variable_region_end and transcript_end >= interval_start and transcript_end <= interval_end:
            return True
    return False

def get_transcript_alignment(aligned_chr:str,read_name:str,read_start:int,cigar:str,transcript_dict:dict[str:list])->tuple[str,int]:
    """
    Method for binning reads into transcript categories based on transcriptome dict:

    transcriptome dict provides intervals:
    
    genome coordinates:                0------1000------2000------3000------4000------5000------6000
    continuous interval coordinates:      {start} ------------ {end} 
    spliced interval coordinates:        {spl1_start}---{spl1_end}---{spl2_start}--{spl2_end}

    anything > 2 intervals w/ > min_splice_len is considered "Other" as there are very very few and its very difficult for a basic overlap program to account for it.
    """
    # 1: get actual transcript position data:
    transcript_position_data = decode_CIGAR(cigar,read_start)
    transcript_end = transcript_position_data[-1][1]
    # TODO: special case: allows transcripts with 1-620 start - orf - 1321 and any amount of gaps up to intervals 2406-2700,2701-3600,3601-5586,5587-5765
    if is_special_2x_case(transcript_position_data):
        input_transcript_is_spliced = False
        # below transcript_position_data[-1][1] negates the need to know how many splices there are for the special case
        ex1_start_pos,ex1_end_pos = transcript_position_data[0][0],transcript_position_data[-1][1]
    elif len(transcript_position_data) > 2:
        return ("Other",transcript_end)
    elif len(transcript_position_data) == 2 and (transcript_position_data[1][0] - transcript_position_data[0][1] > min_splice_len):
        # example return: [(170, 869), (1885, 3536)], account for splices above min_splice_len
        input_transcript_is_spliced = True
        ex1_start_pos,ex1_end_pos,ex2_start_pos,ex2_end_pos = transcript_position_data[0][0],transcript_position_data[0][1],transcript_position_data[1][0],transcript_position_data[1][1]
    elif len(transcript_position_data) == 2 and (transcript_position_data[1][0] - transcript_position_data[0][1] < min_splice_len):
        # option ignores short splices
        input_transcript_is_spliced = False
        ex1_start_pos,ex1_end_pos = transcript_position_data[0][0],transcript_position_data[1][1]
    else:
        # example return: [(170, 869)], single transcript
        input_transcript_is_spliced = False
        ex1_start_pos,ex1_end_pos = transcript_position_data[0][0],transcript_position_data[0][1]
        if debug:
            print("found a normal one",ex1_start_pos,ex1_end_pos)

    # 2. go through transcript dict and check if the aligns to each 
    for key,value in transcript_dict.items():
        chrom = value[0]
        # used for non-spliced reads
        transcript_start_upstream = value[1]
        transcript_start_downstream = value[2]
        transcript_end_upstream = value[3]
        transcript_end_downstream = value[4]
        # used for splicing only to simplify naming
        splice1_start_start = transcript_start_upstream
        splice1_start_end = transcript_start_downstream
        splice1_end_start = transcript_end_upstream
        splice1_end_end = transcript_end_downstream
        splice2_start_start = value[5]
        splice2_start_end = value[6]
        splice2_end_start = value[7]
        splice2_end_end = value[8]
        # check if transcript is spliced
        if input_transcript_is_spliced:
            # only evaluate spliced transcripts:
            if splice2_start_start != 0 and splice2_end_start != 0:
                # <start> <ex1> <splice1> <---splice---> <"end"> <ex2> <splice2_end>
                if ex1_start_pos >= splice1_start_start and ex1_start_pos <= splice1_start_end and ex1_end_pos >= splice1_end_start and ex1_end_pos <= splice1_end_end and ex2_start_pos >= splice2_start_start and ex2_start_pos <= splice2_start_end and ex2_end_pos >= splice2_end_start and ex2_end_pos <= splice2_end_end:
                    return (key,ex2_end_pos)
        else:
            # protect from binning non-spliced reads into in spliced transcripts:
            if splice2_start_start != 0 or splice2_end_start != 0:
                pass
            else:
                try:
                    if aligned_chr == chrom and read_start >= transcript_start_upstream and read_start <= transcript_start_downstream and ex1_end_pos >= transcript_end_upstream and ex1_end_pos <= transcript_end_downstream:
                        return (key,ex1_end_pos)
                except UnboundLocalError:
                    print("[STDERR]: UnboundLocalError in get_transcript_alignment")
                    print("read_name:",aligned_chr,read_name,read_start,transcript_start_upstream,transcript_start_downstream,transcript_end_upstream,transcript_end_downstream,transcript_position_data)
                    exit()
    # finally if not found...
    if debug:
        _print_alignment_type_other(read_name,"Other",read_start,transcript_end)
    return ("Other",transcript_end)

def categorize_reads(sam:str,transcript_dict:dict)->list:
    sample_data = sample_reads_classified(transcript_dict)
    sam_data = filter_reads(sam)
    for sam_key,sam_vals in sam_data.items():
        sam_key_split = sam_key.split('^')
        readName,alignment_start,cigar = str(sam_key_split[0]),int(sam_key_split[1]),str(sam_key_split[2])
        refName,seq = str(sam_vals[0]),str(sam_vals[1])
        sample_data.categorize_hbv_read(refName,readName,alignment_start,seq,cigar,transcript_dict)
    return sample_data

def filter_reads(sam:str)->dict:
    print("[STDOUT]: filtering input sam file: {}".format(sam))
    samp_read_data = {}
    with open(sam,'r') as f:
        while True:
            line = f.readline()
            if line:
                if line.startswith('@'):
                    continue
                else:
                    fields = line.strip().split()
                    readName = str(fields[0])
                    flag = int(fields[1])
                    refName = str(fields[2])
                    alignment_start = int(fields[3])
                    qual = int(fields[4])
                    cigar = str(fields[5])
                    len_on_ref = int(fields[8])
                    seq = str(fields[9])
                    if qual >= min_qual and flag == 0:
                        samp_read_data[readName + '^' + str(alignment_start) + '^' + cigar] = (refName,seq,len_on_ref)
            else:
                break
    print("[STDOUT]: {} reads from sam file {} passed min_qual: {} and flag: 0".format(len(samp_read_data.keys()),sam,min_qual))
    # check for duplicate alignments that both had flag:0 and qual >= min_qual
    # take the longest alignment
    print("[STDOUT]: checking for duplicate alignments and taking the longest read")
    rname_keys = {}
    for k in samp_read_data.keys():
        rname = k.split('^')[0]
        if rname not in rname_keys:
            rname_keys[rname] = [k]
        else:
            rname_keys[rname].append(k)
    
    for k,v in rname_keys.items():
        if len(v) > 1:
            # get the longest read
            longest_read = ''
            longest_read_len = 0
            for read in v:
                read_len = samp_read_data[read][2]
                if read_len > longest_read_len:
                    longest_read = read
                    longest_read_len = read_len
                else:
                    print("[STDOUT]: removing duplicate read:{} as longest_read len:{} and this reads len:{} ".format(read,longest_read_len,read_len))
            # remove all other reads
            for read in v:
                if read != longest_read:
                    print("[STDOUT]: removing duplicate read:{}".format(read))
                    del samp_read_data[read]
    print("[STDOUT]: remaining reads: {}".format(len(samp_read_data.keys())))
    return samp_read_data

def load_transcriptome():
    transcript_dict = {}
    with open(args.transcript_bed,'r') as f:
        while True:
            line = f.readline()
            if line:
                if line.startswith('#'):
                    pass
                else:
                    line = line.strip()
                    fields = line.split("\t")
                    try:
                        if len(fields) > 0:
                            transcript_type,chrom,start_region_start,start_region_end,end_region_start,end_region_end,splice2_start_start,splice2_start_end,splice2_end_start,splice2_end_end = str(fields[5]),str(fields[0]),int(fields[1]),int(fields[2]),int(fields[3]),int(fields[4]),int(fields[6]),int(fields[7]),int(fields[8]),int(fields[9])
                            transcript_dict[transcript_type] = (chrom,start_region_start,start_region_end,end_region_start,end_region_end,splice2_start_start,splice2_start_end,splice2_end_start,splice2_end_end)
                    except IndexError:
                        print("[STDERR]: error in transcriptome bed format, check input file")
                        print(line,"\n",fields)
                        exit()
            else:
                break
    print("[STDOUT]: view of input transcriptome:")
    print("\tTranscript_name\tchrom\tstart_region_start\tstart_region_end\tend_region_start\tend_region_end\tsplice2_start_start\tsplice2_start_end\tsplice2_end_start\tsplice2_end_end")
    for k,v in transcript_dict.items():
        print("\t",k,v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],sep='\t')
    print("")
    return transcript_dict


def get_df_counts(transcripts:list,sample_data:pd.DataFrame,sample_name:str):
    print("[STDOUT]: Counts for each transcript type in sample: {}".format(sample_name))
    if sample_name == "All_Samples" or sample_name == "All_Samples_Normalized":
        if sample_name == "All_Samples":
            tsv_fname = "raw_category_summary.tsv"
        else:
            tsv_fname = "normalized_category_summary.tsv"
        category_sum_file = open(tsv_fname,'w+')
        header_str = "Transcript_name\tcount\tpercent"
        print(header_str)
        category_sum_file.write(header_str + "\n")
        total_len = len(sample_data)
        round_depth = 3
        for transcript in transcripts:
            transcript_str = "{}\t{}\t{}".format(transcript,len(sample_data[sample_data['value'] == transcript]),round(100*(len(sample_data[sample_data['value'] == transcript])/total_len),round_depth))
            print(transcript_str)
            category_sum_file.write(transcript_str + "\n")
        category_sum_file.close()
        print("")
    elif sample_name[-3:] == "raw":
        df_out = pd.DataFrame(columns=['Transcript_name','count'])
        print("Transcript_name\tcount")
        for k,v in sample_data.read_classification_dict.items():
            df_out = pd.concat([df_out,pd.DataFrame({'Transcript_name':k,'count':len(v)},index=[0])],ignore_index=True)
            print(k,len(v))
    elif sample_name[-10:] == "normalized":
        df_out = pd.DataFrame(columns=['Transcript_name','normalized_count'])
        for transcript_type in transcripts:
            count = len(sample_data[sample_data['value'] == transcript_type])
            df_out = pd.concat([df_out,pd.DataFrame({'Transcript_name':transcript_type,'normalized_count':count},index=[0])],ignore_index=True)
            print(transcript_type,count)
    # write outputs for the later options in a non-redundant format
    if sample_name[-3:] == "raw" or sample_name[-10:] == "normalized":
        df_out.reset_index(drop=True,inplace=True)
        df_out.to_csv(sample_name + "_category_summary.tsv",sep='\t',header=True,index=False)
    print("")

def generate_individual_count_plot(sample_combined_data:pd.DataFrame,sample_name:str):
    plt.figure(figsize =(15,10),dpi=300)
    g = sns.countplot(data=sample_combined_data,x='value')
    g.set(xlabel="Transcript Type",ylabel="Raw Count")
    g.set_title(sample_name + " Raw Transcript Counts")
    g.set_xticklabels(g.get_xticklabels(),rotation=45,horizontalalignment='right')
    plt.savefig(sample_name + "_raw_transcript_counts.png",bbox_inches='tight')
    plt.close()
    
def split_df_by_counts(df:pd.DataFrame)->dict[str:pd.DataFrame]:
    """
    split df into low and high count dfs @ arbitrary threshold of 50 reads
    """
    grouped_df = df.groupby('value').count()
    low_count_df = df[df['value'].isin(grouped_df[grouped_df['read_name'] < 50].index)]
    high_count_df = df[df['value'].isin(grouped_df[grouped_df['read_name'] > 50].index)]
    return {"low_count":low_count_df,"high_count":high_count_df}

def generate_combined_count_plot(combined_df:pd.DataFrame):
    """
    Method to split raw counts data and plot based on simple count threshold
    """
    split_dfs = split_df_by_counts(combined_df)
    print("[STDOUT]: plotting raw counts")
    for name,df in split_dfs.items():
        plt.figure(figsize=(30,20),dpi=300)
        if use_custom_palatte:
            g = sns.countplot(data=df,x='value',hue='sample',palette=get_custom_experimental_palatte())
        else:
            g = sns.countplot(data=df,x='value',hue='sample',palette="bright")
        g.set(xlabel="Transcript Type",ylabel="Raw Count")
        g.set_title("All Sample Raw Transcript Counts")
        g.set_xticklabels(g.get_xticklabels(),rotation=45,horizontalalignment='right')
        outf = name + "_samples_compared_raw_transcript_counts.png"
        plt.savefig(outf,bbox_inches='tight')
        plt.close()

def combine_and_melt_frames(transcripts:list,sample_data:object,sample_name:str):
    # prepare for melting by adding column transcript type to respective dfs
    # transcripts = ['HBeAg','pgRNA','Core_WT','HBsAg_L','HBsAg_M','HBsAg_S','X_Canonical','X_Long','Other','HBeAg_Pol_Fusion_0','HBeAg_Pol_Fusion_1','Core_S_fusion','Core_Pol_fusion']
    for transcript_type in transcripts:
        sample_data.read_classification_dict[transcript_type]['transcript_type'] = transcript_type
    ### combine and melt ###
    combined_df = pd.concat((sample_data.read_classification_dict.values()),ignore_index=True)
    combined_df = pd.melt(combined_df,id_vars=['chr','read_name','read_start','read_end','sequence'],value_vars=['transcript_type'])
    combined_df['sample'] = sample_name
    print(combined_df.columns)
    return combined_df

def parse_sample_weights(samp_weights:str):
    print("[STDOUT]: RPM initial sample weights: {}".format(samp_weights))
    weight_list = samp_weights.split(',')
    min_weight = float(min(weight_list))
    transformed_weight_list = [round(min_weight/float(x),4) for x in weight_list]
    print("[STDOUT]: Transformed sample weights: {}".format(transformed_weight_list))
    return transformed_weight_list

def drop_reads_by_weight(df:pd.DataFrame,ordered_samp_weights:list,samfiles:list):
    categories = df.value.unique().tolist()
    cnt = 0
    for filename in samfiles:
        # find only reads coming from correct file
        fname = filename[:-4]
        # drop reads by weight for each category, i.e. downsample equally for each category
        for category in categories:
            # get indexes from reads only from this file & category
            idxlist = df.index[(df['sample'] == fname) & (df['value'] == category)].tolist()
            len_idx_list = len(idxlist)
            downsampled_count = int(round(len_idx_list * float(ordered_samp_weights[cnt]),0))
            if downsampled_count > 0 and downsampled_count < len(idxlist):
                drop_list = random.sample(idxlist,len_idx_list - downsampled_count)
                if debug:
                    print("[STDOUT]: dropping {} reads from category {}, which originally had {} reads based input read weight of  {}".format(len(idxlist)-downsampled_count,category,len(idxlist),ordered_samp_weights[cnt]))
                    print("[STDOUT]: dropping index #'s:")
                    print(drop_list)
                df = df.drop(drop_list)
        # downsample the next category
        cnt += 1
    # output
    return df

def get_custom_experimental_palatte():
    PuBuGnPal = sns.color_palette("PuBuGn",10)
    Pu = sns.color_palette("PuOr",10)
    Rd = sns.color_palette("RdBu", 10)
    Dark2 = sns.color_palette("Dark2",10)
    customPalatte = [Pu[::-1][:3],PuBuGnPal[::-1][:4],Rd[:3],Dark2[:3]]
    customPalatte = list(chain.from_iterable(customPalatte))
    return customPalatte    

def get_and_plot_normalized_counts(output_df:pd.DataFrame,samfiles:list):
    """
    Method to downsample reads by weight and plot normalized counts after splitting by normalized counts
    """
    print("[STDOUT]: plotting normalized counts")
    ordered_samp_weights = parse_sample_weights(args.sample_weights)
    new_df = drop_reads_by_weight(output_df,ordered_samp_weights,samfiles)
    split_dfs = split_df_by_counts(new_df)
    for name,df in split_dfs.items():
        plt.figure(figsize=(30,20),dpi=300)
        if use_custom_palatte:
            g = sns.countplot(data=df,x='value',hue='sample',palette=get_custom_experimental_palatte())
        else:
            g = sns.countplot(data=df,x='value',hue='sample',palette="bright")
        g.set_xticklabels(g.get_xticklabels(),rotation=45,horizontalalignment='right')
        g.set(xlabel="Transcript Type",ylabel="Normalized Transcript Count")
        g.set_title("All Sample Normalized Transcript Counts")
        outf = name + "_samples_compared_normalized_transcript_counts.png"
        plt.savefig(outf,bbox_inches='tight')
        plt.close()
    return new_df  

def return_markers(sample_group_list:list)->list[str]:
    use_markers = []
    marker_palatte = {'A2A004_A3A006_A4A014_HC_D85-141':'o','88A010_95A010_D57':'s','A2A004_A3A006_D351':'d','A2A004_D323':'.'}
    for group in sample_group_list:
        if group in marker_palatte.keys():
            use_markers.append(marker_palatte[group])
        else:
            pass
    return use_markers

def make_facet_grid_of_categories(df:pd.DataFrame,sample_groups:dict[str:tuple],category:str):
    if sample_groups != None:
        # add a new column to df using sample_group key to indicate group...
        for group_name,samples in sample_groups.items():
            for sample in samples:
                df.loc[df['sample'] == sample,'sample_group'] = group_name
    else:
        df['sample_group'] = "All_Samples"

    print("[STDOUT]: making facet grid scatterplots of category:",category,"for sample groups:",sample_groups)
    plt.figure(dpi=300)
    df_new = df[['read_start','read_end','value','sequence','sample_group']]
    if sample_groups != None: 
        # get markers for each group
        col_palatte = {'A2A004_A3A006_A4A014_HC_D85-141':'black','88A010_95A010_D57':'grey','A2A004_A3A006_D351':'blue','A2A004_D323':'green'}    
        use_markers = return_markers(list(df_new['sample_group'].unique()))
        sns.lmplot(data=df_new,fit_reg=False,x="read_start",y="read_end",col="value",hue='sample_group',palette=col_palatte,markers=use_markers,col_wrap=2,ci=None,height=8,aspect=1)
        outf = category + "_sample_grouped_transcript_categorization_grid.png"
    else:
        # setup custom palatte
        scatter_palette = sns.color_palette("hls",len(df_new['value'].unique()))[::-1]
        # plot
        sns.lmplot(data=df_new,fit_reg=False,x="read_start",y="read_end",col="value",hue="value",col_wrap=5,palette=scatter_palette,ci=None,height=8,aspect=1)
        outf = "transcript_categorization_grid.png"
    plt.savefig(outf,bbox_inches='tight')
    plt.close()
    if sample_groups == None:
        # build combined scatter plot from the same data:
        make_combined_scatterplot_of_categories(df_new,scatter_palette)

def make_combined_scatterplot_of_categories(df_new:pd.DataFrame,pal):
    print("[STDOUT]: making combined scatterplot of all categories")
    f, ax = plt.subplots(figsize=(12,12),dpi=300)
    sns.scatterplot(x="read_start", y="read_end",hue="value",palette=pal,linewidth=0,data=df_new,ax=ax)
    ax.set(xlabel="read alignment start (bp)",ylabel="read alignment end (bp)")
    outf = "transcript_categorization_combined_scatter.png"
    plt.savefig(outf,bbox_inches='tight')
    plt.close()

def quantify_categories(df:pd.DataFrame,datatype:str):
    print("[STDOUT]: quantifying categories in",datatype,"dataframe")
    outdf = pd.DataFrame(columns=['sample','category','count'])
    # by sample
    for sample in set(df['sample']):
        # by category
        for catergory_name in set(df['value']):
            if debug:
                print("Sample:",sample,catergory_name,df[(df['value'] == catergory_name) & (df['sample'] == sample)].count()['value'])
            outdf = pd.concat([outdf,pd.DataFrame({'sample':sample,'category':catergory_name,'count':df[(df['value'] == catergory_name) & (df['sample'] == sample)].count()['value']},index=[0])],ignore_index=True)
    if datatype == "raw":
        outdf.to_csv("raw_transcript_counts_x_sample.csv",sep=',',header=True,index=False)
    elif datatype == "normalized":
        outdf.to_csv("normalized_transcript_counts_x_sample.csv",sep=',',header=True,index=False)

def get_bam_of_read_category(sample:str,read_list:list,category:str):
    # get reads in list:S
    read_list = [x.strip() for x in read_list]
    # open combined bam:
    bam = pysam.AlignmentFile(args.combined_bam,'rb')
    name_out = "samp_" + str(sample) + "_" + category + "_reads.bam"
    print("[STDOUT]: creating output bam file:",name_out,"with category:",category,"for sample:",sample)
    bam_out = pysam.AlignmentFile(name_out,'wb',template=bam)
    for read in bam.fetch():
        if read.query_name in read_list:
            bam_out.write(read)
    bam_out.close()
    bam.close()
    # index bam
    print("[STDOUT]: indexing output bam file:",name_out)
    pysam.index(name_out)

if __name__ == '__main__':
    print("*** HBV transcript categorization analysis ***")
    print("[STDOUT]: sam file list: {}".format(args.samfiles))
    print("[STDOUT]: transcript bed: {}".format(args.transcript_bed))
    print("[STDOUT]: minimum quality: {}".format(args.min_qual))
    print("[STDOUT]: debug: {}".format(args.debug))
    print("")
    # def readin hbv transcriptome
    transcriptome_dict = load_transcriptome()
    transcriptome_keys = list(transcriptome_dict.keys())
    transcriptome_keys.append("Other")
    raw_output_df = pd.DataFrame(columns=['chr','read_name','read_start','read_end','variable','value','sequence'])
    for sam in samfiles:
        print("[STDOUT]: categorizing reads in {}".format(sam))
        # should contain categorized dfs of all reads in the sam file:
        samp_data = categorize_reads(sam,transcriptome_dict)
        # get stats and emit a tsv of the counts:
        get_df_counts(transcriptome_keys,samp_data,sam[:-4] + "_raw")
        # combine dataframes:
        combined_df = combine_and_melt_frames(transcriptome_keys,samp_data,sam[:-4])
        # make plots:
        if len(combined_df.index) > 0:
            generate_individual_count_plot(combined_df,sam[:-4])
            # append to output df
            raw_output_df = pd.concat((raw_output_df,combined_df),ignore_index=True)
    # make combined plot:
    print("[STDOUT]: combined dataframe:")
    print(raw_output_df.head())
    # get final df counts
    get_df_counts(transcriptome_keys,raw_output_df,"All_Samples")
    # df count plot:
    generate_combined_count_plot(raw_output_df)
    # make normalized plot of transcript counts:
    norm_output_df = get_and_plot_normalized_counts(raw_output_df,samfiles)
    norm_output_df_for_cats = norm_output_df.copy()
    # get individual normalized counts:
    for sam in samfiles:
        samp_norm_df = norm_output_df[norm_output_df['sample'] == sam[:-4]]
        if len(samp_norm_df.index) > 0:
            get_df_counts(transcriptome_keys,samp_norm_df,sam[:-4] + "_normalized")
    # get all samples normalized counts:    
    get_df_counts(transcriptome_keys,norm_output_df,"All_Samples_Normalized")
    # facet grid of categories:
    # all categories in all samples
    make_facet_grid_of_categories(norm_output_df.copy(),None,"all")
    # select categories in select groups:
    HBxAg_cats = ['X_Canonical','X_Long','X_splice1','X_splice2','X_splice3','X_splice4','X_splice5','X_splice6','X_2X']
    HBsAg_cats = ['HBsAg_L','HBsAg_M','HBsAg_S','HBsAg_L_iDNA','HBsAg_M_iDNA','HBsAg_S_iDNA']
    HBeAg_cats = ['HBeAg_wt','HBeAg_ORF','HBeAg_Cys-','HBeAg_splice1','HBeAg_splice2']
    HBcoreAg_cats = ['Core_wt','Core_ORF','Core-Cys-','Core_splice1','Core_splice2']
    Pol_cats = ['Polymerase','Core_wt']
    # unused
    # 89A008_HC.ccs_polish-flnc.fa.algn.srt.sam
    # 89A008_Day_141.ccs_polish-flnc.fa.algn.srt.sam
    # 89A008_Day_379.ccs_polish-flnc.fa.algn.srt.sam
    sample_grouping = { "A2A004_A3A006_A4A014_HC_D85-141":('A2A004_HC.ccs_polish-flnc.fa.algn.srt','A2A004_d85.ccs_polish-flnc.fa.algn.srt','A3A006_HC.ccs_polish-flnc.fa.algn.srt','A3A006_d141.ccs_polish-flnc.fa.algn.srt','A4A014_HC.ccs_polish-flnc.fa.algn.srt'),"88A010_95A010_D57":('88A010_57.ccs_polish-flnc.fa.algn.srt','95A010_d57.ccs_polish-flnc.fa.algn.srt'),"A2A004_A3A006_D351":('A2A004_Day_351.ccs_polish-flnc.fa.algn.srt','A3A006_Day_351.ccs_polish-flnc.fa.algn.srt'),"A2A004_D323":('A2A004_d323.ccs_polish-flnc.fa.algn.srt')}

    # make combined scatter plot of all categories:
    for cat_list in [HBxAg_cats,HBsAg_cats,HBeAg_cats,HBcoreAg_cats,Pol_cats]:
        # filtered to only categories of interest
        cat_df = norm_output_df_for_cats[norm_output_df_for_cats['value'].isin(cat_list)]
        if len(cat_df) > 0:
            make_facet_grid_of_categories(cat_df,sample_grouping,cat_list[0])
    # generate tsv output of raw counts
    new_col_idx = ['chr','read_start','read_end','value','sample','read_name','sequence']
    raw_output_df = raw_output_df.reindex(columns=new_col_idx)
    raw_output_df.to_csv("raw_read_categorization_x_sample.csv",sep=',',header=True,index=False)
    norm_output_df = norm_output_df.reindex(columns=new_col_idx)
    norm_output_df.to_csv("normalized_read_categorization_x_sample.csv",sep=',',header=True,index=False)
    # quantify categories:
    quantify_categories(raw_output_df,"raw")
    quantify_categories(norm_output_df,"normalized")
    # replot the analysis with bokeh:
    replot_analysis()

    if args.combined_bam != "":
        samples = list(set(raw_output_df['sample']))
        for sample in samples:
            # get bam of reads for a specific category:
            category = "Other"
            sample_read_list = list(raw_output_df[(raw_output_df['value'] == category) & (raw_output_df['sample'] == sample)]['read_name'])
            print("[STDOUT]: found {} reads in category {} for sample {}".format(len(sample_read_list),category,sample))
            if len(sample_read_list) > 0:
                get_bam_of_read_category(sample.replace('.bam',''),sample_read_list,category)
        # merged of all samples with Other reads:
        print("[STDOUT]: getting bam of read category:","Other")
        comb_sample_read_list = list(raw_output_df[raw_output_df['value'] == "Other"]['read_name'])
        get_bam_of_read_category("all_samples",comb_sample_read_list,'Other')
            
