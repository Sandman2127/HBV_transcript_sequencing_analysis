#!/opt/anaconda/envs/py3/bin/python3
import pandas as pd
#import seaborn as sns
import matplotlib.pyplot as plt
from bokeh.models import ColumnDataSource
#from bokeh.palettes import Category20,Bright6,Viridis256
from bokeh.plotting import figure, show
#from bokeh.transform import factor_cmap
from bokeh.io import export_png
from math import pi

sample_color_map = {
    'A2A004_HC':'black',
    'A2A004_d85':'black',
    'A2A004_d323':'blue',
    'A2A004_d351':'blue',
    'A4A014_HC':'black',
    'A3A006_HC':'black',
    'A3A006_d141':'black',
    'A3A006_d351':'blue',
    '88A010_d57':'grey',
    '95A010_d57':'grey',
    '89A008_HC':'orange',
    '89A008_d141':'orange',
    '89A008_d379':'orange',
    701:'black',
    705:'blue',
    709:'orange',
    710:'red',
    712:'purple'
    }

transcript_color_map = {
    'HBeAg_wt': 'lightcoral',
    'HBeAg_ORF': 'salmon',
    'HBeAg_sp6': 'darkorange',
    'HBeAg_sp11': 'coral',
    'HBeAg_sp1_Cys-': 'tomato',
    'HBeAg-S_sp5': 'lightsalmon',
    'HBeAg-S_sp9': 'orangered',
    'POL': 'gold',
    'Core_wt': 'lightpink',
    'Core_ORF': 'hotpink',
    'Core_sp6': 'deeppink',
    'Core_sp11': 'palevioletred',
    'Core-sp1_Cys-': 'crimson',
    'Core-S_sp5': 'pink',
    'Core-S_sp9': 'mediumvioletred',
    'S-Core': 'plum',
    'HBsAg_2X_2': 'lightsteelblue',
    'HBsAg_2X_3': 'cornflowerblue',
    'HBsAg_2X_4': 'royalblue',
    'HBsAg_L': 'skyblue',
    'HBsAg_L_iDNA': 'deepskyblue',
    'HBsAg_M': 'dodgerblue',
    'HBsAg_M_iDNA': 'blue',
    'HBsAg_S': 'mediumblue',
    'HBsAg_S_iDNA': 'navy',
    'X_long': 'mediumorchid',
    'X_canonical': 'rebeccapurple',
    'X_splice1': 'indigo',
    'X_splice2': 'blueviolet',
    'X_splice3': 'darkviolet',
    'X_splice4': 'purple',
    'X1_splice5': 'mediumorchid',
    'X_splice6': 'darkorchid',
    'X_trunc1': 'olive',
    'X_trunc2': 'olivedrab',
    'X_trunc3': 'darkolivegreen',
    'X_2X':'tan',
    'X2_canonical':'peru',
    'X2_long':'burlywood',
    'X2_long_iDNA':'wheat',
    'X2_splice5':'sandybrown',
    'Other':'grey'
    } 

transcript_order = list(transcript_color_map.keys())

sample_rename_map = {
     'A2A004_HC':'A2A004_HC',
     'A2A004_d85':'A2A004_d85',
     'A2A004_d323':'A2A004_d323',
     'A2A004_Day_351':'A2A004_d351',
     'A4A014_HC':'A4A014_HC',
     'A3A006_HC':'A3A006_HC',
     'A3A006_d141':'A3A006_d141',
     'A3A006_Day_351':'A3A006_d351',
     '88A010_57':'88A010_d57',
     '95A010_d57':'95A010_d57',
     '89A008_HC':'89A008_HC',
     '89A008_Day_141':'89A008_d141',
     '89A008_Day_379':'89A008_d379'
    }

stuff_to_combine = {
                    'HBsAg_L':('HBsAg_L1','HBsAg_L2'),
                    'HBsAg_L_iDNA':('HBsAg_L1_iDNA','HBsAg_L2_iDNA'),
                    'HBsAg_M':('HBsAg_M1','HBsAg_M2'),
                    'HBsAg_M_iDNA':('HBsAg_M2_iDNA','HBsAg_M1_iDNA'),
                    'HBsAg_S':('HBsAg_S1','HBsAg_S2'),
                    'HBsAg_S_iDNA':('HBsAg_S1_iDNA','HBsAg_S2_iDNA'),
                    'X_long':('X1_long','X2_long'),
                    'X_canonical':('X1_canonical','X2_canonical'),
                    'X_splice5':('X1_splice5','X2_splice5'),
                    'X_trunc1':('X1_trunc1','X2_trunc1'),
                    'X_trunc2':('X1_trunc2','X2_trunc2'),
                    'X_trunc3':('X1_trunc3','X2_trunc3')
            }
                # Human data x trunc ^^
                # 'X_trunc1':('X1_trunc1','X2_trunc1'),
                # 'X_trunc2':('X1_trunc2','X2_trunc2'),
                # 'X_trunc3':('X1_trunc3','X2_trunc3')
                # Chimp data specific
                #     'X_trunc1':('X_trunc1','X1_trunc1'),
                #     'X_trunc2':('X_trunc2','X1_trunc2'),
                #     'X_trunc3':('X_trunc3','X1_trunc3')
                #    }


def get_combined_category_count(df:pd.DataFrame,sample:str,categories_to_combine:tuple[str,str,str],count_all:bool)->int:
    """
    Combine categories based on
    
    input dataframe: pd.DataFrame
    sample: str
    categories_to_combine: tuple(str,
    """
    cnt = 0
    for cat in categories_to_combine:
        if count_all:
            tupA = tuple(df[df['category'] == cat]['count'])
        else:
            tupA = tuple(df[(df['category'] == cat) & (df['sample'] == sample)]['count'])
        # checks to ensure the tuple has correct length and deals with the situation correctly
        if len(tupA) == 1:
            cnt += int(tupA[0])
        elif len(tupA) == 0:
            cnt += 0
        else:
            exit()
    return cnt

def combine_and_drop_subcategories(input_df:pd.DataFrame,samples:tuple[str,str,str],stuff_to_combine:dict):
    """
    # input df is:
    #	category	count	sample
    # 0	Core_wt	0	95A010_d57
    # 1	X2_long_iDNA	0	95A010_d57
    """
    print("Length of input df:",len(input_df.index))
    for output_category,input_cats in stuff_to_combine.items():
        # do the combinations on the sample level naturally
        for sample in samples:
            category_combined_count = get_combined_category_count(input_df,sample,input_cats,False)
            # add to primary df
            input_df.loc[len(input_df)] = [output_category,category_combined_count,sample] 
    
    print("Length of input df after combination:",len(input_df.index))
    # reset your index and drop the old categories
    input_df.reset_index(inplace=True,drop=True)
    idx_to_drop = []
    # start dropping the old categories:
    for output_category,input_cats in stuff_to_combine.items():
        for cat in input_cats:
            idx_to_drop.extend(list(input_df[input_df['category'] == cat].index))
    # drop these index positions
    input_df.drop(set(idx_to_drop),inplace=True)
    print("Length of input df after combination & old category drop:",len(input_df.index))
    return input_df

def bokeh_bar_plot(sample_combined_data:pd.DataFrame,sample_name:str,transcript_order:list,color_map:dict):
    """
    transform input df to plot single plot based on sample name and make bar plot of categories:
    """
    if sample_name == 'all':
        subsetted_df = sample_combined_data
    else:
        subsetted_df = sample_combined_data[sample_combined_data['sample'] == sample_name]
    # make fig:
    p = figure(x_range=transcript_order,height=1000,width=1400,toolbar_location=None,title=str(sample_name) + " RPM normalized transcript counts")
    # plot by bar:
    for category in transcript_order:
        try:
            y_val = subsetted_df[subsetted_df['category'] == category]['count'].values[0]
        except IndexError:
            y_val = 0
        p.vbar(x=[category], top=[y_val],width=0.6,color=color_map[category])
    # title:
    p.title.text_font_size = '20pt'
    font_type = "helvetica"
    p.xaxis.axis_label_text_font = font_type
    p.yaxis.axis_label_text_font = font_type
    # Configure Axis Major Labels (ticks)
    p.xaxis.major_label_text_font = font_type
    p.yaxis.major_label_text_font = font_type
    # get the axes right
    p.yaxis.axis_label = "Normalized Transcript Count"
    p.yaxis.axis_label_text_font_style = "normal"
    # Setting the text size for the axis labels
    p.yaxis.axis_label_text_font_size = "20pt"
    p.xaxis.axis_label_text_font_size = "15pt"
    # Setting the text size for the axis major labels (the ticks)
    p.yaxis.major_label_text_font_size = "20pt"
    p.xaxis.major_label_text_font_size = "15pt"
    p.xaxis.major_label_orientation = pi/4
    show(p)
    # Export as PNG with a higher resolution
    export_png(p, filename=str(sample_name) + "_bar_plot.png")


def combine_category_variables(orig_category:str,combine_categories:dict):
    try:
        value = transcript_color_map[orig_category]
        return orig_category
    except KeyError:
        for comb_category,constituent_cats in combine_categories.items():
            if orig_category in constituent_cats:
                return comb_category
            else:
                pass
    print("Failed on",orig_category)
    return orig_category

def map_to_color(category:str):
    try:
        return transcript_color_map[category]
    except KeyError:
        print(category)
        return 'black'

def scatter(df,input_line_color:str,sample_comparison_name:str,alpha:float):
    p = figure(width=1200, height=1200, title = str(sample_comparison_name) + " transcript position x type")
    source = ColumnDataSource(df)
    p.scatter("read_start","read_end",source=source,fill_alpha=alpha,size=8,legend_group='category',color=input_line_color,line_color=input_line_color)
    # get fonts same
    font_type = "helvetica"
    p.xaxis.axis_label_text_font = font_type
    p.yaxis.axis_label_text_font = font_type
    # Configure Axis Major Labels (ticks)
    p.xaxis.major_label_text_font = font_type
    p.yaxis.major_label_text_font = font_type
    # Configure Legend
    p.legend.label_text_font = font_type

    # title data:
    p.title.text_font_size = '20pt'
    p.xaxis.axis_label = 'Read Start (bp)'
    p.yaxis.axis_label = 'Read End (bp)'
    p.xaxis.axis_label_text_font_style = "normal"
    p.yaxis.axis_label_text_font_style = "normal"
    # Setting the text size for the axis labels:
    p.yaxis.axis_label_text_font_size = "16pt"
    p.xaxis.axis_label_text_font_size = "16pt"
    # Setting the text size for the axis major labels (the ticks):
    p.yaxis.major_label_text_font_size = "18pt"
    p.xaxis.major_label_text_font_size = "18pt"
    # legend:
    p.legend.location = "top_right"
    p.legend.label_text_font_size = '14pt'  # Adjust the font size of the text in the legend
    p.legend.glyph_height = 28  # Adjust the height of the legend markers
    p.legend.glyph_width = 28  # Adjust the width of the legend markers
    p.legend.spacing = 2  # Adjust spacing between legend entries
    p.legend.padding = 2  # Adjust padding inside the legend box
    p.legend.margin = 2  # Adjust margin around the legend box
    p.add_layout(p.legend[0], 'right')
    # Remove minor grid lines
    # p.xgrid.minor_grid_line_color = None
    # p.ygrid.minor_grid_line_color = None
    # Optionally, remove major grid lines too if needed
    # p.xgrid.visible = False
    # p.ygrid.visible = False
    show(p)
    export_png(p,filename=str(sample_comparison_name) + "_alpha_" + str(alpha)  + "_category_position_x_type.png")

def replot_analysis():
    # start replotting all previous results with updated figure inputs
    norm_transcript_counts = pd.read_csv("normalized_transcript_counts_x_sample.csv")
    if 705 in set(norm_transcript_counts['sample']):
        norm_transcript_counts['sample_new'] = norm_transcript_counts['sample']
    else:
        norm_transcript_counts['sample_new'] = norm_transcript_counts['sample'].apply(lambda sample: sample.replace('.ccs_polish-flnc.fa.algn.srt','')).map(sample_rename_map)
    norm_transcript_counts.drop(columns=['sample'],inplace=True)
    norm_transcript_counts.rename(columns={"sample_new": "sample", "category": "category",'counts':'counts'},inplace=True)
    # get renamed samples
    samples = set(norm_transcript_counts['sample'])
    # recombine and drop subcategory work:
    norm_transcript_counts = combine_and_drop_subcategories(norm_transcript_counts.copy(),samples,stuff_to_combine)
    # get categories
    possible_categories = set(norm_transcript_counts['category'])
    transcripts_in_order = []
    for item in transcript_order:
        if item in possible_categories:
            transcripts_in_order.append(item)
    # build bar plots
    for sample in samples:
        bokeh_bar_plot(norm_transcript_counts,sample,transcripts_in_order,transcript_color_map)
    # get read positioning data and do some transformations
    norm_read_alignment_and_cats_x_sample = pd.read_csv('normalized_read_categorization_x_sample.csv')
    if 705 in set(norm_read_alignment_and_cats_x_sample['sample']):
        norm_read_alignment_and_cats_x_sample['sample_new'] = norm_read_alignment_and_cats_x_sample['sample']
    else:
        norm_read_alignment_and_cats_x_sample['sample_new'] = norm_read_alignment_and_cats_x_sample['sample'].apply(lambda sample: sample.replace('.ccs_polish-flnc.fa.algn.srt','')).map(sample_rename_map)
    norm_read_alignment_and_cats_x_sample['ncategory'] = norm_read_alignment_and_cats_x_sample['value'].apply(lambda category: combine_category_variables(category,stuff_to_combine))
    # drop old columns
    norm_read_alignment_and_cats_x_sample.drop(columns=['value','sample'],inplace=True)
    norm_read_alignment_and_cats_x_sample.rename(columns={'ncategory':'category','sample_new':'sample'},inplace=True)
    # map the category norm_read_alignment_and_cats_x_samples to colors in df
    norm_read_alignment_and_cats_x_sample['transcript_category_color'] = norm_read_alignment_and_cats_x_sample['category'].apply(lambda cat: map_to_color(cat))
    # map the sample names to colors
    norm_read_alignment_and_cats_x_sample['sample_category_color'] = norm_read_alignment_and_cats_x_sample['sample'].map(sample_color_map)
    # make for different alpha's
    for alpha in [0,0.4]:
        for sample in samples:
            df_in = norm_read_alignment_and_cats_x_sample[ norm_read_alignment_and_cats_x_sample['sample'] == sample ]
            scatter(df_in,'transcript_category_color',sample,alpha)
        # one full scatter for all
        scatter(norm_read_alignment_and_cats_x_sample,'transcript_category_color',"All_Samples",alpha)

if __name__ == '__main__':
    replot_analysis()