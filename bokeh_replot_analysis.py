#!/opt/anaconda/envs/py3/bin/python3
import pandas as pd
#import seaborn as sns
import matplotlib.pyplot as plt
from bokeh.models import ColumnDataSource, Label, Span, Legend, LegendItem
#from bokeh.palettes import Category20,Bright6,Viridis256
from bokeh.plotting import figure, show
#from bokeh.transform import factor_cmap
from bokeh.io import export_png
from math import pi

# Dimensions in pixels (6.5 inches * 300 dpi by 2 inches * 300 dpi)
plot_width = 1950
plot_height = 600

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
    'Core_wt': 'black',
    'Core_ORF': 'hotpink',
    'Core_sp6': 'deeppink',
    'Core_sp11': 'palevioletred',
    'Core-sp1_Cys-': 'crimson',
    'Core-S_sp5': 'pink',
    'Core-S_sp9': 'mediumvioletred',
    'S-Core': 'plum',
    'POL': 'gold',
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

sample_title_map = {
    "All_Samples_Chimp":['All chimpanzee (HBeAg+ and HBeAg-) HBV transcripts combined'],
    'A2A004_HC':['HBeAg+ chimpanzee A2A004 transcript counts pre-study (HC)','HBeAg+ chimpanzee A2A004 HBV transcript positions pre-study (HC)'],
    'A2A004_d85':['HBeAg+ chimpanzee A2A004 transcript counts after ETV lead-in (d85)','HBeAg+ chimpanzee A2A004 HBV transcript positions after ETV lead-in (d85)'],
    'A2A004_d323':['HBeAg+ chimpanzee A2A004 transcript counts at end of ARC-520 + ETV treatment (d323)','HBeAg+ chimpanzee A2A004 HBV transcript positions at end of ARC-520 + ETV treatment (d323)'],
    'A2A004_d351':['HBeAg+ chimpanzee A2A004 transcript counts post ARC-520 + ETV treatment (d351)','HBeAg+ chimpanzee A2A004 HBV transcript positions post ARC-520 + ETV treatment (d351)'], 
    '89A008_HC':['HBeAg+ chimpanzee 89A008 HBV transcript counts pre-study (HC)','HBeAg+ chimpanzee 89A008 HBV transcript positions pre-study (HC)'], 
    '89A008_d141':['HBeAg-transitional chimpanzee 89A008 HBV transcript counts after ETV lead-in (d141)','HBeAg-transitional chimpanzee 89A008 HBV transcript positions after ETV lead-in (d141)'],
    '89A008_d379':['HBeAg- chimpanzee 89A008 HBV transcript counts post siRNA + ETV treatment (d379)','HBeAg- chimpanzee 89A008 HBV transcript positions post siRNA + ETV treatment (d379)'], 
    'A3A006_HC':['HBeAg+ chimpanzee A3A006 transcript counts pre-study (HC)','HBeAg+ chimpanzee A3A006 HBV transcript positions pre-study (HC)'], 
    'A3A006_d141':['HBeAg+ chimpanzee A3A006 transcript counts after ETV lead-in (d141)','HBeAg+ chimpanzee A3A006 HBV transcript positions after ETV lead-in (d141)'],
    'A3A006_d351':['HBeAg+ chimpanzee A3A006 transcript counts post ARC-520 + ETV treatment (d351)','HBeAg+ chimpanzee A3A006 HBV transcript positions post ARC-520 + ETV treatment (d351)'], 
    'A4A014_HC':['HBeAg+ chimpanzee A4A014 transcript counts pre-study (HC)','HBeAg+ chimpanzee A4A014 HBV transcript positions pre-study (HC)'], 
    '95A010_d57':['HBeAg- chimpanzee 95A010 HBV transcript counts after ETV lead-in (d57)','HBeAg- chimpanzee 95A010 HBV transcript positions after ETV lead-in (d57)'],
    '88A010_d57':['HBeAg- chimpanzee 88A010 HBV transcript counts after ETV lead-in (d57)','HBeAg- chimpanzee 88A010 HBV transcript positions after ETV lead-in (d57)'],
    "All_Samples_Homo_sapiens":['All (previously HBeAg+ and HBeAg-) end-of-study human patient HBV transcripts combined'],
    '701':['na','na'],
    '705':['HBeAg- patient 705 HBV transcript counts 19.6 months post ARC-520 + ETV treatment','HBeAg- patient 705 HBV transcripts 19.6 months post ARC-520 + ETV treatment'],
    '709':['na','na'],
    '710':['HBeAg- patient 710 (HBeAg+ pre-study) HBV transcript counts 19.7 months post ARC-520 + ETV treatment','HBeAg- patient 710 (HBeAg+ pre-study) HBV transcripts 19.7 months post ARC-520 + ETV treatment '],
    '712':['HBeAg- patient 712 HBV transcript counts 21.6 months post ARC-520 + ETV treatment','HBeAg- patient 712 HBV transcripts 21.6 months post ARC-520 + ETV treatment ']
}

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
    # special request to drop other unique categories X1_splice5,X1_splice6,X2_splice5
    input_df.drop(list(input_df[input_df['category'].isin(('X1_splice5','X1_splice6','X2_splice5'))].index),inplace=True)
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
    p = figure(x_range=transcript_order,height=plot_height,width=plot_width,toolbar_location=None,title=sample_title_map[str(sample_name)][0])
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
    # show the plot
    # show(p)
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

def scatter(df,input_line_color:str,sample_name:str,alpha:float,IS_CHIMP:bool):
    if "All" in str(sample_name):
        title_name =sample_title_map[str(sample_name)][0]
    else:
        title_name = sample_title_map[str(sample_name)][1]
    p = figure(width=plot_width, height=plot_height,title=title_name,toolbar_location=None)
    source = ColumnDataSource(df)
    p.scatter("read_start","read_end",source=source,fill_alpha=alpha,size=9,legend_group='category',color=input_line_color,line_color=input_line_color)
    # get fonts same
    font_type = "helvetica"
    p.xaxis.axis_label_text_font = font_type
    p.yaxis.axis_label_text_font = font_type
    # Configure Axis Major Labels (ticks)
    p.xaxis.major_label_text_font = font_type
    p.yaxis.major_label_text_font = font_type
    # Configure Legend
    p.legend.label_text_font = font_type
    # get the axes right
    p.y_range.start = 0 
    p.y_range.end = 6200
    # p.x_range.start = -60
    # p.x_range.end = 6200
    # title data:
    p.title.text_font_size = '20pt'
    p.xaxis.axis_label = 'Read Start (bp)'
    p.yaxis.axis_label = 'Read End (bp)'
    p.xaxis.axis_label_text_font_style = "normal"
    p.yaxis.axis_label_text_font_style = "normal"
    # Setting the text size for the axis labels:
    p.yaxis.axis_label_text_font_size = "18pt"
    p.xaxis.axis_label_text_font_size = "18pt"
    # Setting the text size for the axis major labels (the ticks):
    p.yaxis.major_label_text_font_size = "18pt"
    p.xaxis.major_label_text_font_size = "18pt"
    # legend:
    #p.legend.orientation = "horizontal"  # , normally vertical
    legend_len = df['category'].nunique()
    #p.legend.nrows = 2  # Adjust the number of rows in the legend
    if legend_len > 22:
        p.legend.ncols = 2
    else:
        p.legend.ncols = 1
    p.legend.location = 'center'
    p.legend.label_text_font_size = '16pt'  # Adjust the font size of the text in the legend
    p.legend.glyph_height = 20  # Adjust the height of the legend markers
    p.legend.glyph_width = 20  # Adjust the width of the legend markers
    p.legend.spacing = 2  # Adjust spacing between legend entries
    p.legend.padding = 2  # Adjust padding inside the legend box
    p.legend.margin = 2  # Adjust margin around the legend box
    p.add_layout(p.legend[0],"right")           #"above")       # 'right')
    x_max = max(df['read_start'])
    if IS_CHIMP:
        p = add_spans_and_labels(p,(10,625),2406,(2778,3807),5588,x_max)
    else:
        p = add_spans_and_labels(p,(10,656),2437,(2804,3871),5652,x_max)
    # Create a dummy renderer (off page) for the custom legend item (invisible glyph)
    custom_renderer = p.scatter([-100],[-100],marker='square',color="red",size=20,alpha=0.1)
    # Create the custom legend item
    legend_item = LegendItem(label="HBV PAS", renderers=[custom_renderer])
    # gather the existing legend and reorder
    existing_legend = p.legend[0]
    existing_legend.items.append(legend_item)
    transcript_order.append("HBV PAS")
    legend_mapping = {item.label['value']: item for item in existing_legend.items}
    ordered_items = [legend_mapping[label] for label in transcript_order if label in legend_mapping]
    existing_legend.items = ordered_items
    # Keep track of seen labels
    seen_labels = set()
    # Filter the legend items
    filtered_items = []
    for item in existing_legend.items:
        label = item.label['value']
        if label == "HBV PAS":
            if label not in seen_labels:
                seen_labels.add(label)
                filtered_items.append(item)
        else:
            filtered_items.append(item)
    # Update the legend items with the filtered list
    existing_legend.items = filtered_items
    # show the plot
    # show(p)
    export_png(p,filename=str(sample_name) + "_alpha_" + str(alpha)  + "_category_position_x_type.png")

def add_spans_and_labels(p:figure,span1_xs:tuple,span1_y:int,span2_xs:tuple,span2_y:int,max_x:int):
    # add spans
    span1 = Span(location=span1_y, dimension='width', line_color='red',line_width=16,line_alpha=0.1)
    span2 = Span(location=span2_y, dimension='width', line_color='red',line_width=16,line_alpha=0.1)
    p.add_layout(span1)
    p.add_layout(span2)
    if max_x > 5000:
        x_label_gain = 25
    else:
        x_label_gain = 10
    # add labels
    p.line(x=span1_xs,y=(span1_y+90,span1_y+90), line_color="black", line_width=2)
    label_hbs_1 = Label(x=span1_xs[0] + int((span1_xs[1]-span1_xs[0])/2) + x_label_gain, y=span1_y+110, text='HBsAg', text_color='black', text_align='center',text_font_size='18pt')
    p.line(x=span2_xs,y=(span2_y+90,span2_y+90), line_color="black", line_width=2)
    label_hbs_2 = Label(x=span2_xs[0] + int((span2_xs[1]-span2_xs[0])/2), y=span2_y+110, text='HBsAg', text_color='black', text_align='center',text_font_size='18pt')
    p.add_layout(label_hbs_1)
    p.add_layout(label_hbs_2)
    return p

def replot_analysis():
    """
                     ************       Make Bar Plots     ************
    """
    norm_transcript_counts = pd.read_csv("normalized_transcript_counts_x_sample.csv")
    if 705 in set(norm_transcript_counts['sample']):
        norm_transcript_counts['sample_new'] = norm_transcript_counts['sample']
    else:
        norm_transcript_counts['sample_new'] = norm_transcript_counts['sample'].apply(lambda sample: sample.replace('.ccs_polish-flnc.fa.algn.srt','')).map(sample_rename_map)
    norm_transcript_counts.drop(columns=['sample'],inplace=True)
    norm_transcript_counts.rename(columns={"sample_new": "sample", "category": "category",'counts':'counts'},inplace=True)
    # get renamed samples
    samples = set(norm_transcript_counts['sample'])
    #print(samples)
    #exit()
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
    
    """
                     ************       Make Scatter Plots     ************
    """    
    # get read positioning data and do some transformations
    norm_read_alignment_and_cats_x_sample = pd.read_csv('normalized_read_categorization_x_sample.csv')
    if 705 in set(norm_read_alignment_and_cats_x_sample['sample']):
        IS_CHIMP_SAMPLE = False
        norm_read_alignment_and_cats_x_sample['sample_new'] = norm_read_alignment_and_cats_x_sample['sample']
    else:
        IS_CHIMP_SAMPLE = True
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
    for alpha in [0.4]:
        for sample in samples:
            df_in = norm_read_alignment_and_cats_x_sample[ norm_read_alignment_and_cats_x_sample['sample'] == sample ]
            scatter(df_in,'transcript_category_color',sample,alpha,IS_CHIMP_SAMPLE)
        # one full scatter for all
        if 705 in set(norm_read_alignment_and_cats_x_sample['sample']):
            sample_comparison_name = "All_Samples_Homo_sapiens"
        else:
            sample_comparison_name = "All_Samples_Chimp"
        scatter(norm_read_alignment_and_cats_x_sample,'transcript_category_color',sample_comparison_name,alpha,IS_CHIMP_SAMPLE)

if __name__ == '__main__':
    replot_analysis()