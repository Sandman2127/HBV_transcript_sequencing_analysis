#!/opt/anaconda/envs/py3/bin/python

import argparse
import pandas as pd

def combine_categories(df):
    # merge HBsAg_L1 & HBsAg_L2 , HBsAg_M1 & HBsAg_M2 , HBsAg_S1 & HBsAg_S2. HBsAg_L1_iDNA & HBsAg_L2_iDNA, HBsAg_M1_iDNA & HBsAg_M2_iDNA, HBsAg_S1_iDNA & HBsAg_S2_iDNA
    df.loc['HBsAg_L'] = df.loc['HBsAg_L1'] + df.loc['HBsAg_L2']
    df.loc['HBsAg_M'] = df.loc['HBsAg_M1'] + df.loc['HBsAg_M2']
    df.loc['HBsAg_S'] = df.loc['HBsAg_S1'] + df.loc['HBsAg_S2']
    df.loc['HBsAg_L_iDNA'] = df.loc['HBsAg_L1_iDNA'] + df.loc['HBsAg_L2_iDNA']
    df.loc['HBsAg_M_iDNA'] = df.loc['HBsAg_M1_iDNA'] + df.loc['HBsAg_M2_iDNA']
    df.loc['HBsAg_S_iDNA'] = df.loc['HBsAg_S1_iDNA'] + df.loc['HBsAg_S2_iDNA']
    df.loc['X_trunc1'] = df.loc['X_trunc1'] + df.loc['X1_trunc1']
    df.loc['X_trunc2'] = df.loc['X_trunc2'] + df.loc['X1_trunc2']
    df.loc['X_trunc3'] = df.loc['X_trunc3'] + df.loc['X1_trunc3']
    df.loc['X_long'] = df.loc['X1_long'] + df.loc['X2_long']
    # drop all rows that were merged
    df.drop([ 
            'HBsAg_L1',
            'HBsAg_L2',
            'HBsAg_M1',
            'HBsAg_M2',
            'HBsAg_S1',
            'HBsAg_S2',
            'HBsAg_L1_iDNA',
            'HBsAg_L2_iDNA',
            'HBsAg_M1_iDNA',
            'HBsAg_M2_iDNA',
            'HBsAg_S1_iDNA',
            'HBsAg_S2_iDNA',
            'X1_long',
            'X2_long',
            'X_splice1',
            'X_splice2',
            'X_splice3',
            'X_splice4',
            'X1_trunc1',
            'X1_trunc2',
            'X1_trunc3',
             ], inplace=True)
    return df

def main():
    parser = argparse.ArgumentParser(description='Calculate normalized transcript counts for each group')
    parser.add_argument('input_csv', help='File containing sample group information')
    args = parser.parse_args()

    HBe_neg = {
                '88A010_57.ccs_polish-flnc.fa.algn.srt',
                '95A010_d57.ccs_polish-flnc.fa.algn.srt',
                '89A008_Day_379.ccs_polish-flnc.fa.algn.srt'
            }
    
    HBe_pos = {
                '89A008_HC.ccs_polish-flnc.fa.algn.srt',
                '89A008_Day_141.ccs_polish-flnc.fa.algn.srt',
                'A2A004_HC.ccs_polish-flnc.fa.algn.srt',
                'A2A004_d85.ccs_polish-flnc.fa.algn.srt',
                'A2A004_d323.ccs_polish-flnc.fa.algn.srt',
                'A2A004_Day_351.ccs_polish-flnc.fa.algn.srt',
                'A3A006_HC.ccs_polish-flnc.fa.algn.srt',
                'A3A006_d141.ccs_polish-flnc.fa.algn.srt',
                'A3A006_Day_351.ccs_polish-flnc.fa.algn.srt',
                'A4A014_HC.ccs_polish-flnc.fa.algn.srt',
               }
    
    category_order = [ 
                    'HBeAg_wt',
                    'HBeAg_ORF',
                    'HBeAg_sp6',
                    'HBeAg_sp11',
                    'HBeAg_sp1_Cys-',
                    'HBeAg-S_sp5',
                    'HBeAg-S_sp9',
                    'POL',
                    'Core_wt',
                    'Core_ORF',
                    'Core_sp6',
                    'Core_sp11',
                    'Core-sp1_Cys-',
                    'Core-S_sp5',
                    'Core-S_sp9',
                    'S-Core',
                    'HBsAg_L',
                    'HBsAg_L_iDNA',
                    'HBsAg_M',
                    'HBsAg_M_iDNA',
                    'HBsAg_S',
                    'HBsAg_S_iDNA',
                    'HBsAg_2X_1',
                    'HBsAg_2X_2',
                    'HBsAg_2X_3',
                    'HBsAg_2X_4',
                    'X1_canonical',
                    'X_long',
                    'X_trunc1',
                    'X_trunc2',
                    'X_trunc3',
                    'X_2X',
                    'Other',
                ]

    csv = pd.read_csv(args.input_csv, sep=',')
    # select HBe positive and negative samples:
    HBe_pos_df = csv[csv['sample'].isin(HBe_pos)]
    HBe_neg_df = csv[csv['sample'].isin(HBe_neg)]
    # lets check if all samples are found in the input csv
    if len(HBe_pos_df['sample'].unique()) != len(HBe_pos):
        print('Error: not all HBe + samples found in input csv')
        exit(1)
    if len(HBe_neg_df['sample'].unique()) != len(HBe_neg):
        print('Error: not all HBe - samples found in input csv')
        exit(1)
    # group:
    HBe_pos_gb = HBe_pos_df.groupby('category')['count'].sum()
    HBe_neg_gb = HBe_neg_df.groupby('category')['count'].sum()
    # convert to df:
    HBe_pos_gb_df = pd.DataFrame(HBe_pos_gb)
    HBe_neg_gb_df = pd.DataFrame(HBe_neg_gb)
    # merge all 1 & 2 categories and drop the original columns
    HBe_pos_gb_df = combine_categories(HBe_pos_gb_df)
    HBe_neg_gb_df = combine_categories(HBe_neg_gb_df)
    # get total sums
    HBe_pos_gb_total = HBe_pos_gb_df['count'].sum()
    HBe_neg_gb_total = HBe_neg_gb_df['count'].sum()
    # calc % for each category in each group
    HBe_pos_gb_df['% total'] = round((HBe_pos_gb_df['count'] / HBe_pos_gb_total) * 100,2)
    HBe_neg_gb_df['% total'] = round((HBe_neg_gb_df['count'] / HBe_neg_gb_total) * 100,2)
    # need to reorder the output order via a list
    HBe_pos_gb_df = HBe_pos_gb_df.reindex(category_order)
    HBe_neg_gb_df = HBe_neg_gb_df.reindex(category_order)
    # rename column to manuscripts name 
    HBe_pos_gb_df.rename(mapper={'X1_canonical':'X_canonical'},inplace=True,axis=0)
    HBe_neg_gb_df.rename(mapper={'X1_canonical':'X_canonical'},inplace=True,axis=0)
    # fill in any NaN values with 0, only occurs of HBsAg_2X_1
    HBe_pos_gb_df.fillna(0, inplace=True)
    HBe_neg_gb_df.fillna(0, inplace=True)
    # print output dfs for each group:
    print()
    print('HBe+',HBe_pos_gb_df)
    print(HBe_pos_gb_df['% total'].sum())
    print(HBe_pos_gb_df['count'].sum())
    print()
    print('HBe-',HBe_neg_gb_df)
    print(HBe_neg_gb_df['% total'].sum())
    print(HBe_neg_gb_df['count'].sum())
    # write each to an independent csv with name by group
    HBe_pos_gb_df.to_csv('HBe_+_normalized_category_counts.csv')
    HBe_neg_gb_df.to_csv('HBe_-_normalized_category_counts.csv')

if __name__ == "__main__":
    main()

# sample,category,count
# 89A008_Day_379.ccs_polish-flnc.fa.algn.srt,HBsAg_L1,0
# 89A008_Day_379.ccs_polish-flnc.fa.algn.srt,Core_wt,0
# 89A008_Day_379.ccs_polish-flnc.fa.algn.srt,HBeAg_ORF,0
# 89A008_Day_379.ccs_polish-flnc.fa.algn.srt,X_splice3,0 
