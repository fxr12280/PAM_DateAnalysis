#高效分析CRISPR/Cas的PAM偏好性的程序，其可以直接使用illumina平台的二代测序的下机数据进行直接分析。
#整个分析过程包括数据下机后的拼接、QC质检、根据barcode进行分类，计数、拟合以及绘制热图。
#本项目提供上述的完整代码以及所需要的tempalte.csv文件的示例，以及完整的代码以及相关的前置库
import glob
import gzip
import itertools
import os
import sys
import warnings
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.optimize import curve_fit
from scipy.stats import linregress
from scipy.stats import pearsonr
from scipy.stats import skew
import subprocess
import csv
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import Counter
def merget(input_forward,input_reverse,output_merged):
    vsearch_command = f"vsearch --fastq_mergepairs {input_forward} --reverse {input_reverse} --fastqout {output_merged}"
    subprocess.run(vsearch_command, shell=True)
def QC_merger(input_fastq,output_dir,output_fastq):
    subprocess.run(["mkdir", "-p", output_dir])
    subprocess.run(["fastqc", input_fastq, "-o", output_dir])
    subprocess.run(["trimmomatic", "SE", "-phred33", input_fastq, output_fastq, "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15"])
class BarcodeClassifier:
    def __init__(self, input_fastq, csv_file, output_directory):
        self.input_fastq = input_fastq
        self.csv_file = csv_file
        self.output_directory = output_directory
        self.barcode_pairs = self.read_barcodes_from_csv() 
        self.barcode_handles = {}  # 

    def read_barcodes_from_csv(self):
        barcode_pairs = {}
        with open(self.csv_file, 'r') as csvfile:
            csvreader = csv.reader(csvfile)
            next(csvreader)  
            for row in csvreader:
                barcode1, barcode2 = row[1].strip(), row[2].strip()
                name = row[3].strip()  # 这里的 name 作为文件名
                barcode_pairs[(barcode1, barcode2)] = name
        return barcode_pairs

    def correct_sequence(self, sequence):
        """计算互补反向序列"""
        base_F = "ATCG"
        base_R = "TAGC"
        complement = {f: r for f, r in zip(base_F, base_R)}
        return ''.join(complement.get(base, base) for base in reversed(sequence))

    def classify_by_barcodes(self):
        """根据 barcode 分类 fastq 文件"""
        unmatched_handle = open(os.path.join(self.output_directory, "unmatched_output.fastq"), "w")

        # 创建输出文件句柄
        for (barcode1, barcode2), filename in self.barcode_pairs.items():
            file_path = os.path.join(self.output_directory, filename)
            self.barcode_handles[(barcode1, barcode2)] = open(file_path, "w")

        # 读取 fastq 文件并进行匹配
        with open(self.input_fastq, "r") as handle:
            for record in SeqIO.parse(handle, "fastq"):
                found_match = False

                # 直接匹配
                for (barcode1, barcode2), filename in self.barcode_pairs.items():
                    if str(record.seq).startswith(barcode1) and str(record.seq).endswith(barcode2):
                        SeqIO.write(record, self.barcode_handles[(barcode1, barcode2)], "fastq")
                        found_match = True
                        break

                # 尝试匹配互补反向序列
                if not found_match:
                    corrected_seq = self.correct_sequence(str(record.seq))
                    for (barcode1, barcode2), filename in self.barcode_pairs.items():
                        if corrected_seq.startswith(barcode1) and corrected_seq.endswith(barcode2):
                            corrected_record = SeqRecord(
                                Seq(corrected_seq),
                                id=record.id,
                                description=record.description,
                                letter_annotations={"phred_quality": record.letter_annotations['phred_quality']}
                            )
                            SeqIO.write(corrected_record, self.barcode_handles[(barcode1, barcode2)], "fastq")
                            found_match = True
                            break

                if not found_match:
                    SeqIO.write(record, unmatched_handle, "fastq")

        for handle in self.barcode_handles.values():
            handle.close()
        unmatched_handle.close()
def fastqcount(barcode_csv ,fastq_dir,  pam_len , spacer, sample_num ,timepoints , pam_orientation):
    try:
        variant_ids = pd.read_csv(barcode_csv)
    except:
        raise Exception("barcode 缺失")
    fastqs = glob.glob(fastq_dir + "/**/*R1*.fastq" ,recursive=True )
    P5_sample_BCs = variant_ids['barcode1'].tolist()
    P5_sample_BC_len = len(P5_sample_BCs[0])
    P7_sample_BCs = variant_ids['barcode2'].tolist()
    P7_sample_BC_len = len(P7_sample_BCs[0])
    nt_complement = dict({'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', '_': '_', '-': '-'})
    nucleotides = ['A', 'T', 'C', 'G']
    total_pam_space = [''.join(p) for p in itertools.product(nucleotides, repeat=pam_len)]
    variant_dict = {}
    for index, row in variant_ids.iterrows():
        variant_dict[str(row['barcode1']) + '_' + str(row['barcode2'])] = row['sample']
    store_all_data = {}
    norm_counts_scale = {}
    for sample in variant_ids['sample'][:sample_num]:
        store_all_data[sample] = {x: [0] * (len(timepoints) - 1) for x in total_pam_space}
    timepoint_fastq = { 'timepoint1_S1': 0, 'timepoint1_S2': 1, 'timepoint1_S3': 2 }
    for fastq in fastqs:
        fastq_name = fastq.split("/")[-1]
        fastq_name = fastq_name.split('_L00')[0]
        timepoint = timepoint_fastq[fastq_name]
        total_reads = 0
        counted_reads = 0
        for seq_record in SeqIO.parse(fastq, "fastq"):
            read_sequence = seq_record.seq
            P5_sample_BC = read_sequence[:P5_sample_BC_len]
            P7_sample_BC = read_sequence[-P7_sample_BC_len:]
            if spacer not in read_sequence:
                continue
            else:
                spacer_loc = read_sequence.find(spacer)
                barcode_pair = P5_sample_BC + '_' + P7_sample_BC
                if barcode_pair in variant_dict.keys():
                    if pam_orientation == "three_prime":
                        PAM = read_sequence[spacer_loc+len(spacer):spacer_loc+len(spacer)+pam_len]
                        try:
                            store_all_data[variant_dict[barcode_pair]][PAM][timepoint] += 1
                            counted_reads += 1
                        except:
                            pass
                    elif pam_orientation == "five_prime":
                        PAM = read_sequence[spacer_loc - pam_len:spacer_loc]
                        try:
                            store_all_data[variant_dict[barcode_pair]][PAM][timepoint] += 1
                            counted_reads += 1
                        except:
                            pass
            total_reads+=1                
        print(counted_reads)
    output_file = os.path.join(fastq_dir, "PAMDA_1_raw_counts.csv")
    with open(output_file, mode='w' ,newline = "") as f_out:
        writer = csv.writer(f_out)
        header = ['Sample', 'PAM'] + [f"Raw_Counts_{x}" for x in range(1, len(timepoints))]
        writer.writerow(header)

        for sample, pam_counts in store_all_data.items():
            for pam, counts in pam_counts.items():
                writer.writerow([sample, pam] + counts)
def fastqcount(barcode_csv ,fastq_dir,  pam_len , spacer, sample_num ,timepoints , pam_orientation):
    try:
        variant_ids = pd.read_csv(barcode_csv)
    except:
        raise Exception("barcode 缺失")
    fastqs = glob.glob(fastq_dir + "/**/*R1*.fastq" ,recursive=True )
    P5_sample_BCs = variant_ids['barcode1'].tolist()
    P5_sample_BC_len = len(P5_sample_BCs[0])
    P7_sample_BCs = variant_ids['barcode2'].tolist()
    P7_sample_BC_len = len(P7_sample_BCs[0])
    nt_complement = dict({'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', '_': '_', '-': '-'})
    nucleotides = ['A', 'T', 'C', 'G']
    total_pam_space = [''.join(p) for p in itertools.product(nucleotides, repeat=pam_len)]
    variant_dict = {}
    for index, row in variant_ids.iterrows():
        variant_dict[str(row['barcode1']) + '_' + str(row['barcode2'])] = row['sample']
    store_all_data = {}
    norm_counts_scale = {}
    for sample in variant_ids['sample'][:sample_num]:
        store_all_data[sample] = {x: [0] * (len(timepoints) - 1) for x in total_pam_space}
    timepoint_fastq = { 'timepoint1_S1': 0, 'timepoint1_S2': 1, 'timepoint1_S3': 2 }
    for fastq in fastqs:
        fastq_name = fastq.split("/")[-1]
        fastq_name = fastq_name.split('_L00')[0]
        timepoint = timepoint_fastq[fastq_name]
        total_reads = 0
        counted_reads = 0
        for seq_record in SeqIO.parse(fastq, "fastq"):
            read_sequence = seq_record.seq
            P5_sample_BC = read_sequence[:P5_sample_BC_len]
            P7_sample_BC = read_sequence[-P7_sample_BC_len:]
            if spacer not in read_sequence:
                continue
            else:
                spacer_loc = read_sequence.find(spacer)
                barcode_pair = P5_sample_BC + '_' + P7_sample_BC
                if barcode_pair in variant_dict.keys():
                    if pam_orientation == "three_prime":
                        PAM = read_sequence[spacer_loc+len(spacer):spacer_loc+len(spacer)+pam_len]
                        try:
                            store_all_data[variant_dict[barcode_pair]][PAM][timepoint] += 1
                            counted_reads += 1
                        except:
                            pass
                    elif pam_orientation == "five_prime":
                        PAM = read_sequence[spacer_loc - pam_len:spacer_loc]
                        try:
                            store_all_data[variant_dict[barcode_pair]][PAM][timepoint] += 1
                            counted_reads += 1
                        except:
                            pass
            total_reads+=1                
        print(counted_reads)
    output_file = os.path.join(fastq_dir, "PAMDA_1_raw_counts.csv")
    with open(output_file, mode='w' ,newline = "") as f_out:
        writer = csv.writer(f_out)
        header = ['Sample', 'PAM'] + [f"Raw_Counts_{x}" for x in range(1, len(timepoints))]
        writer.writerow(header)

        for sample, pam_counts in store_all_data.items():
            for pam, counts in pam_counts.items():
                writer.writerow([sample, pam] + counts)
def rawcount2normcount(pam_len,input_csv , timepoints , control_sample ,fastq_dir):
    control_sample_timepoint_fastq  = 1
    top_n = 5
    nucleotides = ['A', 'T', 'C', 'G']
    total_pam_space = [''.join(p) for p in itertools.product(nucleotides, repeat=pam_len)]
                       
    df_input = pd.read_csv(input_csv)

    count_columns = df_input.columns.values[2:]
    df_list = []
    for sample in df_input["Sample"].unique().tolist():
        temp_df = df_input[df_input['Sample'] == sample].sort_values(by=['PAM'])
        df_list.append(temp_df)
    df = pd.concat(df_list)
    for i in range(1,len(timepoints)):
        df['Norm_Counts_' + str(i)] = df['Raw_Counts_' + str(i)] / df.groupby(['Sample'])['Raw_Counts_' + str(i)].transform("sum").replace(0, np.nan)
    
    control_dict = {}
    for PAM in total_pam_space:
         control_dict[PAM] = 0
    for index, row in df.iterrows():
        if row['Sample'] == control_sample:
            control_dict[row['PAM']] = row['Norm_Counts_' +str(control_sample_timepoint_fastq)]
    norm_counts_0 = []
    for index, row in df.iterrows():
        norm_counts_0.append(control_dict[row['PAM']])
    df['Norm_Counts_0'] = norm_counts_0
    df = df[df['Sample'] != control_sample]
    sample_last = None
    uptrends = {}
    x = range(len(timepoints))
    for index, row in df.iterrows():
        sample_current = row['Sample']
        y = [row['Norm_Counts_' + str(i)] for i in range(len(timepoints))]
        slope = linregress(x, y)
        if sample_current in uptrends:
            uptrends[sample_current].append([slope[0], y])
        else:
            uptrends[sample_current] = [[slope[0], y]]
        sample_last = sample_current
    uptrend_corrections = {}
    for u in uptrends:
        uptrends[u] = sorted(uptrends[u])
        top_n_entries = [x[1] for x in uptrends[u][-top_n:]]
        top_n_entries_reformat = map(list, zip(*top_n_entries))
        top_n_entries_median = [np.median(x) for x in top_n_entries_reformat]
        uptrend_corrections[u] = top_n_entries_median
    sample_last = None
    for index, row in df.iterrows():
        sample_current = row['Sample']

        for i in range(len(timepoints)):
            df.loc[index, 'Norm_Counts_' + str(i)] = row['Norm_Counts_' + str(i)] / \
                                                     uptrend_corrections[sample_current][i]
            df.loc[index, 'Norm_Counts_' + str(i)] = row['Norm_Counts_' + str(i)] / \
                                                     row['Norm_Counts_0']
        sample_last = sample_current
    output_file = os.path.join(fastq_dir, "PAMDA_2_norm_counts.csv")
    df.to_csv(output_file ,index=False )    
 def normcount2rate(fastq_dir,input_csv , timepoints):
    read_sum_minimum=10
    tps_sum_minimum=2
    init_rate_est=[0.0001, 0.001, 0.01]
    df = pd.read_csv(input_csv)
    def func(x , a, b):
        return a*np.exp(-b*x)
    ks = []
    timepoint_indices = []
    for i in range(len(timepoints)):
        timepoint_indices.append([i , timepoints[i]])
    previous_row = 'n/a'
    for index, row in df.iterrows():

        current_row = row['Sample']

        tps = [x[1] for x in timepoint_indices]
        obs_raw = [0.00000000001] + [row['Raw_Counts_' + str(x[0])]
                                     for x in timepoint_indices if str(x[0]) != '0']
        obs_norm = [row['Norm_Counts_' + str(x[0])] for x in timepoint_indices]

        zero_indices = [i for i, e in enumerate(obs_raw) if e == 0]
        for zero_index in sorted(zero_indices)[::-1]:
            del tps[zero_index]
            del obs_norm[zero_index]
            del obs_raw[zero_index]

        if sum(obs_raw) >= read_sum_minimum and len(obs_norm) >= tps_sum_minimum:

        
            min_search = []

        
            for j in init_rate_est:
                p0 = [1.0, j]
        
                try:
                    popt, pcov = curve_fit(func, tps, obs_norm, p0=p0)
                    pred = [func(x, popt[0], popt[1]) for x in tps]
        
                    error = sum([(x[0] - x[1]) ** 2 for x in zip(pred, obs_norm)])

                    min_search.append([error, list(popt)])
                except:
                    continue
            if len(min_search) != 0:
                ks.append(sorted(min_search)[0][1][1])
            else:
                ks.append('NaN')

        else:
            ks.append('NaN')
    min_k = min([100 if ((isinstance(x, float) and x <= 0) or x == 'NaN') else x for x in ks])
    df['Rate_Constant_k'] = [x if ((isinstance(x, float) and x > 0) or x == 'NaN') else min_k for x in ks]
    output_file = os.path.join(fastq_dir, "PAMDA_3_rates.csv")
    df.to_csv(output_file ,index=False )
def rate2heatmap(fastq_dir,barcode_csv,pam_len,input_csv):

    pam1_nucleotide_rank={1: 'A', 2: 'C', 3: 'G', 4: 'T'}
    pam2_nucleotide_rank={1: 'A', 2: 'C', 3: 'G', 4: 'T'}
    heatmap_fixed_min=False
    heatmap_fixed_max=False
    plt.switch_backend('agg')

    variant_ids = pd.read_csv(barcode_csv)  
    variants = variant_ids['sample'].tolist()
    variant_name_dict = {}
    for index, row in variant_ids.iterrows():
        variant_name_dict[row['sample']] = row['description']

    def translate_pam(pam_as_numbers, tranlsate_dict):
        pam_as_bases = ''
        for n in str(pam_as_numbers):
            pam_as_bases += tranlsate_dict[int(n)]
        return pam_as_bases

    # Reformatting input dataframe
    df_input = pd.read_csv(input_csv)

    split_pam_index = np.floor_divide(pam_len, 2)
    pam1_index_rank = [x for x in range(0, split_pam_index)][::-1]
    pam2_index_rank = [x for x in range(split_pam_index, pam_len)]

    df_input['PAM_pt1'] = [x[:split_pam_index] for x in df_input['PAM'].tolist()]
    df_input['PAM_pt2'] = [x[split_pam_index:] for x in df_input['PAM'].tolist()]


    numbers = ['1', '2', '3', '4']
    pam_space = [''.join(p) for p in itertools.product(numbers, repeat=pam_len)]

    # Sort PAMs according to rules
    pam1_ids = []
    pam2_ids = []
    for pam in pam_space:
        pam1_ids.append(int(pam[:split_pam_index]))
        pam2_ids.append(int(pam[split_pam_index:]))

    pam1_ids = sorted([int(x) for x in set(pam1_ids)])
    pam2_ids = sorted([int(x) for x in set(pam2_ids)])

    columns = []
    indices = []
    tmp = ''
    for pam in pam1_ids:
        tmp_dict = {i: j for i, j in zip(pam1_index_rank, list(str(pam)))}
        for i in sorted(tmp_dict):
            tmp += pam1_nucleotide_rank[int(tmp_dict[i])]
        indices.append(tmp)
        tmp = ''

    tmp = ''
    for pam in pam2_ids:
        tmp_dict = {i: j for i, j in zip(pam2_index_rank, list(str(pam)))}
        for i in sorted(tmp_dict):
            tmp += pam2_nucleotide_rank[int(tmp_dict[i])]
        columns.append(tmp)
        tmp = ''

    def save_heatmap():
        """
        save heatmaps
        """

        plt.savefig(fastq_dir+ '/PAMDA_HEATMAP_%s_%s.pdf' %(variant, variant_name_dict[variant]))
        df_output.to_csv(fastq_dir+ '/PAMDA_%s_%s.csv' %(variant, variant_name_dict[variant]))

    def plot_nice_heatmap(df_output, spacer=None):
        """
        plot pretty heatmaps with color-formatted axes labels (only for 2-5 nt PAMs)
        """
        heatmap_min = None
        heatmap_max = None
        if heatmap_fixed_min:
            heatmap_min = heatmap_fixed_min
        if heatmap_fixed_max:
            heatmap_max = heatmap_fixed_max
        cbar_label = 'Log10(rate)'
        # heatmap plotting
        axes = [4 * (pam_len - split_pam_index), 4 * (split_pam_index)]
        fig, ax1 = plt.subplots(1, figsize=(axes[0], axes[1]))
        sns.heatmap(df_output,
                    ax=ax1,
                    vmin=heatmap_min,
                    vmax=heatmap_max,
                    square=True,
                    cmap='Blues',
                    cbar=True,
                    cbar_kws={'shrink': axes[1] / axes[0] / 2,
                              'label': cbar_label,
                              'aspect': 8},
                    linewidth=0.2,
                    linecolor="White",
                    xticklabels=False,
                    yticklabels=False)
        heatmap_size = fig.get_size_inches()

        # format the axis labels

        # colors of the bases of the axis labels
        colors = {'A': '#feb2b1', 'C': '#14c7fe', 'T': '#14f485', 'G': '#f8ffa3'}
        # scaling dict manually scales the colored axis labels
        # dict structure: {pam_len:{split_index:[x_width,x_height,y_width,y_height]}}
        scaling_dict = {2: {1: [1, 3.5, 4, 1.07]},
                        3: {1: [1, 1.7, 1, 2.15], 2: [1, 2, 4, 3.53]},
                        4: {1: [1, 0.7, 0.3, 5.5], 2: [1, 1.7, 1, 4.3], 3: [1, 1, 4, 14.1]},
                        5: {1: [1, 0.25, 0.08, 17], 2: [1, 0.6, 0.3, 11.5],
                            3: [1, 1, 1, 14.15], 4: [1, 0.3, 3, 56.5]}}
        x_text = [[columns[n][m] for n in range(len(columns))] for m in range(len(columns[0]))]
        x_text = x_text[::-1]
        x_color = [[colors[x_text[n][m]] for m in range(len(x_text[0]))] for n in range(len(x_text))]
        xtable = ax1.table(cellText=x_text, cellColours=x_color,
                           cellLoc='center', loc='top')
        for key, cell in xtable.get_celld().items():
            cell.set_linewidth(0)
        y_text = [[indices[n][m] for m in range(len(indices[0]))] for n in range(len(indices))]
        y_color = [[colors[y_text[n][m]] for m in range(len(y_text[0]))] for n in range(len(y_text))]
        y_widths = [0.06 for n in enumerate(y_text)]
        ytable = ax1.table(cellText=y_text, cellColours=y_color, colWidths=y_widths,
                           cellLoc='center', loc='left')
        for key, cell in ytable.get_celld().items():
            cell.set_linewidth(0)
        xtable.set_fontsize(8)
        ytable.set_fontsize(8)
        xtable.scale(scaling_dict[pam_len][split_pam_index][0], scaling_dict[pam_len][split_pam_index][1])
        ytable.scale(scaling_dict[pam_len][split_pam_index][2],heatmap_size[1] / scaling_dict[pam_len][split_pam_index][3])
        plt.tight_layout(pad=6)

        # save the heatmap
        save_heatmap()
        plt.close()

    new_columns = []
    for variant in df_input['Sample'].unique().tolist():
        for column in columns:
            new_columns.append(str(column) + '-' + str(variant))
        df_output = pd.DataFrame(columns=new_columns, index=indices)

    
        for row in df_output.index:
            for column in df_output.columns:
                pam2 = str(column.split('-')[0].strip('\n'))
                sample = str(column.split('-')[1].strip('\n'))
                rate_avg = 0

                rate_avg += float(df_input['Rate_Constant_k'][(df_input['PAM_pt1'] == str(row)) &
                                                                (df_input['PAM_pt2'] == pam2) &
                                                                (df_input['Sample'] == sample) ])
                rate_avg = rate_avg / len(spacer)
                df_output.loc[row, column] = np.log10(rate_avg)

        df_output = df_output[df_output.columns].astype(float)
        plot_nice_heatmap(df_output)
        new_columns = []
def PAM_DateAnalysis(fastq_dir , pam_len, spacer , sample_num ,timepoints , pam_orientation , control_sample):
    input_forward = fastq_dir + "/F.fq"
    input_reverse = fastq_dir + "/R.fq"
    output_merged = fastq_dir + "/merged.fasta"
    merget(input_forward,input_reverse,output_merged)
    input_fastq = output_merged
    output_dir = fastq_dir + "/fastqc_output"
    output_fastq = fastq_dir + "/output_trimmed.fastq"
    QC_merger(input_fastq,output_dir,output_fastq)
    input_fastq = output_fastq
    barcode_csv = fastq_dir + "/tempalte_CSV.csv"
    output_directory = fastq_dir
    classifier = BarcodeClassifier(input_fastq, barcode_csv, output_directory)
    classifier.classify_by_barcodes()
    fastqcount(barcode_csv ,fastq_dir,  pam_len , spacer, sample_num ,timepoints , pam_orientation)
    fastqcount_csv = fastq_dir + "/PAMDA_1_raw_counts.csv" 
    rawcount2normcount(pam_len,fastqcount_csv , timepoints , control_sample ,fastq_dir)
    rawcount2normcount_csv = fastq_dir + "/PAMDA_2_norm_counts.csv"
    normcount2rate(fastq_dir, rawcount2normcount_csv , timepoints)
    rawcount2normcount_csv = fastq_dir + "/PAMDA_3_rates.csv"
    rate2heatmap(fastq_dir,barcode_csv,pam_len,rawcount2normcount_csv)
    
    
    
                             
    
