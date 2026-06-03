# 前ver.: py_4GTF_HSGRCh38_20250812_1.py


import argparse
import pandas as pd
import numpy as np
import csv
import itertools
from multiprocessing import Pool
import multiprocessing as mp
import glob
import sys
import os
import time
import datetime
import re
#import matplotlib.pyplot as plt

start_time = time.time()
dt_4GTF_start = datetime.datetime.now()
print( "dt_4GTF_start", dt_4GTF_start)


aAP = argparse.ArgumentParser()

aAP.add_argument(
    "-b", "--blast_result", metavar="BLASTRESULT_TXT", required=True,
    help="the input file - BLAST result (required)"
)
aAP.add_argument(
    "-c", "--chr_converter", metavar="LABEL_TSV", required=True,
    help="a converter of Acc. No. <-> chr No. (required)"
)
aAP.add_argument(
    "-s", "--species_label", metavar=str, required=True,
    help="species label (HS, PTR, GGO, ......) (required)"
)
aAP.add_argument(
    "-n", "--native_align", metavar="TSV", default=None,
    help="optional; cross-species analysis, the native locus exon tsv required"
)
aAP.add_argument(
    "-t", "--TEdup_threshold", type=int, default=100, metavar="N",
    help="optional; TE duplication threshold (default: 100)"
)

p_args = aAP.parse_args()
p_args_br = p_args.blast_result
p_args_cc = p_args.chr_converter
p_args_sl = p_args.species_label
p_args_na = p_args.native_align# if p_args.native_align else None
p_args_Tt = p_args.TEdup_threshold


#folder_path = sys.argv[1]  # 引数でフォルダパスを指定
#file_list = glob.glob(f"{folder_path}/*605242*.txt")  # .txtファイルをリスト化
#file_list = glob.glob(f"{folder_path}/*469061*.txt")  # .txtファイルをリスト化
#file_list = glob.glob(f"{folder_path}/*612828*.txt")  # .txtファイルをリスト化
#file_list = glob.glob(f"{folder_path}/*525723*.txt")  # .txtファイルをリスト化
#file_list = glob.glob(f"{folder_path}/*537473*.txt")  # .txtファイルをリスト化
#file_list = glob.glob(f"{folder_path}/*611957*.txt")  # .txtファイルをリスト化
#file_list = glob.glob(f"{folder_path}/*637485*.txt")  # .txtファイルをリスト化
#file_list = glob.glob(f"{folder_path}/*256646*.txt")  # .txtファイルをリスト化
#file_list = glob.glob(f"{folder_path}/*.txt")  # .txtファイルをリスト化

## ファイルを順番に処理
#for p_args_br in file_list:
#    # ファイル名のベース部分だけを取り出す（パスを除く）
#    filename_only = os.path.basename(p_args_br)
#
#    # 出力ディレクトリをenst_folderに指定
#    output_dir = "enst_folder"

###args = sys.argv

#with open('x_onlyperfect_testgenes.tsv', 'r') as tl:
#    tl_list = tl.read().split("\n")
#    print(tl_list)
#    print(tl_list[:-1])



#def forGTF(dups_txt_file): #p_args_br

    # データの読み込み
    #dups_txt_file = "merge_blastresult_00000602_210406.txt" #[p_args_br]
#df = pd.read_table(dups_txt_file, #p_args_br

dtype_map = {
  "query acc.ver": str,
  "subject acc.ver": str,
  "% identity": float,
  "alignment length": int,
  "mismatches": int,
  "gap opens": int,
  "q_head": int,
  "q_tail": int,
  "s_head": int,
  "s_tail": int,
  "evalue": str,
  "bit score": float,
}

df = pd.read_table(p_args_br,
        header     = None,
        index_col  = False,
        names      = ["query acc.ver", "subject acc.ver", "% identity",
                      "alignment length", "mismatches", "gap opens",
                      "q_head", "q_tail", "s_head", "s_tail",
                      "evalue", "bit score"],
        dtype      = dtype_map,
        low_memory = False,
        converters = {"subject acc.ver":str})
#print("df",df)

if df.empty:                      # ← DataFrameが空（行数0）ならTrue
    print("empty DataFrame - de novo?")
    query_type = "unknown - de novo?"
    print(p_args_br.split("."))
    pd.DataFrame([{"query_tra_ID": p_args_br.split(".")[0],#sss.iloc[0,0],
                       "query_gene_name": "(unable_to_retrieve)",
                       "query_type" : query_type,
                       f"{p_args_sl}_total_CN": 0,
                       f"{p_args_sl}_SD_CN"   : 0,
                       f"{p_args_sl}_retro_CN": 0,
                       f"({p_args_sl}_single)": 0,
                       f"({p_args_sl}_TE)"    : 0,
                       "repeat?": "No"}]).to_csv("statistics_{}_{}.tsv".format(p_args_br, "(unable_to_retrieve)"), sep="\t", index=False, mode="w")

    raise SystemExit(0)


elapsed_l77 = time.time() - start_time
dt_l77 = datetime.datetime.now()
print(f"DataFrame l-77, 読み込み完了: {elapsed_l77:.2f} 秒, 現在: {dt_l77}")
t_l77 = time.time()

first_field = df.iat[0, 0]  # 先頭行・先頭カラム
ffsplit = first_field.split("___")
ff6split = ffsplit[6].split(":")

#print(ffsplit, ff6split)

if "transcript_biotype:retained_intron" in first_field:
    # 遺伝子名を抽出
    sss = df["query acc.ver"].str.split("___", expand=True)  # 参照: https://note.nkmk.me/python-pandas-split-extract/
    s2split = sss[2].str.split(":", expand=True)
    s2split[[3,4,5]] = s2split[[3,4,5]].astype("int")
    chr_No = s2split[2]
    q_inv = s2split.loc[0,5]
    s6split = sss[6].str.split(":", expand=True)
    pd.DataFrame([{"query_tra_ID": ffsplit[0],
                   "query_gene_name": ff6split[1],
                   "query_type" : "retained_intron",
                   f"{p_args_sl}_total_CN": -10,
                   f"{p_args_sl}_SD_CN"   : -10,
                   f"{p_args_sl}_retro_CN": -10,
                   f"({p_args_sl}_single)": -10,
                   f"({p_args_sl}_TE)"    : -10,
                   "repeat?": "No"}]).to_csv("statistics_{}_{}.tsv".format(p_args_br, s6split.iloc[0,1]), sep="\t", index=False, mode="w")


    print(p_args_br)
    sys.exit("Skipped - This BLAST-result has been generated from retained_intron query.")


# ラベルについて。
# start/end : 単純に番地が若い方の端点がstart
# head/tail : exonの向きを考慮し、先頭塩基の番地をhead
# inversionが-なら、start == tail になるということ。

# 列の追加
df.insert(1,"start",0)
df.insert(2,"end",0)
df.insert(3,"inversion",0)
df.insert(4,"ances",0)
#df.insert(5,"rank",0)
df.insert(5,"q_whole_start",0)
df.insert(6,"q_whole_end",0)
df.insert(7,"repeat",-1)
df.insert(8,"rank",0)
df.insert(9,"from",0)
#df.insert(10,"retro","no")

elapsed_l120 = time.time() - start_time
t_l120 = time.time()
t_l77_120 = t_l120 - t_l77
dt_l120 = datetime.datetime.now()
print(f"l77-120, done: {t_l77_120:.2f} 秒, 4GTF_startから: {elapsed_l120:.2f} 秒, 現在: {dt_l120}")

#print("insert done")

df["name"] = df["query acc.ver"]

dfr = df.reindex(columns=["query acc.ver", "name", "start", "end", "inversion",
                          "ances", "repeat", "rank", "from", "q_whole_start", "q_whole_end",
                          "subject acc.ver", "% identity", "alignment length", "mismatches", "gap opens",
                          "q_head", "q_tail", "s_head", "s_tail", "evalue", "bit score"])

#print("reindex done")

df_cc = pd.read_csv(p_args_cc, header=None, sep=r"\s+", names=["Acc. No.", "chr. No."])
#print(df_cc)

# 個人差を除去
chr_No_list = df_cc.iloc[:, -1].astype(str).tolist()
#print(chr_No_list)

#number_series = ["{}".format(i) for i in range(1, 23)]
#number_series += ["X", "Y"]
#print(number_series)
#print(number_series[0])
#print(number_series[21])
#print(number_series[22])
#print(number_series[23])
#print(type(number_series[0]))
#print(type(number_series[21]))
#print(type(number_series[22]))
#print(type(number_series[23]))

# 変換辞書を作成（Acc. No. -> chr. No.）
acc2chr = df_cc.set_index("Acc. No.")["chr. No."].astype(str).to_dict()
#print(acc2chr)

dfr["subject acc.ver"] = dfr["subject acc.ver"].replace(acc2chr)

#dfr = dfr.replace({"subject acc.ver": {"NC_000001.11": "1", "NC_000002.12": "2", "NC_000003.12": "3", "NC_000004.12": "4",
#    "NC_000005.10": "5", "NC_000006.12": "6", "NC_000007.14": "7", "NC_000008.11": "8", "NC_000009.12": "9", "NC_000010.11": "10", "NC_000011.10": "11",
#    "NC_000012.12": "12", "NC_000013.11": "13", "NC_000014.9": "14", "NC_000015.10": "15", "NC_000016.10": "16", "NC_000017.11": "17", "NC_000018.10": "18",
#    "NC_000019.10": "19", "NC_000020.11": "20", "NC_000021.9": "21", "NC_000022.11": "22", "NC_000023.11": "X", "NC_000024.10": "Y"}})

#elapsed_l152 = time.time() - start_time
#print(f"DataFrame 読み込み完了: {elapsed_l152:.2f} 秒")

allowed_chr = set(map(str, chr_No_list))
df = dfr[dfr["subject acc.ver"].astype(str).isin(allowed_chr)]

#df = dfr[dfr["subject acc.ver"] in number_series] ## 多分以下の処理はこれで代用できると思う。エラー履いたら戻してくれ
#df = dfr[(dfr["subject acc.ver"] == "1")|(dfr["subject acc.ver"] == "2")|
#          (dfr["subject acc.ver"] == "3")|(dfr["subject acc.ver"] == "4")|
#          (dfr["subject acc.ver"] == "5")|(dfr["subject acc.ver"] == "6")|
#          (dfr["subject acc.ver"] == "7")|(dfr["subject acc.ver"] == "8")|
#          (dfr["subject acc.ver"] == "9")|(dfr["subject acc.ver"] == "10")|
#          (dfr["subject acc.ver"] == "11")|(dfr["subject acc.ver"] == "12")|
#          (dfr["subject acc.ver"] == "13")|(dfr["subject acc.ver"] == "14")|
#          (dfr["subject acc.ver"] == "15")|(dfr["subject acc.ver"] == "16")|
#          (dfr["subject acc.ver"] == "17")|(dfr["subject acc.ver"] == "18")|
#          (dfr["subject acc.ver"] == "19")|(dfr["subject acc.ver"] == "20")|
#          (dfr["subject acc.ver"] == "21")|(dfr["subject acc.ver"] == "22")|
#          (dfr["subject acc.ver"] == "X")|(dfr["subject acc.ver"] == "Y")]





df = df.reset_index(drop = True)
#print("individuals_del done")
#print(df)

elapsed_l178 = time.time() - start_time
t_l178 = time.time()
t_l120_178 = t_l178 - t_l120
dt_l178 = datetime.datetime.now()
print(f"l120-178, done: {t_l120_178:.2f} 秒, 4GTF_startから: {elapsed_l178:.2f} 秒, 現在: {dt_l178}")


# 遺伝子名を抽出
sss = df["name"].str.split("___", expand=True)  # 参照: https://note.nkmk.me/python-pandas-split-extract/
#s2split = df["name"].str.split(":", expand=True)
s2split = sss[2].str.split(":", expand=True)
s2split[[3,4,5]] = s2split[[3,4,5]].astype("int")
chr_No = s2split[2]
q_inv = s2split.loc[0,5]
s6split = sss[6].str.split(":", expand=True)
df["name"] = sss.iloc[:,0] + "_" + s6split.iloc[:,1] + "_" ### NOTCH2_ENST00000256646.7

try:
    _qinv_int = int(q_inv)
except Exception:
    _qinv_int = 1
q_inv_sign = "+" if _qinv_int == 1 else "-"


# GTFのinversion/forward/backwardカラムを作成
df["q_whole_start"] = s2split.iloc[:, 3]
df["q_whole_end"] = s2split.iloc[:, 4]
idx_forward = df["s_head"] < df["s_tail"]
idx_backward = df["s_head"] > df["s_tail"]
df.loc[idx_forward,"inversion"] = "+"
df.loc[idx_backward,"inversion"] = "-"
#print(idx_forward)
## 厳密にはidxではないけど該当する行にTrueが入っているからdf[idx_forward]で条件に合致する行だけを選択していることになる
#print("dgfsdfgs")
#print(df[idx_forward])
#print(df.loc[idx_forward,"s_start"])
df.loc[idx_forward,"start"] = df.loc[idx_forward,"s_head"]
df.loc[idx_forward,"end"] = df.loc[idx_forward,"s_tail"]
df.loc[idx_backward,"start"] = df.loc[idx_backward,"s_tail"]
df.loc[idx_backward,"end"] = df.loc[idx_backward,"s_head"]


# 染色体上での並び順、クエリごと、
# の順に並べ変えてインデックス振り直し(て表を表示)
predfsv = df.sort_values(by = ["subject acc.ver", "start"],
                         ascending = [True, True])
dfsv = predfsv.reset_index(drop = True)
print("l-269", dfsv, dfsv[dfsv["ances"]==1])

TE_dfsv = pd.DataFrame()
TE_flag = "No"


# 共有ブロックスイッチ & エイリアス（初期化）
__RUN_HEAVY = False
__MODE_SV = False
__MODE_NATIVEALIGN = False
__DF_ALIAS = None


#0 クエリ領域内のhitを集めてくる
if p_args_na is None: # クエリtranscriptとゲノムDBが同生物種の場合、itself_{}.tsvを使わずにクエリ領域hitを特定する
    df_ances_likely = dfsv[(dfsv["subject acc.ver"] == chr_No.loc[0])&
                           (dfsv["start"] >= dfsv["q_whole_start"])&
                           (dfsv["end"] <= dfsv["q_whole_end"])&
                           (dfsv["inversion"] == q_inv_sign)&
                           (dfsv["% identity"] > 90)]
    #                       (dfsv["% identity"] > 95)]
    #cond_id = df_ances_likely["% identity"] >= 90
    #if cond_id.any():
    #    df_ances_likely = df_ances_likely[cond_id]
    df_ances = df_ances_likely.reset_index()
    #print("l-298",df_ances_likely)
    dfsv.loc[dfsv[(dfsv["subject acc.ver"] == chr_No.loc[0])&
                  (dfsv["start"] >= dfsv["q_whole_start"])&
                  (dfsv["end"] <= dfsv["q_whole_end"])&
                  (dfsv["inversion"] == q_inv_sign)&
                  (dfsv["% identity"] > 90)].index,"ances"] = 2 #メモっとく
                  #(dfsv["% identity"] > 95)].index,"ances"] = 2 #メモっとく

##### dfsv.to_csv("test1_{}.tsv".format(dups_txt_file), sep="\t", index=False, mode="w") # p_args_br

#1 集めたhitについて、存在領域をそれぞれ集める
    ances_likely_area = []
    for index,item in df_ances.iterrows():
        ances_likely_area.append(list(range(min(df_ances.loc[index,"start"], df_ances.loc[index,"end"]),
                        max(df_ances.loc[index,"start"], df_ances.loc[index,"end"])+1)))


#2 領域どうしの被りを見る
#2-0 他のどれとも被らないもの
    yosen = []
    repeat_type = []
    for i in range(len(ances_likely_area)):
        yosen_i = []
        r_t_i = []
        for j in range(len(ances_likely_area)):
            if i == j: #i,jが同じ時はスキップ
                yosen_i.append(100)
            else:
                commons_ij = set(ances_likely_area[i]) & set(ances_likely_area[j])
                if len(commons_ij) == 0:
                    yosen_i.append(0)
#2-1-1 ほぼ完全に一致(repeat)。端点についてはプラマイ５塩基まで許容。
                elif abs(min(ances_likely_area[i]) - min(ances_likely_area[j])) <= 5 and abs(max(ances_likely_area[i]) - max(ances_likely_area[j])) <= 5:
                    yosen_i.append(11)
                    if r_t_i == []:
                        r_t_i.append(i)
                    r_t_i.append(j)
#2-1-2 被るが、両端が共に外側
                elif min(ances_likely_area[i]) <= min(ances_likely_area[j]) and max(ances_likely_area[i]) > max(ances_likely_area[j]):
                    yosen_i.append(1)
                elif min(ances_likely_area[i]) < min(ances_likely_area[j]) and max(ances_likely_area[i]) >= max(ances_likely_area[j]):
                    yosen_i.append(1)
##2-1-2 完全に一致(repeat)
#            elif min(ances_likely_area[i]) == min(ances_likely_area[j]) and max(ances_likely_area[i]) == max(ances_likely_area[j]):
#                yosen_i.append(11)
#                if r_t_i == []:
#                    r_t_i.append(i)
#                r_t_i.append(j)
#2-2 被るが、両端が共に内側
                elif min(ances_likely_area[i]) > min(ances_likely_area[j]) and max(ances_likely_area[i]) < max(ances_likely_area[j]):
                    yosen_i.append(3)
#2-3 被り、かつ、一端は外側だがもう一端は内側
                else:
                    diff_ij = set(ances_likely_area[i]) - set(ances_likely_area[j])
                    diff_ji = set(ances_likely_area[j]) - set(ances_likely_area[i])
                    if diff_ij > diff_ji:
                        yosen_i.append(21)
                    elif diff_ij < diff_ji:
                        yosen_i.append(23)
                    else:
                        yosen_i.append(22)
        yosen.append(yosen_i)
        if r_t_i != []:
            r_t_i.sort()
            repeat_type.append(r_t_i)

#print("repeat_type is ...")
#print(repeat_type)
    set_r_t = list(map(list, set(map(tuple, repeat_type))))
    repeat_num = list(itertools.chain.from_iterable(set_r_t))
#print("set_r_t is ...")
#print(set_r_t)

#3 決勝
    y_winner = []
    ances_num =[]
    for k in range(len(ances_likely_area)):
        if 11 in yosen[k]:
            y_winner.append(k)
    #elif 22 in yosen[k]:
    #    y_winner.append(k)
        else:
            y_winner.append(-1)




    if set_r_t != []:
        for l in range(len(set_r_t)):
            for index,item in df_ances.iterrows():
                if index in set_r_t[l]:
                    df_ances.loc[index,"repeat"] = l

### repeatに関しては、一つ選んで残り全てloseとして扱うことにする。
    if set_r_t != []:
        for m in range(len(set_r_t)):
            df_m = df_ances[df_ances["repeat"] == m+1]
            min_list = list(df_m["mismatches"][df_m["mismatches"] == df_m["mismatches"].min()].index)
            non_min_list = list(df_m["mismatches"][df_m["mismatches"] != df_m["mismatches"].min()].index)
            for item in non_min_list:
                y_winner[item] = "lose"
            #if len(min_list) != 1:
            #    for item in min_list[1:]:
            #        y_winner[item] = "lose" ### ここは妥協。同じならどれか一つ適当に選べば良い。


#print("y_winner is ...")
#print(y_winner)




    for n in range(len(ances_likely_area)):
        if 3 in yosen[n]:
            y_winner[n] = "lose"
        elif 23 in yosen[n]:
            y_winner[n] = "lose"




    for o in range(len(ances_likely_area)):
        if type(y_winner[o]) == int:
            ances_num.append(o)

#print("ances_num is ...")
#print(ances_num)





    for index,item in df_ances.iterrows():
        if index in ances_num:
            df_ances.loc[index,"ances"] = 1

    #print("df_ances is ...")
    #print(df_ances)




    df_ances_2 = df_ances[df_ances["ances"] == 1]
    df_ances_2.rename(columns={"index": "index1"}, inplace=True)
#print("l313: df_ances_2 is ...")
#print(df_ances_2)
    if q_inv == -1:
    #predf_a2 = df_ances_2.sort_values(by = ["start", "q_head"],
    #                                  ascending = [False, True])

    #df_ances_2 = predf_a2.reset_index()
        predf_a2 = df_ances_2.sort_values(by = ["index1"],
                                      ascending = [False])
    else:
        predf_a2 = df_ances_2.sort_values(by = ["index1"],
                                      ascending = [True])

    df_ances_2 = predf_a2.reset_index()

#print("l328: again:df_ances_2 is ...")
#print(df_ances_2)

    non_ances = []
    for index,item in df_ances_2.iterrows():
        if index != df_ances_2.index.max():
            if df_ances_2.loc[index,"start"] != df_ances_2.loc[index+1,"start"]: #自身と1行下が同じリピートではない時
                if index != df_ances_2.index.min():
                    confirmed_num = df_ances_2.loc[:(index-1)][df_ances_2.loc[:(index-1),"ances"] == 1].index.max()
                    if abs(df_ances_2.loc[index,"q_head"] - df_ances_2.loc[confirmed_num,"q_tail"]) > 30: #15:278550x 10?
                        non_ances.append(index)
                        df_ances_2.loc[index,"ances"] = 11 #元は1が入っているので、11とかに書き換えて、スキップされた印にする。
                            #if index + TE_counter == df_ances_2.index.max(): #TE_counterぶん下がって最終行になる時はbreak
                             #   break
                            #else:
                             #   TE_counter += 1
               # else:
                    #print("TE index is ...")
                    #print(index)
                    #print(df_ances_2.loc[index])
                #    df_ances_2.loc[index,"ances"] = 11 #元は1が入っているので、11とかに書き換えて、スキップされた印にする。
                    #print("again:df_ances_2.loc[index] is ...")
                    #print(df_ances_2.loc[index])


            else: #自身と1行下が同じリピートの時
                repeat_counter = 1
                while repeat_counter != 0:
                    if df_ances_2.loc[index,"start"] == df_ances_2.loc[index+repeat_counter+1,"start"]:
                        repeat_counter += 1
                        if index + repeat_counter + 1 == df_ances_2.index.max():
                            break
                    else:
                        break
            #if index + repeat_counter - 1 != df_ances_2.index.max():
                if index != df_ances_2.index.min(): #つまり0行めじゃなければ.
                    confirmed_num = df_ances_2.loc[:(index-1)][df_ances_2.loc[:(index-1),"ances"] == 1].index.max()
                    for p in range(repeat_counter + 1):
                        if abs(df_ances_2.loc[index+p,"q_head"] - df_ances_2.loc[confirmed_num,"q_tail"]) > 10:
                            df_ances_2.loc[index+p,"ances"] = 12 #元は1、スキップ印として12を入れる
                else:
                    for q in range(repeat_counter - 1):
                        df_ances_2.loc[index+1+q,"ances"] = 12 #元は1、スキップ印として12を入れる




#print("non_ances is ...")
#print(non_ances)





    df_ances_rank = df_ances_2[df_ances_2["ances"] == 1]
    df_ances_rank.rename(columns={"index":"index2"}, inplace=True)
    df_ances_rank = df_ances_rank.reset_index()
    rank_num = pd.RangeIndex(start=1, stop=(df_ances_rank["ances"].sum())+1, step=1)

#print("rank_num is ...")
#print(rank_num)
    df_ances_rank["rank"] = rank_num

#print("df_ances_rank is ...")
#print(df_ances_rank)



    df_ances_rank.set_index("index", inplace=True)
    df_ances_rank.rename(columns={"index":"index_x"}, inplace=True)
    df_ances_rank.rename(columns={"index2":"index"}, inplace=True)
#print("l418: df_ances_rank is ...")
#print(df_ances_rank)
    df_ances_2.loc[df_ances_2[df_ances_2["ances"] == 1].index] = df_ances_rank

#print("l422: df_ances_2 is ...")
#print(df_ances_2)



    if q_inv == -1:
    #predf_a = df_ances_2.sort_values(by = ["start", "q_head"],
    #                                  ascending = [True, False])
        predf_a = df_ances_2.sort_values(by = ["index1"],
                                      ascending = [True])
    #print("predf_a is ...")
    #print(predf_a)
        df_ances_2 = predf_a.reset_index(drop = True)


    df_ances_2.set_index("index", inplace=True)
#print("again2:df_ances_2 is ...")
#print(df_ances_2)
#df_ances_2 = df_ances_2.drop(["index"], axis=1, inplace=True)
    df_ances.rename(columns={"index":"index1"}, inplace=True)

#print("again2:df_ances is ...")
#print(df_ances)
    df_ances.loc[df_ances[df_ances["ances"] == 1].index] = df_ances_2


#print("again2:df_ances is ...")
#print(df_ances)
    df_ances.set_index("index1", inplace=True)



    dfsv.loc[dfsv[dfsv["ances"] == 2].index] = df_ances
    dfsv.loc[dfsv[dfsv["ances"] == 1].index,"name"] = dfsv.loc[dfsv[dfsv["ances"] == 1].index,"name"] + "native_align"

#
#dfsv_itself = dfsv[dfsv["ances"] == 1]
#dfsv_itself.reset_index(drop=True)
#for index,item in dfsv_itself.iterrows():
#    if index != (max(rank_num)-1):
#        dfsv_itself.loc[index+1,"start"] - dfsv_itself.loc[index,"end"]
    dfsv = dfsv.dropna(how='any').reset_index(drop=True)
    int_cols = ["start", "end", "ances", "repeat", "rank", "from", "q_whole_start", "q_whole_end",
                "alignment length", "mismatches", "gap opens", "q_head", "q_tail", "s_head", "s_tail"]
    dfsv[int_cols] = dfsv[int_cols].astype('Int64')


    elapsed_l518 = time.time() - start_time
    print(f"DataFrame 読み込み完了: {elapsed_l518:.2f} 秒")

    elapsed_l518 = time.time() - start_time
    t_l518 = time.time()
    t_l178_518 = t_l518 - t_l178
    dt_l518 = datetime.datetime.now()
    print(f"l178-518, done: {t_l178_518:.2f} 秒, 4GTF_startから: {elapsed_l518:.2f} 秒, 現在: {dt_l518}")


    repeat_df = dfsv[(dfsv["repeat"]>0) & (dfsv["rank"]==0)]
### 同じ塩基配列のエキソンを繰り返し抱える、NBPFタイプの遺伝子については解析対象外とし、以降の操作を行わず、仮のGTFをそのまま出力
    if not repeat_df.empty:
        repeat_df.to_csv("Repeat_{}_{}.tsv".format(p_args_br, s6split.iloc[0,1]), sep="\t", index=False, mode="w")
        query_type = "multi-exon"
        pd.DataFrame([{"query_tra_ID": sss.iloc[0,0],
                       "query_gene_name": s6split.iloc[0,1],
                       "query_type" : query_type,
                       f"{p_args_sl}_total_CN": -1,
                       f"{p_args_sl}_SD_CN"   : -1,
                       f"{p_args_sl}_retro_CN": -1,
                       f"({p_args_sl}_single)": -1,
                       f"({p_args_sl}_TE)"    : -1,
                       "repeat?": "Yes"}]).to_csv("statistics_{}_{}.tsv".format(p_args_br, s6split.iloc[0,1]), sep="\t", index=False, mode="w")

    else:
        # ここでは heavy 実行の“準備だけ”を行う（本体は後段で1回だけ実行）
        __RUN_HEAVY = True
        __MODE_SV = True
        __MODE_NATIVEALIGN = False
        __DF_ALIAS = dfsv

else: #クエリtranscriptとゲノムDBが異生物種の場合、itself_{}.tsvを読み込んで使う

    # itself_{}.tsvの読み込み
    #with open(p_args_na, 'r') as nal:
    #    native_align_list = nal.read().split("\n")
    #print(native_align_list) # [itself_ENSTxxxxx......, itself_ENSTxxxxx......, ......, ......]
    #fnsplit = p_args_br.split("_") # fnsplit[0] == ENST00000xxxxx.......
    #native_align_name = [fn for fn in native_align_list if fnsplit[0] in fn]
    #print("native_align_name is ...")
    #print(native_align_name)
    #df_native_align = pd.read_table(native_align_name[0])#,
    df_native_align = pd.read_table(p_args_na)#,
    if df_native_align.loc[0, "rank"] != 1:
        df_native_align = df_native_align.sort_values(by = ["rank"])
                                    # ascending = [True, True])
        df_native_align = df_native_align.reset_index(drop = True)
    #print("df_native_align", df_native_align)
    #fnsplit = p_args_na.split("_")
    #print("fnsplit", fnsplit) #
    #enst_tokens = [x for x in fnsplit if x.startswith("ENST")]
    #ENST_ID = {m.group(0) for x in enst_tokens for m in [re.match(r'^ENST\d+\.\d+', x)] if m}
    #print("ENST_ID", ENST_ID)
    #native_align_name = [fn for fn in native_align_list if any(a in fn for a in ENST_ID)]
    #native_align_name = [fn for fn in native_align_list if fnsplit[0] in fn]
    #print("native_align_name is ...")
    #print(native_align_name)
    #df_native_align = pd.read_table(native_align_name[0])#,
    #if df_native_align.loc[0, "rank"] != 1:
    #    df_native_align = df_native_align.sort_values(by = ["rank"])
    #                                # ascending = [True, True])
    #    df_native_align = df_native_align.reset_index(drop = True)

    __RUN_HEAVY = True
    __MODE_SV = False
    __MODE_NATIVEALIGN = True
    __DF_ALIAS = df_native_align


### NBPFタイプ以外の遺伝子

if __RUN_HEAVY:
    # エイリアスとモードフラグ
    MODE_SV = __MODE_SV
    MODE_NATIVEALIGN = __MODE_NATIVEALIGN
    print("__DF_ALIAS", __DF_ALIAS)
    ### ancesエキソンの個数が１、すなわちクエリがシングルエキソンの場合
    if __DF_ALIAS["ances"].value_counts()[1] == 1:
        print("single")
        query_type = "single-exon"
        if MODE_SV:
            dfsv[dfsv["ances"] == 1] = dfsv[dfsv["ances"] == 1].replace({"from": {0: 1}})
        dfsv_nonan = dfsv[dfsv["rank"] != 1]
        #print("l-644", dfsv[dfsv["ances"] == 1])
        #dfsv.to_csv("testif_{}.tsv".format(p_args_br), sep="\t", index=False, mode="w")

        print("l-672", dfsv["q_head"], dfsv["q_head"].dtype)
        print("NA count:", dfsv["q_head"].isna().sum())

        ### TE部位を特定して除去する
        # NumPy 配列に変換し形状を (行数,1) に
        v_q_head = dfsv["q_head"].to_numpy().reshape(-1, 1)
        v_q_tail = dfsv["q_tail"].to_numpy().reshape(-1, 1)
        # 集計したい整数の範囲
        if MODE_SV:
            full_l = np.arange(1, dfsv.loc[dfsv["ances"] == 1, "q_tail"].item() + 1)
        else:
            full_l = np.arange(1, df_native_align.loc[df_native_align["ances"] == 1, "q_tail"].item() + 1)
        print("full_l : ", full_l)

        # ブロードキャストで各行 × 各整数 のマスク配列 (行数×itself塩基数)
        mask = (full_l >= v_q_head) & (full_l <= v_q_tail)
        print("mask : ", mask)

        # 行方向に True の数を数えて、各整数の頻度を得る
        freq = mask.sum(axis=0)

        # pandas Series／DataFrame にまとめる
        freq_series = pd.Series(freq, index=full_l, name='Frequency')
        freq_df = freq_series.to_frame()  # 必要に応じて DataFrame 化

        # 結果表示（例: 最初と最後の数行だけ）
        print(freq_df.head())
        print(freq_df.tail())
        print(freq_df)

        #plt.figure()
        #plt.plot(freq_series.index, freq_series.values)
        #plt.xlabel('Value')
        #plt.ylabel('Frequency')
        #plt.title('Frequency of integers')
        #plt.tight_layout()
        #plt.show()

        # 1) 閾値の設定とマスク作成
        dup_threshold = p_args_Tt # 101回以上のコピーはTEとして扱う
        t_mask = freq_series > dup_threshold

        # 2) 値とマスク情報を DataFrame にまとめ、連続セグメントにラベル付け
        df_temp = pd.DataFrame({
            'value': freq_series.index,
            'mask': t_mask
            })
        df_temp['segment'] = df_temp['mask'].ne(df_temp['mask'].shift()).cumsum()

        # 3) mask=True（山の部分）に絞り、segment 単位で value をリスト化
        groups = (
            df_temp[df_temp['mask']]
            .groupby('segment')['value']
            .apply(list)
            .tolist()
            )

        print("山ごとの整数リスト:", groups)
        for n in range(len(groups)):
            print("閾値回数以上出現するエキソン : ", n+1, ". - ", groups[n])

        if groups != []:
            for index,item in dfsv_nonan.iterrows():
                hit_area = list(range(dfsv_nonan.loc[index, "q_head"], dfsv_nonan.loc[index,"q_tail"]+1))
                for g in range(len(groups)):
                    common_area_nn = set(hit_area) & set(groups[g])
                    if len(common_area_nn) / len(hit_area) > 0.9: # hit配列の9割以上がTE成分
                        dfsv_nonan.loc[index, "from"] += 100000000

            TE_dfsv = dfsv_nonan[(dfsv_nonan["ances"] != 1) & (dfsv_nonan["from"] >= 100000000)]
            #print(TE_dfsv)
            #mistakes? TE_dfsv = TE_dfsv[(TE_dfsv["ances"] != 1) & (TE_dfsv["rank"] >= 100000000)]
            TE_dfsv.to_csv("TE_{}_{}.tsv".format(p_args_br, s6split.iloc[0,1]), sep="\t", index=False, mode="w")

            #mistakes? dfsv = dfsv[dfsv["rank"] < 100000000]
            #mistakes? dfsv = dfsv[(dfsv["from"] == 10000) | (dfsv["% identity"] >= 95) | (dfsv["ances"] == 1)]
        else:
            TE_dfsv = pd.DataFrame()
            with open("TE_{}_{}.tsv".format(p_args_br, s6split.iloc[0,1]), mode="w") as file:
                file.write("no TEs in this result" + "\n")

        dfsv_nonan2 = dfsv_nonan[(dfsv_nonan["ances"] == 1) | (dfsv_nonan["from"] < 100000000)]
        #dfsv_old1 = dfsv_nonan
        #dfsv_old1.to_csv("test751_{}.tsv".format(p_args_br), sep="\t", index=False, mode="w")
        #print(dfsv_nonan2)
        #dfsv_nonan.reset_index(drop=True)
        dfsv_nonan2 = dfsv_nonan2.reset_index()
        #print("l-759", dfsv_nonan2)
        dfsv_nonan2.rename(columns={"index":"idx2"}, inplace=True)
        #dfsv.to_csv("test758_{}.tsv".format(p_args_br), sep="\t", index=False, mode="w")

        elapsed_l655 = time.time() - start_time
        #t_l120 = time.time()
        #t_l77_120 = t_l120 - t_l77
        dt_l655 = datetime.datetime.now()
        print(f"l655, 4GTF_startから: {elapsed_l655:.2f} 秒, 現在: {dt_l655}")


        # 染色体ごとに分けて、距離が近いものを遺伝子領域としてまとめていく
        chr_ndarray = dfsv_nonan2["subject acc.ver"].unique()
        chr_ndarray_num = dfsv_nonan2["subject acc.ver"].nunique()

        from_se_list_box = []

        for v in range(chr_ndarray_num):
            dfsv_chr = dfsv_nonan2[(dfsv_nonan2["subject acc.ver"] == chr_ndarray[v]) & (dfsv_nonan2["rank"] == 0)]
            dfsv_chr = dfsv_chr.reset_index()
            dfsv_chr.rename(columns={"index":"idx"}, inplace=True)
            dfsv_chr = dfsv_chr.replace({"rank": {0: 1}, "from": {0: 5000}})
            #print("dfsv_chr['rank'].value_counts() is ...")
            #print(dfsv_chr["rank"].value_counts())
            #dfsv_chr["rank"] = 1
            #dfsv_chr[dfsv_chr["rank"] == 10000]["rank"] = 1
            dfsv_chr["name"] += "chr"
            dfsv_chr["name"] += dfsv_chr["subject acc.ver"]

            #dfsv_chr.to_csv("test_chr_{0}_{1}.tsv".format(p_args_br,v), sep="\t", index=False, mode="w")

            paralog_num = 1
            single_num = 1
            rank_in_paralog = 2


            gene_gap = 200000#150000

            from_se_list = []
            #print("l-715, dfsv_chr is ......", dfsv_chr)
            for index,item in dfsv_chr.iterrows():
                if dfsv_chr.loc[index, "rank"] == 1: #1個上とまとまらない場合 → 次とまとまればmulti-ex-paralog、まとまらなければretro or single
                    if index == len(dfsv_chr) - 1: #最終行の場合
                        dfsv_chr.loc[index, "name"] = dfsv_chr.loc[index, "name"] + "_retro_" + str(single_num)#retroかsingleか決まらなくない？>retroで仮置き、ラストでSD型の有無で判定する
                        dfsv_chr.loc[index, "from"] = 10000 # 元は-1

                    else:#最終行ではない場合
                        hit_q_area_here = list(range(dfsv_chr.loc[index,"q_head"], dfsv_chr.loc[index,"q_tail"]+1))
                        hit_q_area_next = list(range(dfsv_chr.loc[index+1,"q_head"], dfsv_chr.loc[index+1,"q_tail"]+1))
                        if dfsv_chr.loc[index, "inversion"] == "+":
                            if dfsv_chr.loc[index+1, "inversion"] == dfsv_chr.loc[index, "inversion"] \
                                    and dfsv_chr.loc[index+1, "start"] - dfsv_chr.loc[index, "end"] < gene_gap \
                                    and len(set(hit_q_area_here) & set(hit_q_area_next)) / len(hit_q_area_here) < 0.5 \
                                    and len(set(hit_q_area_here) & set(hit_q_area_next)) / len(hit_q_area_next) < 0.5 \
                                    and dfsv_chr.loc[index+1, "q_head"] - dfsv_chr.loc[index, "q_tail"] > -30:
                                dfsv_chr.loc[index, "name"] = dfsv_chr.loc[index, "name"] + "_paralog_" + str(paralog_num)
                                dfsv_chr.loc[index+1, "name"] = dfsv_chr.loc[index+1, "name"] + "_paralog_" + str(paralog_num)
                                dfsv_chr.loc[index+1, "rank"] = rank_in_paralog
                                rank_in_paralog += 1
                            else: # retroではなくsingleの可能性だってある
                                dfsv_chr.loc[index, "name"] = dfsv_chr.loc[index, "name"] + "_retro_" + str(single_num)
                                dfsv_chr.loc[index, "from"] = 10000#
                                single_num += 1

                        else: # inversion が - の場合
                            if dfsv_chr.loc[index+1, "inversion"] == dfsv_chr.loc[index, "inversion"] \
                                    and dfsv_chr.loc[index+1, "start"] - dfsv_chr.loc[index, "end"] < gene_gap \
                                    and len(set(hit_q_area_here) & set(hit_q_area_next)) / len(hit_q_area_here) < 0.5 \
                                    and len(set(hit_q_area_here) & set(hit_q_area_next)) / len(hit_q_area_next) < 0.5 \
                                    and dfsv_chr.loc[index, "q_head"] - dfsv_chr.loc[index+1, "q_tail"] > -30:
                                dfsv_chr.loc[index, "name"] = dfsv_chr.loc[index, "name"] + "_paralog_" + str(paralog_num)
                                dfsv_chr.loc[index+1, "name"] = dfsv_chr.loc[index+1, "name"] + "_paralog_" + str(paralog_num)
                                dfsv_chr.loc[index+1, "rank"] = rank_in_paralog
                                rank_in_paralog += 1
                            else:
                                dfsv_chr.loc[index, "name"] = dfsv_chr.loc[index, "name"] + "_retro_" + str(single_num)
                                dfsv_chr.loc[index, "from"] = 10000#
                                single_num += 1

                else: #一個上の行のエキソンとまとまった場合 → multi-ex-paralogで確定
                    if index != len(dfsv_chr) - 1: #最終行ではない場合
                        if dfsv_chr.loc[index, "inversion"] == "+":###
                            if dfsv_chr.loc[index+1, "inversion"] == dfsv_chr.loc[index, "inversion"] \
                                    and dfsv_chr.loc[index+1, "start"] - dfsv_chr.loc[index, "end"] < gene_gap:
                                dfsv_chr.loc[index+1, "name"] = dfsv_chr.loc[index+1, "name"] + "_paralog_" + str(paralog_num)
                                dfsv_chr.loc[index+1, "rank"] = rank_in_paralog
                                dfsv_chr.loc[index+1, "from"] = 5000
                                rank_in_paralog += 1
                            else:
                                single_num += 1
                                paralog_num += 1
                                rank_in_paralog = 2
                                if dfsv_chr.loc[index, "inversion"] == "-":
                                    rank_num = dfsv_chr.loc[index, "rank"]
                                    for w in range(rank_num):
                                        dfsv_chr.loc[index - w, "rank"] = w + 1
                    #最終行の場合はそのまま終了

                        else:### inversion が - の場合
                            if dfsv_chr.loc[index+1, "inversion"] == dfsv_chr.loc[index, "inversion"] \
                                    and dfsv_chr.loc[index+1, "start"] - dfsv_chr.loc[index, "end"] < gene_gap:
                                dfsv_chr.loc[index+1, "name"] = dfsv_chr.loc[index+1, "name"] + "_paralog_" + str(paralog_num)
                                dfsv_chr.loc[index+1, "rank"] = rank_in_paralog
                                rank_in_paralog += 1
                            else:
                                single_num += 1
                                paralog_num += 1
                                rank_in_paralog = 2
                                if dfsv_chr.loc[index, "inversion"] == "-":
                                    rank_num = dfsv_chr.loc[index, "rank"]
                                    for w in range(rank_num):
                                        dfsv_chr.loc[index - w, "rank"] = w + 1
                    #最終行の場合はそのまま終了


                from_se_list.append([dfsv_chr.loc[index, "q_head"], dfsv_chr.loc[index, "q_tail"]])

            from_se_list_box.append(from_se_list)
            #print(from_se_list)
            #print(from_se_list_box)
            #dfsv_chr = dfsv_chr[dfsv_chr["q_tail"]>=150]
            dfsv_chr.set_index("idx", inplace=True)


            #df_ances.rename(columns={"index":"index1"}, inplace=True)
            #print("l-798, dfsv_chr is ......", dfsv_chr)
            dfsv_nonan2.loc[dfsv_nonan2[(dfsv_nonan2["subject acc.ver"] == chr_ndarray[v]) & (dfsv_nonan2["rank"] == 0)].index] = dfsv_chr
            #print("l-800, dfsv_nonan2 is ......", dfsv_nonan2)
            dfsv_nonan2_tmp = dfsv_nonan2.set_index("idx2", drop=True)
            #print("l-885, dfsv_nonan2_tmp is ......", dfsv_nonan2_tmp)
            dfsv_nonan.loc[dfsv_nonan2_tmp.index] = dfsv_nonan2_tmp
            #print("l-887, dfsv_nonan is ......", dfsv_nonan)
            dfsv.loc[dfsv_nonan.index] = dfsv_nonan
            #print("l-802, dfsv is ......", dfsv, dfsv[dfsv["rank"]==1])
            #dfsv.loc[dfsv[(dfsv["subject acc.ver"] == chr_ndarray[v]) & (dfsv["rank"] == 0)].index] = dfsv_chr
            #dfsv.loc[dfsv[(dfsv["subject acc.ver"] == chr_ndarray[v]) & ((dfsv["rank"] == 0) | (dfsv["rank"] == 10000))].index] = dfsv_chr

            #dfsv = dfsv[dfsv["q_tail"] >= 150]

        #dfsv = dfsv.dropna(how='any').reset_index(drop=True)
        int_cols = ["start", "end", "ances", "repeat", "rank", "from", "q_whole_start", "q_whole_end",
                    "alignment length", "mismatches", "gap opens", "q_head", "q_tail", "s_head", "s_tail"]
        dfsv[int_cols] = dfsv[int_cols].astype('Int64')




    ### ancesの番号が１でない、すなわちクエリがマルチエキソンの場合
    else:
        query_type = "multi-exon"
        print("multi")
        #print(dfsv[dfsv["ances"]==1])

        if __MODE_NATIVEALIGN:
            df_ances_rank = df_native_align
        # クエリのexon領域を確定させてリストに入れる
        ances_area = []
        for index,item in df_ances_rank.iterrows():
            #exon_area = []
            exon_area = list(range(df_ances_rank.loc[index,"q_head"], df_ances_rank.loc[index,"q_tail"]+1))
            ances_area.append(exon_area)


        #print("ances_area is ...")
        #print(ances_area)
        #print("ances_area[0] is ...")
        #print(ances_area[0])


        # hitについて、queryのexonと比較し、共通部分を取得する
        common_box = []
        for index,item in dfsv.iterrows():
            common_area = []
            if dfsv.loc[index,"rank"] > 0: # query自身の場合はスキップ
                dfsv.loc[index,"from"] = dfsv.loc[index,"rank"]
                for r in range(len(ances_area)):
                    common_area_r = []
                    common_area.append(common_area_r)

            else: # query自身でないhitの場合の処理
                hit_area = list(range(dfsv.loc[index,"q_head"], dfsv.loc[index,"q_tail"]+1))
                #for r in len(df_ances_rank["ances"].sum()):
                for r in range(len(ances_area)):
                    common_area_r = set(hit_area) & set(ances_area[r])
                    common_area.append(list(common_area_r))
            common_box.append(common_area)

        #print("common_box is ...")
        #print(common_box)


        # hit一つと、queryに含まれる各exonとの共通部分を比較して、どのexonと一致するか（exon番号がどれになるのか）を決定するにあたって、
        # query内の複数のexonの成分を持つhit（＝共通部分の長さが11以上である配列、が複数個存在したhit）については、レトロポゾン型ということにする
        common_judge_box = []
        nzcj_box = [] # 確か、non_zero_common_judgeの頭文字
        for s in range(len(common_box)):
            common_judge = []
            for t in range(len(common_box[s])):
                if len(common_box[s][t]) <= 10: #元は10
                    common_judge.append(0)
                else:
                    common_judge.append(len(common_box[s][t]))
            common_judge_box.append(common_judge)
            com_num = np.count_nonzero(common_judge)
            nzcj = np.nonzero(common_judge)
            nzcj_box.append(nzcj)
            #print("common_judge is ...")
            #print(common_judge)
            if com_num > 2:
                #print("nzcj is ...")
                #print(nzcj)
                retro_mark = "retro_"
                for u in range(len(nzcj[0])): # retro_+1+2+3+4、とかがわかるようにする
                    retro_mark += ("+" + (nzcj[0][u]+1).astype("str"))
                #dfsv.loc[s,"retro"] = retro_mark
                dfsv.loc[s,"name"] += retro_mark
                dfsv.loc[s,"from"] = 10000
                dfsv.loc[s,"rank"] = 1 # 元は1でやってた
            elif com_num == 2:
                #true_from = common_judge.index(max(common_judge))
                true_from = common_judge.index(sorted(common_judge)[-1])
                false_from = common_judge.index(sorted(common_judge)[-2])
                # 「ヒットとの一致塩基数が少ない方のクエリ」との一致領域が短すぎる（一致塩基数が35bp以下）場合 -> レトロ型ではなく1エキソン # 20bpから変更
                if sorted(common_judge)[-2] <= 35:
                    dfsv.loc[s,"from"] = true_from + 1
                    dfsv.loc[s,"rank"] = 0
                    common_judge[false_from] = 0
                    common_judge_box[s] = common_judge
                    nzcj = np.nonzero(common_judge)
                    nzcj_box[s] = nzcj
                # そうでない場合 -> ちゃんとレトロ型（exon数：2）
                else:
                    #print("nzcj is ...")
                    #print(nzcj)
                    retro_mark = "retro_"
                    for u in range(len(nzcj[0])):
                        retro_mark += ("+" + (nzcj[0][u]+1).astype("str"))
                    #dfsv.loc[s,"retro"] = retro_mark
                    dfsv.loc[s,"name"] += retro_mark
                    dfsv.loc[s,"from"] = 10000
                    dfsv.loc[s,"rank"] = 1 # 元は1でやってた


        #print("dfsv is ...")
        #print(dfsv)
        dfsv_retro = dfsv[dfsv["from"] == 10000]
        #print("dfsv_retro is ...")
        #print(dfsv_retro)
        dfsv_retro = dfsv_retro.reset_index()
        dfsv_retro.rename(columns={"index":"idx"}, inplace=True)
        #print("again:dfsv_retro is ...")
        #print(dfsv_retro)
        for index,item in dfsv_retro.iterrows():
            dfsv_retro.loc[index,"name"] += ("_" + str(index+1))

        dfsv_retro.set_index("idx", inplace=True)
        dfsv.loc[dfsv[dfsv["from"] == 10000].index] = dfsv_retro

        #print("nzcj_boxi & dfsvis ...")
        #print(nzcj_box)
        #print(nzcj_box[0])
        #print(dfsv)

        #print("common_judge_box is ...")
        #print(common_judge_box)
        #print("common_judge is ...")
        #print(common_judge)

        #com_num = np.count_nonzero(common_judge)
        #print("com_num is ...")
        #print(com_num)

        #dfsv.to_csv("testtest_{}.tsv".format(p_args_br), sep="\t", index=False, mode="w")



        # それぞれのhitが、どのexonから飛んできたのか、を判定する
        for index,item in dfsv.iterrows():
            if dfsv.loc[index,"rank"] == 0:
                #print(nzcj_box[index][0])
                dfsv.loc[index,"from"] = nzcj_box[index][0] + 1

        #print("dfsv is ...")
        #print(dfsv)

        #dfsv.to_csv("test705_{}.tsv".format(p_args_br), sep="\t", index=False, mode="w")

        dfsvvc = dfsv["from"].value_counts().sort_index()
        dfsvvc_list = dfsvvc.to_list()
        dfsvvc_dict = dfsvvc.to_dict()
        #print(dfsvvc)
        #print(dfsvvc_list)
        #print(dfsvvc_dict)
        dfsvvc.to_csv("exon_CN_stat_{}_{}.tsv".format(p_args_br, s6split.iloc[0,1]), sep="\t", mode="w")

        sorted_dict = dict(sorted(dfsvvc_dict.items(), key=lambda item: item[0]))
        sorted_dict_items = list(sorted_dict.items())
        #print(sorted_dict)
        #print(sorted_dict_items)
        #print(len(sorted_dict_items)-1)
        #print(range(len(sorted_dict_items)-1))

        TE_dfsv = pd.DataFrame()
        TE_flag = "No"
        # 重複数が最小のエキソンでさえ100回より多く重複していたら、後回しにする
        if sorted(dfsvvc_list)[0] > 100:
        #if sorted_dict_items[0][1] > 100:
            print("over 100 dups")
            list_100overexons = []
            TE_dfsv = dfsv[dfsv["ances"] != 1]
            TE_dfsv.to_csv("TE_{}_{}.tsv".format(p_args_br, s6split.iloc[0,1]), sep="\t", index=False, mode="w")
            TE_flag = "Yes"
        else:
            list_100overexons = []
            for i in range(len(sorted_dict_items)-1):
                xdup_diff = abs(sorted_dict_items[i+1][1] - sorted_dict_items[i][1])
                print(i)
                print(list_100overexons)
                if xdup_diff > 100:
                    # xdup_diffが大きくなる方のキーを取得
                    larger_key = sorted_dict_items[i][0] if sorted_dict_items[i][1] > sorted_dict_items[i + 1][1] else sorted_dict_items[i + 1][0]
                    print("差分が100を越えたエキソン", i, larger_key)
                    list_100overexons.append(larger_key)
                    #print("差分が100を越えたエキソン", i, sorted_dict_items[max(sorted_dict_items[i][1], sorted_dict_items[i+1][1])])
                    #list_100overexons.append(max(sorted_dict_items[i][1], sorted_dict_items[i+1][1]))
                elif i != 0 and list_100overexons != [] and list_100overexons[-1] == sorted_dict_items[i][0] and sorted_dict_items[i+1][1] > 100:
                    print("差分100以下だがTE疑惑のエキソン", i, sorted_dict_items[i+1][0])
                    list_100overexons.append(sorted_dict_items[i+1][0])
            if MODE_NATIVEALIGN:
                list_100overexons.extend((df_native_align.loc[df_native_align["from"] > 100000000, "from"] - 100000000).tolist())
                print("l-1100,", list_100overexons)
            # 重複回数が100回を超えているexonについては、TEだとして一旦除去する
            list_100overexons = list(set(list_100overexons))
            print(list_100overexons)
            if list_100overexons != []:
                for item in list_100overexons:
                    if item < 100000000:
                        item_100000000 = item + 100000000
                    dfsv_TE = dfsv[dfsv["from"] == item]
                    dfsv_TE["from"] = item_100000000
                    dfsv.loc[dfsv[dfsv["from"] == item].index] = dfsv_TE
                list_100overexons = [item + 100000000 for item in list_100overexons if item < 100000000]
                list_100overexons = list(dict.fromkeys(list_100overexons))
                list_100overexons.sort()
                print("l-773", list_100overexons)

                #dfsv_TE.to_csv("test_dfsv_TE_738_{}.tsv".format(p_args_br), sep="\t", index=False, mode="w")
                #dfsv.to_csv("test739_{}.tsv".format(p_args_br), sep="\t", index=False, mode="w")

            #mistakes? TE_dfsv = dfsv[dfsv["from"] != 10000] # レトロ型コピーを除く
            #mistakes? TE_dfsv = TE_dfsv[(TE_dfsv["ances"] != 1) & (TE_dfsv["from"] >= 100000000)]
                TE_dfsv = dfsv[(dfsv["ances"] != 1) & (dfsv["from"] >= 100000000)]
                #print(TE_dfsv)
            #mistakes? TE_dfsv = TE_dfsv[(TE_dfsv["ances"] != 1) & (TE_dfsv["rank"] >= 100000000)]
                TE_dfsv.to_csv("TE_{}_{}.tsv".format(p_args_br, s6split.iloc[0,1]), sep="\t", index=False, mode="w")

            #mistakes? dfsv = dfsv[dfsv["rank"] < 100000000]
            #mistakes? dfsv = dfsv[(dfsv["from"] == 10000) | (dfsv["% identity"] >= 95) | (dfsv["ances"] == 1)]
            else:
                with open("TE_{}_{}.tsv".format(p_args_br, s6split.iloc[0,1]), mode="w") as file:
                    file.write("no TEs in this result" + "\n")
            dfsv = dfsv[(dfsv["ances"] == 1) | (dfsv["from"] < 100000000)]
            dfsv_old1 = dfsv
            #dfsv_old1.to_csv("test751_{}.tsv".format(p_args_br), sep="\t", index=False, mode="w")
            #print(dfsv)
            dfsv.reset_index(drop=True)
            #dfsv.to_csv("test758_{}.tsv".format(p_args_br), sep="\t", index=False, mode="w")




        # 染色体ごとに分けて、距離が近いものを遺伝子領域としてまとめていく
        chr_ndarray = dfsv["subject acc.ver"].unique()
        chr_ndarray_num = dfsv["subject acc.ver"].nunique()

        #print("chr_ndarray is ...")
        #print(chr_ndarray)
        #print(type(chr_ndarray))


        exon_gap_list = []
        exon_len = []
        if MODE_SV:
            df_ances_4cul_distance = dfsv[dfsv["ances"] == 1].sort_values(by = ["rank"], ascending = [True]).reset_index(drop = True) # invにかかわらずrank昇順にする
        else:
            df_ances_4cul_distance = df_native_align.sort_values(by = ["rank"], ascending = [True]).reset_index(drop = True) # invにかかわらずrank昇順にする
        #print(df_ances_4cul_distance)
        for index,item in df_ances_4cul_distance.iterrows():
            if index != len(df_ances_4cul_distance) - 1:
                exon_gap_list.append(abs(df_ances_4cul_distance.loc[index+1, "start"] - df_ances_4cul_distance.loc[index, "end"]))
            else:
                exon_len = df_ances_4cul_distance["q_tail"].max()
        print(exon_gap_list)
        print(exon_len)

        for v in range(chr_ndarray_num):
            Te_check_box = []
            #dfsv_chr = dfsv[(dfsv["subject acc.ver"] == chr_ndarray[v]) & ((dfsv["rank"] == 0) | (dfsv["rank"] == 10000))]
            dfsv_chr = dfsv[(dfsv["subject acc.ver"] == chr_ndarray[v]) & (dfsv["rank"] == 0)]
            #print(chr_ndarray[v], dfsv_chr)
            if dfsv_chr.empty == False:
                # TEの中から正しいエキソンを救出して戻す
                #tmp_dfsv_chr[(dfsv_chr["from"] == list_100overexons[0]-1) & (dfsv_chr["from"] == list_100overexons[0]+1)]

                dfsv_chr = dfsv_chr.reset_index()
                dfsv_chr.rename(columns={"index":"idx"}, inplace=True)
                dfsv_chr = dfsv_chr.replace({"rank": {0: 1}})
                #print("dfsv_chr is ......")
                #print(dfsv_chr)
                #print("dfsv_chr['rank'].value_counts() is ...")
                #print(dfsv_chr["rank"].value_counts())
                #dfsv_chr["rank"] = 1
                #dfsv_chr[dfsv_chr["rank"] == 10000]["rank"] = 1
                dfsv_chr["name"] += "chr"
                dfsv_chr["name"] += dfsv_chr["subject acc.ver"]

                #dfsv_chr.to_csv("test_chr_{0}_{1}.tsv".format(p_args_br,v), sep="\t", index=False, mode="w")

                from_list = []
                paralog_num = 1
                single_num = 1
                rank_in_paralog = 2
                locus_len = sum(exon_gap_list, exon_len)*1.2 # 1.2倍まで許容ということにする


                for index,item in dfsv_chr.iterrows():
                    if from_list == []: #一個上の行のエキソンとまとまらなかった場合
                        from_list.append(dfsv_chr.loc[index, "from"])
                        if index == len(dfsv_chr) - 1: #最終行の場合
                            dfsv_chr.loc[index, "name"] = dfsv_chr.loc[index, "name"] + "_single_" + str(single_num)
                            dfsv_chr.loc[index, "from"] = dfsv_chr.loc[index, "from"] * (-1)
                        else:
                            if dfsv_chr.loc[index, "inversion"] == "+":
                                exon_num_diff = dfsv_chr.loc[index+1, "from"] - dfsv_chr.loc[index, "from"]
                                intron_len = 0
                                print("exon_num_diff, exon_gap_list", exon_num_diff, exon_gap_list)
                                for e_n_d in range(exon_num_diff):
                                    print(index, dfsv_chr.loc[index, "from"])
                                    #print(exon_gap_list[dfsv_chr.loc[index, "from"]])
                                    intron_len += exon_gap_list[dfsv_chr.loc[index, "from"] - 1 + e_n_d]
                                    print("intron_len", intron_len)
                                if dfsv_chr.loc[index+1, "inversion"] == dfsv_chr.loc[index, "inversion"] \
                                        and dfsv_chr.loc[index+1, "start"] - dfsv_chr.loc[index, "end"] < locus_len + intron_len \
                                        and dfsv_chr.loc[index+1, "from"] not in from_list \
                                        and dfsv_chr.loc[index+1, "from"] > dfsv_chr.loc[index, "from"]: #一個下の行のエキソンとまとまる場合 #and (dfsv_chr.loc[index+1, "from"] not in from_list or dfsv_chr.loc[index+1, "from"] >= dfsv_chr.loc[index, "from"]): #一個下の行のエキソンとまとまる場合
                                    dfsv_chr.loc[index, "name"] = dfsv_chr.loc[index, "name"] + "_paralog_" + str(paralog_num)
                                    dfsv_chr.loc[index+1, "name"] = dfsv_chr.loc[index+1, "name"] + "_paralog_" + str(paralog_num)
                                    dfsv_chr.loc[index+1, "rank"] = rank_in_paralog # 1行下は2番エキソンで、rank_in_paralogも2なので。
                                    rank_in_paralog += 1
                                    from_list.append(dfsv_chr.loc[index+1, "from"])
                                else:
                                    dfsv_chr.loc[index, "name"] = dfsv_chr.loc[index, "name"] + "_single_" + str(single_num)
                                    dfsv_chr.loc[index, "from"] = dfsv_chr.loc[index, "from"] * (-1)
                                    from_list = []
                                    single_num += 1

                            else:
                                exon_num_diff = dfsv_chr.loc[index, "from"] - dfsv_chr.loc[index+1, "from"]
                                intron_len = 0
                                for e_n_d in range(exon_num_diff):
                                    intron_len += exon_gap_list[dfsv_chr.loc[index+1, "from"] - 1 + e_n_d]
                                if dfsv_chr.loc[index+1, "inversion"] == dfsv_chr.loc[index, "inversion"] \
                                        and dfsv_chr.loc[index+1, "start"] - dfsv_chr.loc[index, "end"] < locus_len + intron_len \
                                        and dfsv_chr.loc[index+1, "from"] not in from_list \
                                        and dfsv_chr.loc[index+1, "from"] < dfsv_chr.loc[index, "from"]: #一個下の行のエキソンとまとまる場合 #and (dfsv_chr.loc[index+1, "from"] not in from_list or dfsv_chr.loc[index+1, "from"] <= dfsv_chr.loc[index, "from"]): #一個下の行のエキソンとまとまる場合
                                    dfsv_chr.loc[index, "name"] = dfsv_chr.loc[index, "name"] + "_paralog_" + str(paralog_num)
                                    dfsv_chr.loc[index+1, "name"] = dfsv_chr.loc[index+1, "name"] + "_paralog_" + str(paralog_num)
                                    dfsv_chr.loc[index+1, "rank"] = rank_in_paralog
                                    rank_in_paralog += 1
                                    from_list.append(dfsv_chr.loc[index+1, "from"])
                                else:
                                    dfsv_chr.loc[index, "name"] = dfsv_chr.loc[index, "name"] + "_single_" + str(single_num)
                                    dfsv_chr.loc[index, "from"] = dfsv_chr.loc[index, "from"] * (-1)
                                    from_list = []
                                    single_num += 1

                    else: #一個上の行のエキソンとまとまった場合
                        if index != len(dfsv_chr) - 1: #最終行ではない場合
                            if dfsv_chr.loc[index, "inversion"] == "+":
                                exon_num_diff = dfsv_chr.loc[index+1, "from"] - dfsv_chr.loc[index, "from"]
                                intron_len = 0
                                for e_n_d in range(exon_num_diff):
                                    intron_len += exon_gap_list[dfsv_chr.loc[index, "from"] - 1 + e_n_d]
                                if dfsv_chr.loc[index+1, "inversion"] == dfsv_chr.loc[index, "inversion"] \
                                        and dfsv_chr.loc[index+1, "start"] - dfsv_chr.loc[index, "end"] < locus_len + intron_len \
                                        and dfsv_chr.loc[index+1, "from"] not in from_list \
                                        and dfsv_chr.loc[index+1, "from"] > dfsv_chr.loc[index, "from"]: #一個下の行のエキソンとまとまる場合 #and (dfsv_chr.loc[index+1, "from"] not in from_list or dfsv_chr.loc[index+1, "from"] >= dfsv_chr.loc[index, "from"]): #一個下の行のエキソンとまとまる場合
                                    dfsv_chr.loc[index+1, "name"] = dfsv_chr.loc[index+1, "name"] + "_paralog_" + str(paralog_num)
                                    dfsv_chr.loc[index+1, "rank"] = rank_in_paralog
                                    rank_in_paralog += 1
                                    from_list.append(dfsv_chr.loc[index+1, "from"])
                                else:
                                    from_list = []
                                    single_num += 1
                                    paralog_num += 1
                                    rank_in_paralog = 2
                                    if dfsv_chr.loc[index, "inversion"] == "-":
                                        rank_num = dfsv_chr.loc[index, "rank"]
                                        for w in range(rank_num):
                                            dfsv_chr.loc[index - w, "rank"] = w + 1

                            else:
                                exon_num_diff = dfsv_chr.loc[index, "from"] - dfsv_chr.loc[index+1, "from"]
                                intron_len = 0
                                for e_n_d in range(exon_num_diff):
                                    intron_len += exon_gap_list[dfsv_chr.loc[index+1, "from"] - 1 + e_n_d]
                                if dfsv_chr.loc[index+1, "inversion"] == dfsv_chr.loc[index, "inversion"] \
                                        and dfsv_chr.loc[index+1, "start"] - dfsv_chr.loc[index, "end"] < locus_len + intron_len \
                                        and dfsv_chr.loc[index+1, "from"] not in from_list \
                                        and dfsv_chr.loc[index+1, "from"] < dfsv_chr.loc[index, "from"]: #一個下の行のエキソンとまとまる場合 #and (dfsv_chr.loc[index+1, "from"] not in from_list or dfsv_chr.loc[index+1, "from"] <= dfsv_chr.loc[index, "from"]): #一個下の行のエキソンとまとまる場合
                                    dfsv_chr.loc[index+1, "name"] = dfsv_chr.loc[index+1, "name"] + "_paralog_" + str(paralog_num)
                                    dfsv_chr.loc[index+1, "rank"] = rank_in_paralog
                                    rank_in_paralog += 1
                                    from_list.append(dfsv_chr.loc[index+1, "from"])
                                else:
                                    from_list = []
                                    single_num += 1
                                    paralog_num += 1
                                    rank_in_paralog = 2
                                    if dfsv_chr.loc[index, "inversion"] == "-":
                                        rank_num = dfsv_chr.loc[index, "rank"]
                                        for w in range(rank_num):
                                            dfsv_chr.loc[index - w, "rank"] = w + 1
                        else:
                            if dfsv_chr.loc[index, "inversion"] == "-":
                                rank_num = dfsv_chr.loc[index, "rank"]
                                for w in range(rank_num):
                                    dfsv_chr.loc[index - w, "rank"] = w + 1



                paralog_num = dfsv_chr[dfsv_chr["name"].str.contains("paralog")]["name"].drop_duplicates().count()
                #if paralog_num == 1:
                #    paralog_num = 0
                print("l-924 paralog_num")
                print(paralog_num)

                if list_100overexons != []:
                    name_unique = dfsv_chr["name"].unique()
                    name_unique_num = dfsv_chr["name"].nunique()
                    TE_dfsv_ri = TE_dfsv.reset_index()
                    TE_dfsv = TE_dfsv_ri
                    TE_dfsv.rename(columns={"index":"idx"}, inplace=True)
                    #print(TE_dfsv)
                    #print(TE_dfsv.columns)
                    #TE_dfsv.drop(columns=['Unnamed: 0'], inplace=True)
                    TE_dfsv.set_index('idx', inplace=True)
                    #print("TE_dfsv is ......")
                    #print(TE_dfsv)
                    #dfsv_chr.set_index("idx", inplace=True)
                    # パラログごとに分けて処理していく
                    for n_u_m in range(name_unique_num):
                        from_list_paralog = []
                        #print(dfsv_chr)
                        dfsv_chr_paralog = dfsv_chr[dfsv_chr["name"] == name_unique[n_u_m]]
                        #dfsv_chr_paralog = dfsv_chr_paralog.reset_index()
                        #dfsv_chr_paralog.rename(columns={"index":"idx2"}, inplace=True)
                        #print("l-936 dfsv_chr_paralog is ......")
                        #print(dfsv_chr_paralog)
                        #print(dfsv_chr_paralog.columns)
                        #print(dfsv_chr_paralog.index)


                        dfsv_chr_paralog.set_index("idx", inplace=True)
                        #print("dfsv_chr_paralog is ......")
                        #print(dfsv_chr_paralog)
                        # 以降、dfsv_paralogが2行以上ある（＝マルチエキソン）場合のみ処理が行われるようにする。
                        # dfsv_chr_paralogが1行しかない場合は、singleのまま終わらせる。
                        if len(dfsv_chr_paralog) >= 2:
                            from_list_paralog = sorted(abs(dfsv_chr_paralog["from"]).to_list())
                            print(from_list_paralog)
                            print(list_100overexons)
                            if len(list_100overexons)==1: # 100回以上重複しているTE疑惑エキソンが1個しかない時
                                #print(dfsv_chr)
                                dfsv_chr.set_index("idx", inplace=True)
                                print("l-996 list_100overexons[0]", list_100overexons[0], list_100overexons)
                                if list_100overexons[0] - 100000000 < from_list_paralog[0]: # TE疑惑エキソンが先頭エキソンである場合
                                    closest_TE_plus1  = min(from_list_paralog, key=lambda x: abs(x - (list_100overexons[0] - 100000000 + 1)))

                                    #print(dfsv_chr_paralog.loc[dfsv_chr_paralog.index.min(), "inversion"])
                                    if dfsv_chr_paralog.loc[dfsv_chr_paralog.index.min(), "inversion"] == "+": # 左端が開放なので、locus_lenだけ手前からということにする
                                        exon_candidate_area_start = dfsv_chr_paralog[abs(dfsv_chr_paralog["from"]) == closest_TE_plus1]["start"].iloc[0] - locus_len
                                        exon_candidate_area_end   = dfsv_chr_paralog[abs(dfsv_chr_paralog["from"]) == closest_TE_plus1]["start"].iloc[0]
                                        # ---- ここから追加（def なし, 隣接パラログで探索窓クランプ in-place） ----
                                        cur_name = name_unique[n_u_m]
                                        _parts = cur_name.split("_")
                                        _is_paralog = (len(_parts) >= 2 and _parts[-2] == "paralog" and _parts[-1].isdigit())
                                        if _is_paralog:
                                            _group_prefix = "_".join(_parts[:-1])
                                            _blocks = (dfsv_chr[dfsv_chr["name"].str.startswith(_group_prefix + "_")]
                                                       .groupby("name", as_index=False)
                                                       .agg(block_start=("start", "min"), block_end=("end", "max"))
                                                       .sort_values("block_start")
                                                       .reset_index(drop=True))
                                            _idx_list = _blocks.index[_blocks["name"] == cur_name].tolist()
                                            if len(_idx_list) > 0:
                                                _i = int(_idx_list[0])
                                                if _i == 0:
                                                    _left_limit = int(_blocks.loc[_i, "block_start"] - locus_len)
                                                else:
                                                    _left_limit = int(_blocks.loc[_i - 1, "block_end"])
                                                if _i == len(_blocks) - 1:
                                                    _right_limit = int(_blocks.loc[_i, "block_end"] + locus_len)
                                                else:
                                                    _right_limit = int(_blocks.loc[_i + 1, "block_start"])
                                                _eca_start = max(exon_candidate_area_start, _left_limit)
                                                _eca_end   = min(exon_candidate_area_end, _right_limit)
                                                if _eca_start < _eca_end:
                                                    exon_candidate_area_start = _eca_start
                                                    exon_candidate_area_end   = _eca_end
                                        # ---- ここまで追加 ----

                                        dfsv_chr_paralog_eca = TE_dfsv[(TE_dfsv["subject acc.ver"] == chr_ndarray[v]) & (TE_dfsv["inversion"] == "+") & (TE_dfsv["start"] > exon_candidate_area_start) & (TE_dfsv["end"] < exon_candidate_area_end)]
                                    else: # 右端が開放なので、locus_lenだけ奥まで、ということにする
                                        exon_candidate_area_start = dfsv_chr_paralog[abs(dfsv_chr_paralog["from"]) == closest_TE_plus1]["end"].iloc[0]
                                        exon_candidate_area_end   = dfsv_chr_paralog[abs(dfsv_chr_paralog["from"]) == closest_TE_plus1]["end"].iloc[0] + locus_len
                                        # ---- ここから追加（def なし, 隣接パラログで探索窓クランプ in-place） ----
                                        cur_name = name_unique[n_u_m]
                                        _parts = cur_name.split("_")
                                        _is_paralog = (len(_parts) >= 2 and _parts[-2] == "paralog" and _parts[-1].isdigit())
                                        if _is_paralog:
                                            _group_prefix = "_".join(_parts[:-1])
                                            _blocks = (dfsv_chr[dfsv_chr["name"].str.startswith(_group_prefix + "_")]
                                                       .groupby("name", as_index=False)
                                                       .agg(block_start=("start", "min"), block_end=("end", "max"))
                                                       .sort_values("block_start")
                                                       .reset_index(drop=True))
                                            _idx_list = _blocks.index[_blocks["name"] == cur_name].tolist()
                                            if len(_idx_list) > 0:
                                                _i = int(_idx_list[0])
                                                if _i == 0:
                                                    _left_limit = int(_blocks.loc[_i, "block_start"] - locus_len)
                                                else:
                                                    _left_limit = int(_blocks.loc[_i - 1, "block_end"])
                                                if _i == len(_blocks) - 1:
                                                    _right_limit = int(_blocks.loc[_i, "block_end"] + locus_len)
                                                else:
                                                    _right_limit = int(_blocks.loc[_i + 1, "block_start"])
                                                _eca_start = max(exon_candidate_area_start, _left_limit)
                                                _eca_end   = min(exon_candidate_area_end, _right_limit)
                                                if _eca_start < _eca_end:
                                                    exon_candidate_area_start = _eca_start
                                                    exon_candidate_area_end   = _eca_end
                                        # ---- ここまで追加 ----

                                        dfsv_chr_paralog_eca = TE_dfsv[(TE_dfsv["subject acc.ver"] == chr_ndarray[v]) & (TE_dfsv["inversion"] == "-") & (TE_dfsv["start"] > exon_candidate_area_start) & (TE_dfsv["end"] < exon_candidate_area_end)]
                                    #print("dfsv_chr_paralog_eca is ......")
                                    #print(dfsv_chr_paralog_eca)
                                    #print("dfsv_chr_paralog is ......")
                                    #print(dfsv_chr_paralog)
                                    if dfsv_chr_paralog_eca.empty == False: # 記述していないが、空だったら処理を終える
                                        #max_bit_score = dfsv_chr_paralog["bit score"].max().index()
                                        max_bit_score_idx = dfsv_chr_paralog_eca["bit score"].idxmax()
                                        #print(max_bit_score_idx)
                                        mbsi_check = list(dfsv_chr_paralog_eca["bit score"][dfsv_chr_paralog_eca["bit score"] == dfsv_chr_paralog_eca["bit score"].max()].index)
                                        #if len(mbsi_check) > 1:
                                            #max_length_idx =
                                        #print(dfsv_chr_paralog_eca["bit score"].dtype)
                                        #print(type(max_bit_score))
                                        #True_exon = dfsv_chr_paralog_eca[dfsv_chr_paralog_eca["bit score"] == max_bit_score]
                                        True_exon = dfsv_chr_paralog_eca.loc[max_bit_score_idx]
                                        #print("True_exon is ......")
                                        #print(True_exon)
                                        #print("Te_check_box is ......")
                                        #print(Te_check_box)
                                        if True_exon.name not in Te_check_box:
                                            Te_check_box.append(True_exon.name)
                                            if "_single_" in name_unique[n_u_m]:
                                                dfsv_chr_paralog["from"] = (-1) * dfsv_chr_paralog["from"]
                                            True_exon["name"] = name_unique[n_u_m]
                                            #print(True_exon)

                                            tmp_dcp_concat = pd.concat([dfsv_chr_paralog, True_exon.set_axis(dfsv_chr_paralog.columns).to_frame().T]).sort_index()
                                            dfsv_chr_paralog = tmp_dcp_concat
                                            #print("l-1042", dfsv_chr_paralog)
                                            #print("l-1043", dfsv_chr)
                                            dfsv_chr.update(dfsv_chr_paralog)
                                            dfsv_chr = pd.concat([dfsv_chr, dfsv_chr_paralog[~dfsv_chr_paralog.index.isin(dfsv_chr.index)]], axis=0)

                                            #tmp_dc_concat = pd.concat([dfsv_chr, dfsv_chr_paralog]).sort_index()
                                            #dfsv_chr = tmp_dc_concat
                                            #print("l-1049", dfsv_chr)

                                            #tmp_dcp_concat = pd.concat([dfsv_chr_paralog, True_exon.set_axis(dfsv_chr_paralog.columns).to_frame().T]).sort_index()
                                            #dfsv_chr_paralog = tmp_dcp_concat
                                            #print(dfsv_chr_paralog)

                                            #tmp_dc_concat = pd.concat([dfsv_chr, True_exon.set_axis(dfsv_chr.columns).to_frame().T]).sort_index()
                                            #dfsv_chr = tmp_dc_concat

                                            #tmp_dfsv_concat = pd.concat([dfsv, True_exon.set_axis(dfsv.columns).to_frame().T]).sort_index()
                                            #dfsv = tmp_dfsv_concat
                                            #print(dfsv)






                                elif list_100overexons[0] - 100000000 > from_list_paralog[-1]: # TE疑惑エキソンが最終エキソンである場合
                                    closest_TE_minus1 = min(from_list_paralog, key=lambda x: abs(x - (list_100overexons[0] - 100000000 - 1)))
                                    #print(closest_TE_minus1)

                                    #print(dfsv_chr_paralog.loc[dfsv_chr_paralog.index.min(), "inversion"])
                                    if dfsv_chr_paralog.loc[dfsv_chr_paralog.index.min(), "inversion"] == "+": # 右端が開放なので、locus_lenだけ奥までということにする
                                        exon_candidate_area_start = dfsv_chr_paralog[abs(dfsv_chr_paralog["from"]) == closest_TE_minus1]["end"].iloc[0]
                                        exon_candidate_area_end   = dfsv_chr_paralog[abs(dfsv_chr_paralog["from"]) == closest_TE_minus1]["end"].iloc[0] + locus_len
                                        # ---- ここから追加（def なし, 隣接パラログで探索窓クランプ in-place） ----
                                        cur_name = name_unique[n_u_m]
                                        _parts = cur_name.split("_")
                                        _is_paralog = (len(_parts) >= 2 and _parts[-2] == "paralog" and _parts[-1].isdigit())
                                        if _is_paralog:
                                            _group_prefix = "_".join(_parts[:-1])
                                            _blocks = (dfsv_chr[dfsv_chr["name"].str.startswith(_group_prefix + "_")]
                                                       .groupby("name", as_index=False)
                                                       .agg(block_start=("start", "min"), block_end=("end", "max"))
                                                       .sort_values("block_start")
                                                       .reset_index(drop=True))
                                            _idx_list = _blocks.index[_blocks["name"] == cur_name].tolist()
                                            if len(_idx_list) > 0:
                                                _i = int(_idx_list[0])
                                                if _i == 0:
                                                    _left_limit = int(_blocks.loc[_i, "block_start"] - locus_len)
                                                else:
                                                    _left_limit = int(_blocks.loc[_i - 1, "block_end"])
                                                if _i == len(_blocks) - 1:
                                                    _right_limit = int(_blocks.loc[_i, "block_end"] + locus_len)
                                                else:
                                                    _right_limit = int(_blocks.loc[_i + 1, "block_start"])
                                                _eca_start = max(exon_candidate_area_start, _left_limit)
                                                _eca_end   = min(exon_candidate_area_end, _right_limit)
                                                if _eca_start < _eca_end:
                                                    exon_candidate_area_start = _eca_start
                                                    exon_candidate_area_end   = _eca_end
                                        # ---- ここまで追加 ----

                                        #print("e_c_a_start")
                                        dfsv_chr_paralog_eca = TE_dfsv[(TE_dfsv["subject acc.ver"] == chr_ndarray[v]) & (TE_dfsv["inversion"] == "+") & (TE_dfsv["start"] > exon_candidate_area_start) & (TE_dfsv["end"] < exon_candidate_area_end)]
                                    else: # 左端が開放なので、locus_lenだけ手前から、ということにする
                                        exon_candidate_area_start = dfsv_chr_paralog[abs(dfsv_chr_paralog["from"]) == closest_TE_minus1]["start"].iloc[0] - locus_len
                                        exon_candidate_area_end   = dfsv_chr_paralog[abs(dfsv_chr_paralog["from"]) == closest_TE_minus1]["start"].iloc[0]
                                        # ---- ここから追加（def なし, 隣接パラログで探索窓クランプ in-place） ----
                                        cur_name = name_unique[n_u_m]
                                        _parts = cur_name.split("_")
                                        _is_paralog = (len(_parts) >= 2 and _parts[-2] == "paralog" and _parts[-1].isdigit())
                                        if _is_paralog:
                                            _group_prefix = "_".join(_parts[:-1])
                                            _blocks = (dfsv_chr[dfsv_chr["name"].str.startswith(_group_prefix + "_")]
                                                       .groupby("name", as_index=False)
                                                       .agg(block_start=("start", "min"), block_end=("end", "max"))
                                                       .sort_values("block_start")
                                                       .reset_index(drop=True))
                                            _idx_list = _blocks.index[_blocks["name"] == cur_name].tolist()
                                            if len(_idx_list) > 0:
                                                _i = int(_idx_list[0])
                                                if _i == 0:
                                                    _left_limit = int(_blocks.loc[_i, "block_start"] - locus_len)
                                                else:
                                                    _left_limit = int(_blocks.loc[_i - 1, "block_end"])
                                                if _i == len(_blocks) - 1:
                                                    _right_limit = int(_blocks.loc[_i, "block_end"] + locus_len)
                                                else:
                                                    _right_limit = int(_blocks.loc[_i + 1, "block_start"])
                                                _eca_start = max(exon_candidate_area_start, _left_limit)
                                                _eca_end   = min(exon_candidate_area_end, _right_limit)
                                                if _eca_start < _eca_end:
                                                    exon_candidate_area_start = _eca_start
                                                    exon_candidate_area_end   = _eca_end
                                        # ---- ここまで追加 ----

                                        dfsv_chr_paralog_eca = TE_dfsv[(TE_dfsv["subject acc.ver"] == chr_ndarray[v]) & (TE_dfsv["inversion"] == "-") & (TE_dfsv["start"] > exon_candidate_area_start) & (TE_dfsv["end"] < exon_candidate_area_end)]
                                    #print("l-1087 dfsv_chr_paralog_eca is ......")
                                    #print(dfsv_chr_paralog_eca, dfsv_chr_paralog_eca.empty)
                                    #print("l-1088 dfsv_chr_paralog is ......")
                                    #print(dfsv_chr_paralog)
                                    if dfsv_chr_paralog_eca.empty == False: # 記述していないが、空だったら処理を終える
                                        #max_bit_score = dfsv_chr_paralog["bit score"].max().index()
                                        max_bit_score_idx = dfsv_chr_paralog_eca["bit score"].idxmax()
                                        #print(max_bit_score_idx)
                                        mbsi_check = list(dfsv_chr_paralog_eca["bit score"][dfsv_chr_paralog_eca["bit score"] == dfsv_chr_paralog_eca["bit score"].max()].index)
                                        #if len(mbsi_check) > 1:
                                            #max_length_idx =
                                        #print(dfsv_chr_paralog_eca["bit score"].dtype)
                                        #print(type(max_bit_score))
                                        #True_exon = dfsv_chr_paralog_eca[dfsv_chr_paralog_eca["bit score"] == max_bit_score]
                                        True_exon = dfsv_chr_paralog_eca.loc[max_bit_score_idx]
                                        if True_exon.name not in Te_check_box:
                                            Te_check_box.append(True_exon.name)
                                            if "_single_" in name_unique[n_u_m]:
                                                dfsv_chr_paralog["from"] = (-1) * dfsv_chr_paralog["from"]
                                            True_exon["name"] = name_unique[n_u_m]
                                            #print(True_exon)

                                            tmp_dcp_concat = pd.concat([dfsv_chr_paralog, True_exon.set_axis(dfsv_chr_paralog.columns).to_frame().T]).sort_index()
                                            dfsv_chr_paralog = tmp_dcp_concat
                                            #print("l-1093", dfsv_chr_paralog)
                                            #print("l-1094", dfsv_chr)
                                            dfsv_chr.update(dfsv_chr_paralog)
                                            dfsv_chr = pd.concat([dfsv_chr, dfsv_chr_paralog[~dfsv_chr_paralog.index.isin(dfsv_chr.index)]], axis=0)

                                            #tmp_dc_concat = pd.concat([dfsv_chr, dfsv_chr_paralog]).sort_index()
                                            #dfsv_chr = tmp_dc_concat
                                            #print("l-1100", dfsv_chr)

                                            #tmp_dcp_concat = pd.concat([dfsv_chr_paralog, True_exon.set_axis(dfsv_chr_paralog.columns).to_frame().T]).sort_index()
                                            #dfsv_chr_paralog = tmp_dcp_concat
                                            #print(dfsv_chr_paralog)

                                            #tmp_dc_concat = pd.concat([dfsv_chr, True_exon.set_axis(dfsv_chr.columns).to_frame().T]).sort_index()
                                            #dfsv_chr = tmp_dc_concat

                                            #tmp_dfsv_concat = pd.concat([dfsv, True_exon.set_axis(dfsv.columns).to_frame().T]).sort_index()
                                            #dfsv = tmp_dfsv_concat
                                            #print(dfsv)


                                else: # TE疑惑エキソンが中間エキソンである場合
                                    closest_TE_minus1 = min(from_list_paralog, key=lambda x: abs(x - (list_100overexons[0] - 100000000 - 1)))
                                    closest_TE_plus1  = min(from_list_paralog, key=lambda x: abs(x - (list_100overexons[0] - 100000000 + 1)))
                                    print(closest_TE_minus1)
                                    print(closest_TE_plus1)

                                    print(dfsv_chr_paralog.loc[dfsv_chr_paralog.index.min(), "inversion"])
                                    if dfsv_chr_paralog.loc[dfsv_chr_paralog.index.min(), "inversion"] == "+":
                                        exon_candidate_area_start = dfsv_chr_paralog[abs(dfsv_chr_paralog["from"]) == closest_TE_minus1]["end"].iloc[0]
                                        exon_candidate_area_end   = dfsv_chr_paralog[abs(dfsv_chr_paralog["from"]) == closest_TE_plus1]["start"].iloc[0]
                                        # ---- ここから追加（def なし, 隣接パラログで探索窓クランプ in-place） ----
                                        cur_name = name_unique[n_u_m]
                                        _parts = cur_name.split("_")
                                        _is_paralog = (len(_parts) >= 2 and _parts[-2] == "paralog" and _parts[-1].isdigit())
                                        if _is_paralog:
                                            _group_prefix = "_".join(_parts[:-1])
                                            _blocks = (dfsv_chr[dfsv_chr["name"].str.startswith(_group_prefix + "_")]
                                                       .groupby("name", as_index=False)
                                                       .agg(block_start=("start", "min"), block_end=("end", "max"))
                                                       .sort_values("block_start")
                                                       .reset_index(drop=True))
                                            _idx_list = _blocks.index[_blocks["name"] == cur_name].tolist()
                                            if len(_idx_list) > 0:
                                                _i = int(_idx_list[0])
                                                if _i == 0:
                                                    _left_limit = int(_blocks.loc[_i, "block_start"] - locus_len)
                                                else:
                                                    _left_limit = int(_blocks.loc[_i - 1, "block_end"])
                                                if _i == len(_blocks) - 1:
                                                    _right_limit = int(_blocks.loc[_i, "block_end"] + locus_len)
                                                else:
                                                    _right_limit = int(_blocks.loc[_i + 1, "block_start"])
                                                _eca_start = max(exon_candidate_area_start, _left_limit)
                                                _eca_end   = min(exon_candidate_area_end, _right_limit)
                                                if _eca_start < _eca_end:
                                                    exon_candidate_area_start = _eca_start
                                                    exon_candidate_area_end   = _eca_end
                                        # ---- ここまで追加 ----

                                        dfsv_chr_paralog_eca = TE_dfsv[(TE_dfsv["subject acc.ver"] == chr_ndarray[v]) & (TE_dfsv["inversion"] == "+") & (TE_dfsv["start"] > exon_candidate_area_start) & (TE_dfsv["end"] < exon_candidate_area_end)]
                                    else:
                                        exon_candidate_area_start = dfsv_chr_paralog[abs(dfsv_chr_paralog["from"]) == closest_TE_plus1]["end"].iloc[0]
                                        exon_candidate_area_end   = dfsv_chr_paralog[abs(dfsv_chr_paralog["from"]) == closest_TE_minus1]["start"].iloc[0]
                                        # ---- ここから追加（def なし, 隣接パラログで探索窓クランプ in-place） ----
                                        cur_name = name_unique[n_u_m]
                                        _parts = cur_name.split("_")
                                        _is_paralog = (len(_parts) >= 2 and _parts[-2] == "paralog" and _parts[-1].isdigit())
                                        if _is_paralog:
                                            _group_prefix = "_".join(_parts[:-1])
                                            _blocks = (dfsv_chr[dfsv_chr["name"].str.startswith(_group_prefix + "_")]
                                                       .groupby("name", as_index=False)
                                                       .agg(block_start=("start", "min"), block_end=("end", "max"))
                                                       .sort_values("block_start")
                                                       .reset_index(drop=True))
                                            _idx_list = _blocks.index[_blocks["name"] == cur_name].tolist()
                                            if len(_idx_list) > 0:
                                                _i = int(_idx_list[0])
                                                if _i == 0:
                                                    _left_limit = int(_blocks.loc[_i, "block_start"] - locus_len)
                                                else:
                                                    _left_limit = int(_blocks.loc[_i - 1, "block_end"])
                                                if _i == len(_blocks) - 1:
                                                    _right_limit = int(_blocks.loc[_i, "block_end"] + locus_len)
                                                else:
                                                    _right_limit = int(_blocks.loc[_i + 1, "block_start"])
                                                _eca_start = max(exon_candidate_area_start, _left_limit)
                                                _eca_end   = min(exon_candidate_area_end, _right_limit)
                                                if _eca_start < _eca_end:
                                                    exon_candidate_area_start = _eca_start
                                                    exon_candidate_area_end   = _eca_end
                                        # ---- ここまで追加 ----

                                        print(exon_candidate_area_start)
                                        print(exon_candidate_area_end)
                                        dfsv_chr_paralog_eca = TE_dfsv[(TE_dfsv["subject acc.ver"] == chr_ndarray[v]) & (TE_dfsv["inversion"] == "-") & (TE_dfsv["start"] > exon_candidate_area_start) & (TE_dfsv["end"] < exon_candidate_area_end)]
                                    #print("dfsv_chr_paralog_eca is ......")
                                    #print(dfsv_chr_paralog_eca)
                                    #print("dfsv_chr_paralog is ......")
                                    #print(dfsv_chr_paralog)
                                    if dfsv_chr_paralog_eca.empty == False: # 記述していないが、空だったら処理を終える
                                        #max_bit_score = dfsv_chr_paralog["bit score"].max().index()
                                        max_bit_score_idx = dfsv_chr_paralog_eca["bit score"].idxmax()
                                        print(max_bit_score_idx)
                                        mbsi_check = list(dfsv_chr_paralog_eca["bit score"][dfsv_chr_paralog_eca["bit score"] == dfsv_chr_paralog_eca["bit score"].max()].index)
                                        #if len(mbsi_check) > 1:
                                            #max_length_idx =
                                        #print(dfsv_chr_paralog_eca["bit score"].dtype)
                                        #print(type(max_bit_score))
                                        #True_exon = dfsv_chr_paralog_eca[dfsv_chr_paralog_eca["bit score"] == max_bit_score]
                                        True_exon = dfsv_chr_paralog_eca.loc[max_bit_score_idx]
                                        if True_exon.name not in Te_check_box:
                                            Te_check_box.append(True_exon.name)
                                            if "_single_" in name_unique[n_u_m]:
                                                dfsv_chr_paralog["from"] = (-1) * dfsv_chr_paralog["from"]
                                            True_exon["name"] = name_unique[n_u_m]
                                            #print(True_exon)

                                            tmp_dcp_concat = pd.concat([dfsv_chr_paralog, True_exon.set_axis(dfsv_chr_paralog.columns).to_frame().T]).sort_index()
                                            dfsv_chr_paralog = tmp_dcp_concat
                                            #print("l-1137", dfsv_chr_paralog)
                                            #print("l-1138", dfsv_chr)
                                            dfsv_chr.update(dfsv_chr_paralog)
                                            dfsv_chr = pd.concat([dfsv_chr, dfsv_chr_paralog[~dfsv_chr_paralog.index.isin(dfsv_chr.index)]], axis=0)

                                            #tmp_dc_concat = pd.concat([dfsv_chr, dfsv_chr_paralog]).sort_index()
                                            #dfsv_chr = tmp_dc_concat
                                            #print("l-1144", dfsv_chr)

                                            #tmp_dcp_concat = pd.concat([dfsv_chr_paralog, True_exon.set_axis(dfsv_chr_paralog.columns).to_frame().T]).sort_index()
                                            #dfsv_chr_paralog = tmp_dcp_concat
                                            #print(dfsv_chr_paralog)

                                            #tmp_dc_concat = pd.concat([dfsv_chr, True_exon.set_axis(dfsv_chr.columns).to_frame().T]).sort_index()
                                            #dfsv_chr = tmp_dc_concat

                                            #tmp_dfsv_concat = pd.concat([dfsv, True_exon.set_axis(dfsv.columns).to_frame().T]).sort_index()
                                            #dfsv = tmp_dfsv_concat
                                            #print(dfsv)



                            else: # 100回以上重複しているTE疑惑エキソンが2個以上ある時
                                dfsv_chr.set_index("idx", inplace=True)
                                #print("l-1122", dfsv_chr)
                                print(list_100overexons)
                                eca_group = []
                                eca_group_box = []
                                for item in range(len(list_100overexons)): # TE疑惑エキソンをグループ化して、グループの数、グループ内エキソンの数を把握する
                                    if item == 0:
                                        eca_group.append(list_100overexons[item])

                                    elif list_100overexons[item] == list_100overexons[item-1] + 1: # TE疑惑エキソンが連続する場合は同じecaグループ
                                        eca_group.append(list_100overexons[item])
                                        if item == len(list_100overexons) - 1:
                                            eca_group_box.append(eca_group)

                                    else: # TE疑惑エキソンが連続しなくなったら、別グループにする
                                        eca_group_box.append(eca_group)
                                        eca_group_box = []
                                        eca_group.append(list_100overexons[item])
                                        if item == len(list_100overexons) - 1:
                                            eca_group_box.append(eca_group)
                                print(eca_group_box)

                                for item in range(len(eca_group_box)): # TE疑惑エキソングループごとに処理する
                                    if len(eca_group_box[item]) == 1: # グループ内エキソンが1個しかない場合

                                        if eca_group_box[item][0] < from_list_paralog[0]: # TE疑惑エキソンが先頭エキソンである場合
                                            closest_TE_plus1  = min(from_list_paralog, key=lambda x: abs(x - (eca_group_box[item][0] - 100000000 + 1)))

                                            print(dfsv_chr_paralog.loc[dfsv_chr_paralog.index.min(), "inversion"])
                                            if dfsv_chr_paralog.loc[dfsv_chr_paralog.index.min(), "inversion"] == "+": # 左端が開放なので、locus_lenだけ手前からということにする
                                                exon_candidate_area_start = dfsv_chr_paralog[abs(dfsv_chr_paralog["from"]) == closest_TE_plus1]["start"].iloc[0] - locus_len
                                                exon_candidate_area_end   = dfsv_chr_paralog[abs(dfsv_chr_paralog["from"]) == closest_TE_plus1]["start"].iloc[0]
                                                # ---- ここから追加（def なし, 隣接パラログで探索窓クランプ in-place） ----
                                                cur_name = name_unique[n_u_m]
                                                _parts = cur_name.split("_")
                                                _is_paralog = (len(_parts) >= 2 and _parts[-2] == "paralog" and _parts[-1].isdigit())
                                                if _is_paralog:
                                                    _group_prefix = "_".join(_parts[:-1])
                                                    _blocks = (dfsv_chr[dfsv_chr["name"].str.startswith(_group_prefix + "_")]
                                                               .groupby("name", as_index=False)
                                                               .agg(block_start=("start", "min"), block_end=("end", "max"))
                                                               .sort_values("block_start")
                                                               .reset_index(drop=True))
                                                    _idx_list = _blocks.index[_blocks["name"] == cur_name].tolist()
                                                    if len(_idx_list) > 0:
                                                        _i = int(_idx_list[0])
                                                        if _i == 0:
                                                            _left_limit = int(_blocks.loc[_i, "block_start"] - locus_len)
                                                        else:
                                                            _left_limit = int(_blocks.loc[_i - 1, "block_end"])
                                                        if _i == len(_blocks) - 1:
                                                            _right_limit = int(_blocks.loc[_i, "block_end"] + locus_len)
                                                        else:
                                                            _right_limit = int(_blocks.loc[_i + 1, "block_start"])
                                                        _eca_start = max(exon_candidate_area_start, _left_limit)
                                                        _eca_end   = min(exon_candidate_area_end, _right_limit)
                                                        if _eca_start < _eca_end:
                                                            exon_candidate_area_start = _eca_start
                                                            exon_candidate_area_end   = _eca_end
                                                # ---- ここまで追加 ----

                                                dfsv_chr_paralog_eca = TE_dfsv[(TE_dfsv["subject acc.ver"] == chr_ndarray[v]) & (TE_dfsv["inversion"] == "+") & (TE_dfsv["start"] > exon_candidate_area_start) & (TE_dfsv["end"] < exon_candidate_area_end)]
                                            else: # 右端が開放なので、locus_lenだけ奥まで、ということにする
                                                exon_candidate_area_start = dfsv_chr_paralog[abs(dfsv_chr_paralog["from"]) == closest_TE_plus1]["end"].iloc[0]
                                                exon_candidate_area_end   = dfsv_chr_paralog[abs(dfsv_chr_paralog["from"]) == closest_TE_plus1]["end"].iloc[0] + locus_len
                                                # ---- ここから追加（def なし, 隣接パラログで探索窓クランプ in-place） ----
                                                cur_name = name_unique[n_u_m]
                                                _parts = cur_name.split("_")
                                                _is_paralog = (len(_parts) >= 2 and _parts[-2] == "paralog" and _parts[-1].isdigit())
                                                if _is_paralog:
                                                    _group_prefix = "_".join(_parts[:-1])
                                                    _blocks = (dfsv_chr[dfsv_chr["name"].str.startswith(_group_prefix + "_")]
                                                               .groupby("name", as_index=False)
                                                               .agg(block_start=("start", "min"), block_end=("end", "max"))
                                                               .sort_values("block_start")
                                                               .reset_index(drop=True))
                                                    _idx_list = _blocks.index[_blocks["name"] == cur_name].tolist()
                                                    if len(_idx_list) > 0:
                                                        _i = int(_idx_list[0])
                                                        if _i == 0:
                                                            _left_limit = int(_blocks.loc[_i, "block_start"] - locus_len)
                                                        else:
                                                            _left_limit = int(_blocks.loc[_i - 1, "block_end"])
                                                        if _i == len(_blocks) - 1:
                                                            _right_limit = int(_blocks.loc[_i, "block_end"] + locus_len)
                                                        else:
                                                            _right_limit = int(_blocks.loc[_i + 1, "block_start"])
                                                        _eca_start = max(exon_candidate_area_start, _left_limit)
                                                        _eca_end   = min(exon_candidate_area_end, _right_limit)
                                                        if _eca_start < _eca_end:
                                                            exon_candidate_area_start = _eca_start
                                                            exon_candidate_area_end   = _eca_end
                                                # ---- ここまで追加 ----

                                                dfsv_chr_paralog_eca = TE_dfsv[(TE_dfsv["subject acc.ver"] == chr_ndarray[v]) & (TE_dfsv["inversion"] == "-") & (TE_dfsv["start"] > exon_candidate_area_start) & (TE_dfsv["end"] < exon_candidate_area_end)]
                                            #print("dfsv_chr_paralog_eca is ......")
                                            #print(dfsv_chr_paralog_eca)
                                            #print("dfsv_chr_paralog is ......")
                                            #print(dfsv_chr_paralog)
                                            if dfsv_chr_paralog_eca.empty == False: # 記述していないが、空だったら処理を終える
                                                #max_bit_score = dfsv_chr_paralog["bit score"].max().index()
                                                max_bit_score_idx = dfsv_chr_paralog_eca["bit score"].idxmax()
                                                print(max_bit_score_idx)
                                                mbsi_check = list(dfsv_chr_paralog_eca["bit score"][dfsv_chr_paralog_eca["bit score"] == dfsv_chr_paralog_eca["bit score"].max()].index)
                                                #if len(mbsi_check) > 1:
                                                    #max_length_idx =
                                                #print(dfsv_chr_paralog_eca["bit score"].dtype)
                                                #print(type(max_bit_score))
                                                #True_exon = dfsv_chr_paralog_eca[dfsv_chr_paralog_eca["bit score"] == max_bit_score]
                                                True_exon = dfsv_chr_paralog_eca.loc[max_bit_score_idx]
                                                if True_exon.name not in Te_check_box:
                                                    Te_check_box.append(True_exon.name)
                                                    if "_single_" in name_unique[n_u_m]:
                                                        dfsv_chr_paralog["from"] = (-1) * dfsv_chr_paralog["from"]
                                                    True_exon["name"] = name_unique[n_u_m]
                                                    #print(True_exon)

                                                    tmp_dcp_concat = pd.concat([dfsv_chr_paralog, True_exon.set_axis(dfsv_chr_paralog.columns).to_frame().T]).sort_index()
                                                    dfsv_chr_paralog = tmp_dcp_concat
                                                    #print("l-1207", dfsv_chr_paralog)
                                                    #print("l-1208", dfsv_chr)
                                                    dfsv_chr.update(dfsv_chr_paralog)
                                                    dfsv_chr = pd.concat([dfsv_chr, dfsv_chr_paralog[~dfsv_chr_paralog.index.isin(dfsv_chr.index)]], axis=0)

                                                    #tmp_dc_concat = pd.concat([dfsv_chr, dfsv_chr_paralog]).sort_index()
                                                    #dfsv_chr = tmp_dc_concat
                                                    #print("l-1214", dfsv_chr)

                                                    #tmp_dcp_concat = pd.concat([dfsv_chr_paralog, True_exon.set_axis(dfsv_chr_paralog.columns).to_frame().T]).sort_index()
                                                    #dfsv_chr_paralog = tmp_dcp_concat
                                                    #print(dfsv_chr_paralog)

                                                    #tmp_dc_concat = pd.concat([dfsv_chr, True_exon.set_axis(dfsv_chr.columns).to_frame().T]).sort_index()
                                                    #dfsv_chr = tmp_dc_concat

                                                    #tmp_dfsv_concat = pd.concat([dfsv, True_exon.set_axis(dfsv.columns).to_frame().T]).sort_index()
                                                    #dfsv = tmp_dfsv_concat
                                                    #print(dfsv)



                                        elif eca_group_box[item][0] > from_list_paralog[-1]: # TE疑惑エキソンが最終エキソンである場合
                                            closest_TE_minus1 = min(from_list_paralog, key=lambda x: abs(x - (eca_group_box[item][0] - 100000000 - 1)))
                                            print(closest_TE_minus1)

                                            print(dfsv_chr_paralog.loc[dfsv_chr_paralog.index.min(), "inversion"])
                                            if dfsv_chr_paralog.loc[dfsv_chr_paralog.index.min(), "inversion"] == "+": # 右端が開放なので、locus_lenだけ奥までということにする
                                                exon_candidate_area_start = dfsv_chr_paralog[abs(dfsv_chr_paralog["from"]) == closest_TE_minus1]["end"].iloc[0]
                                                exon_candidate_area_end   = dfsv_chr_paralog[abs(dfsv_chr_paralog["from"]) == closest_TE_minus1]["end"].iloc[0] + locus_len
                                                # ---- ここから追加（def なし, 隣接パラログで探索窓クランプ in-place） ----
                                                cur_name = name_unique[n_u_m]
                                                _parts = cur_name.split("_")
                                                _is_paralog = (len(_parts) >= 2 and _parts[-2] == "paralog" and _parts[-1].isdigit())
                                                if _is_paralog:
                                                    _group_prefix = "_".join(_parts[:-1])
                                                    _blocks = (dfsv_chr[dfsv_chr["name"].str.startswith(_group_prefix + "_")]
                                                               .groupby("name", as_index=False)
                                                               .agg(block_start=("start", "min"), block_end=("end", "max"))
                                                               .sort_values("block_start")
                                                               .reset_index(drop=True))
                                                    _idx_list = _blocks.index[_blocks["name"] == cur_name].tolist()
                                                    if len(_idx_list) > 0:
                                                        _i = int(_idx_list[0])
                                                        if _i == 0:
                                                            _left_limit = int(_blocks.loc[_i, "block_start"] - locus_len)
                                                        else:
                                                            _left_limit = int(_blocks.loc[_i - 1, "block_end"])
                                                        if _i == len(_blocks) - 1:
                                                            _right_limit = int(_blocks.loc[_i, "block_end"] + locus_len)
                                                        else:
                                                            _right_limit = int(_blocks.loc[_i + 1, "block_start"])
                                                        _eca_start = max(exon_candidate_area_start, _left_limit)
                                                        _eca_end   = min(exon_candidate_area_end, _right_limit)
                                                        if _eca_start < _eca_end:
                                                            exon_candidate_area_start = _eca_start
                                                            exon_candidate_area_end   = _eca_end
                                                # ---- ここまで追加 ----

                                                print("e_c_a_start")
                                                print(exon_candidate_area_start)
                                                print(exon_candidate_area_end)
                                                dfsv_chr_paralog_eca = TE_dfsv[(TE_dfsv["subject acc.ver"] == chr_ndarray[v]) & (TE_dfsv["inversion"] == "+") & (TE_dfsv["start"] > exon_candidate_area_start) & (TE_dfsv["end"] < exon_candidate_area_end)]
                                            else: # 左端が開放なので、locus_lenだけ手前から、ということにする
                                                exon_candidate_area_start = dfsv_chr_paralog[abs(dfsv_chr_paralog["from"]) == closest_TE_minus1]["start"].iloc[0] - locus_len
                                                exon_candidate_area_end   = dfsv_chr_paralog[abs(dfsv_chr_paralog["from"]) == closest_TE_minus1]["start"].iloc[0]
                                                # ---- ここから追加（def なし, 隣接パラログで探索窓クランプ in-place） ----
                                                cur_name = name_unique[n_u_m]
                                                _parts = cur_name.split("_")
                                                _is_paralog = (len(_parts) >= 2 and _parts[-2] == "paralog" and _parts[-1].isdigit())
                                                if _is_paralog:
                                                    _group_prefix = "_".join(_parts[:-1])
                                                    _blocks = (dfsv_chr[dfsv_chr["name"].str.startswith(_group_prefix + "_")]
                                                               .groupby("name", as_index=False)
                                                               .agg(block_start=("start", "min"), block_end=("end", "max"))
                                                               .sort_values("block_start")
                                                               .reset_index(drop=True))
                                                    _idx_list = _blocks.index[_blocks["name"] == cur_name].tolist()
                                                    if len(_idx_list) > 0:
                                                        _i = int(_idx_list[0])
                                                        if _i == 0:
                                                            _left_limit = int(_blocks.loc[_i, "block_start"] - locus_len)
                                                        else:
                                                            _left_limit = int(_blocks.loc[_i - 1, "block_end"])
                                                        if _i == len(_blocks) - 1:
                                                            _right_limit = int(_blocks.loc[_i, "block_end"] + locus_len)
                                                        else:
                                                            _right_limit = int(_blocks.loc[_i + 1, "block_start"])
                                                        _eca_start = max(exon_candidate_area_start, _left_limit)
                                                        _eca_end   = min(exon_candidate_area_end, _right_limit)
                                                        if _eca_start < _eca_end:
                                                            exon_candidate_area_start = _eca_start
                                                            exon_candidate_area_end   = _eca_end
                                                # ---- ここまで追加 ----

                                                dfsv_chr_paralog_eca = TE_dfsv[(TE_dfsv["subject acc.ver"] == chr_ndarray[v]) & (TE_dfsv["inversion"] == "-") & (TE_dfsv["start"] > exon_candidate_area_start) & (TE_dfsv["end"] < exon_candidate_area_end)]
                                            #print("dfsv_chr_paralog_eca is ......")
                                            #print(dfsv_chr_paralog_eca)
                                            #print("dfsv_chr_paralog is ......")
                                            #print(dfsv_chr_paralog)
                                            if dfsv_chr_paralog_eca.empty == False: # 記述していないが、空だったら処理を終える
                                                #max_bit_score = dfsv_chr_paralog["bit score"].max().index()
                                                max_bit_score_idx = dfsv_chr_paralog_eca["bit score"].idxmax()
                                                print(max_bit_score_idx)
                                                mbsi_check = list(dfsv_chr_paralog_eca["bit score"][dfsv_chr_paralog_eca["bit score"] == dfsv_chr_paralog_eca["bit score"].max()].index)
                                                #if len(mbsi_check) > 1:
                                                    #max_length_idx =
                                                #print(dfsv_chr_paralog_eca["bit score"].dtype)
                                                #print(type(max_bit_score))
                                                #True_exon = dfsv_chr_paralog_eca[dfsv_chr_paralog_eca["bit score"] == max_bit_score]
                                                True_exon = dfsv_chr_paralog_eca.loc[max_bit_score_idx]
                                                if True_exon.name not in Te_check_box:
                                                    Te_check_box.append(True_exon.name)
                                                    if "_single_" in name_unique[n_u_m]:
                                                        dfsv_chr_paralog["from"] = (-1) * dfsv_chr_paralog["from"]
                                                    True_exon["name"] = name_unique[n_u_m]
                                                    #print(True_exon)

                                                    tmp_dcp_concat = pd.concat([dfsv_chr_paralog, True_exon.set_axis(dfsv_chr_paralog.columns).to_frame().T]).sort_index()
                                                    dfsv_chr_paralog = tmp_dcp_concat
                                                    #print("l-1330", dfsv_chr_paralog)
                                                    #print("l-1331", dfsv_chr)
                                                    dfsv_chr.update(dfsv_chr_paralog)
                                                    dfsv_chr = pd.concat([dfsv_chr, dfsv_chr_paralog[~dfsv_chr_paralog.index.isin(dfsv_chr.index)]], axis=0)

                                                    #tmp_dc_concat = pd.concat([dfsv_chr, dfsv_chr_paralog]).sort_index()
                                                    #dfsv_chr = tmp_dc_concat
                                                    #print("l-1337", dfsv_chr)

                                                    #tmp_dcp_concat = pd.concat([dfsv_chr_paralog, True_exon.set_axis(dfsv_chr_paralog.columns).to_frame().T]).sort_index()
                                                    #dfsv_chr_paralog = tmp_dcp_concat
                                                    #print(dfsv_chr_paralog)

                                                    #tmp_dc_concat = pd.concat([dfsv_chr, True_exon.set_axis(dfsv_chr.columns).to_frame().T]).sort_index()
                                                    #dfsv_chr = tmp_dc_concat

                                                    #tmp_dfsv_concat = pd.concat([dfsv, True_exon.set_axis(dfsv.columns).to_frame().T]).sort_index()
                                                    #dfsv = tmp_dfsv_concat
                                                    #print("l-1236", dfsv)



                                        else: # TE疑惑エキソンが中間エキソンである場合
                                            closest_TE_minus1 = min(from_list_paralog, key=lambda x: abs(x - (eca_group_box[item][0] - 100000000 - 1)))
                                            closest_TE_plus1  = min(from_list_paralog, key=lambda x: abs(x - (eca_group_box[item][0] - 100000000 + 1)))

                                            print(dfsv_chr_paralog.loc[dfsv_chr_paralog.index.min(), "inversion"])
                                            if dfsv_chr_paralog.loc[dfsv_chr_paralog.index.min(), "inversion"] == "+":
                                                exon_candidate_area_start = dfsv_chr_paralog[abs(dfsv_chr_paralog["from"]) == closest_TE_minus1]["end"].iloc[0]
                                                exon_candidate_area_end   = dfsv_chr_paralog[abs(dfsv_chr_paralog["from"]) == closest_TE_plus1]["start"].iloc[0]
                                                # ---- ここから追加（def なし, 隣接パラログで探索窓クランプ in-place） ----
                                                cur_name = name_unique[n_u_m]
                                                _parts = cur_name.split("_")
                                                _is_paralog = (len(_parts) >= 2 and _parts[-2] == "paralog" and _parts[-1].isdigit())
                                                if _is_paralog:
                                                    _group_prefix = "_".join(_parts[:-1])
                                                    _blocks = (dfsv_chr[dfsv_chr["name"].str.startswith(_group_prefix + "_")]
                                                               .groupby("name", as_index=False)
                                                               .agg(block_start=("start", "min"), block_end=("end", "max"))
                                                               .sort_values("block_start")
                                                               .reset_index(drop=True))
                                                    _idx_list = _blocks.index[_blocks["name"] == cur_name].tolist()
                                                    if len(_idx_list) > 0:
                                                        _i = int(_idx_list[0])
                                                        if _i == 0:
                                                            _left_limit = int(_blocks.loc[_i, "block_start"] - locus_len)
                                                        else:
                                                            _left_limit = int(_blocks.loc[_i - 1, "block_end"])
                                                        if _i == len(_blocks) - 1:
                                                            _right_limit = int(_blocks.loc[_i, "block_end"] + locus_len)
                                                        else:
                                                            _right_limit = int(_blocks.loc[_i + 1, "block_start"])
                                                        _eca_start = max(exon_candidate_area_start, _left_limit)
                                                        _eca_end   = min(exon_candidate_area_end, _right_limit)
                                                        if _eca_start < _eca_end:
                                                            exon_candidate_area_start = _eca_start
                                                            exon_candidate_area_end   = _eca_end
                                                # ---- ここまで追加 ----

                                                dfsv_chr_paralog_eca = TE_dfsv[(TE_dfsv["subject acc.ver"] == chr_ndarray[v]) & (TE_dfsv["inversion"] == "+") & (TE_dfsv["start"] > exon_candidate_area_start) & (TE_dfsv["end"] < exon_candidate_area_end)]
                                            else:
                                                exon_candidate_area_start = dfsv_chr_paralog[abs(dfsv_chr_paralog["from"]) == closest_TE_plus1]["end"].iloc[0]
                                                exon_candidate_area_end   = dfsv_chr_paralog[abs(dfsv_chr_paralog["from"]) == closest_TE_minus1]["start"].iloc[0]
                                                # ---- ここから追加（def なし, 隣接パラログで探索窓クランプ in-place） ----
                                                cur_name = name_unique[n_u_m]
                                                _parts = cur_name.split("_")
                                                _is_paralog = (len(_parts) >= 2 and _parts[-2] == "paralog" and _parts[-1].isdigit())
                                                if _is_paralog:
                                                    _group_prefix = "_".join(_parts[:-1])
                                                    _blocks = (dfsv_chr[dfsv_chr["name"].str.startswith(_group_prefix + "_")]
                                                               .groupby("name", as_index=False)
                                                               .agg(block_start=("start", "min"), block_end=("end", "max"))
                                                               .sort_values("block_start")
                                                               .reset_index(drop=True))
                                                    _idx_list = _blocks.index[_blocks["name"] == cur_name].tolist()
                                                    if len(_idx_list) > 0:
                                                        _i = int(_idx_list[0])
                                                        if _i == 0:
                                                            _left_limit = int(_blocks.loc[_i, "block_start"] - locus_len)
                                                        else:
                                                            _left_limit = int(_blocks.loc[_i - 1, "block_end"])
                                                        if _i == len(_blocks) - 1:
                                                            _right_limit = int(_blocks.loc[_i, "block_end"] + locus_len)
                                                        else:
                                                            _right_limit = int(_blocks.loc[_i + 1, "block_start"])
                                                        _eca_start = max(exon_candidate_area_start, _left_limit)
                                                        _eca_end   = min(exon_candidate_area_end, _right_limit)
                                                        if _eca_start < _eca_end:
                                                            exon_candidate_area_start = _eca_start
                                                            exon_candidate_area_end   = _eca_end
                                                # ---- ここまで追加 ----

                                                print("exon_candidate_area_start", exon_candidate_area_start)
                                                print("exon_candidate_area_end", exon_candidate_area_end)
                                                dfsv_chr_paralog_eca = TE_dfsv[(TE_dfsv["subject acc.ver"] == chr_ndarray[v]) & (TE_dfsv["inversion"] == "-") & (TE_dfsv["start"] > exon_candidate_area_start) & (TE_dfsv["end"] < exon_candidate_area_end)]
                                            #print("dfsv_chr_paralog_eca is ......")
                                            #print(dfsv_chr_paralog_eca)
                                            #print("dfsv_chr_paralog is ......")
                                            #print(dfsv_chr_paralog)
                                            if dfsv_chr_paralog_eca.empty == False: # 記述していないが、空だったら処理を終える
                                                #max_bit_score = dfsv_chr_paralog["bit score"].max().index()
                                                max_bit_score_idx = dfsv_chr_paralog_eca["bit score"].idxmax()
                                                print(max_bit_score_idx)
                                                mbsi_check = list(dfsv_chr_paralog_eca["bit score"][dfsv_chr_paralog_eca["bit score"] == dfsv_chr_paralog_eca["bit score"].max()].index)
                                                #if len(mbsi_check) > 1:
                                                    #max_length_idx =
                                                #print(dfsv_chr_paralog_eca["bit score"].dtype)
                                                #print(type(max_bit_score))
                                                #True_exon = dfsv_chr_paralog_eca[dfsv_chr_paralog_eca["bit score"] == max_bit_score]
                                                True_exon = dfsv_chr_paralog_eca.loc[max_bit_score_idx]
                                                if True_exon.name not in Te_check_box:
                                                    Te_check_box.append(True_exon.name)
                                                    if "_single_" in name_unique[n_u_m]:
                                                        dfsv_chr_paralog["from"] = (-1) * dfsv_chr_paralog["from"]
                                                    True_exon["name"] = name_unique[n_u_m]
                                                    #print(True_exon)

                                                    tmp_dcp_concat = pd.concat([dfsv_chr_paralog, True_exon.set_axis(dfsv_chr_paralog.columns).to_frame().T]).sort_index()
                                                    dfsv_chr_paralog = tmp_dcp_concat
                                                    #print("l-1303", dfsv_chr_paralog)
                                                    #print("l-1304", dfsv_chr)
                                                    dfsv_chr.update(dfsv_chr_paralog)
                                                    dfsv_chr = pd.concat([dfsv_chr, dfsv_chr_paralog[~dfsv_chr_paralog.index.isin(dfsv_chr.index)]], axis=0)

                                                    #tmp_dc_concat = pd.concat([dfsv_chr, dfsv_chr_paralog]).sort_index()
                                                    #dfsv_chr = tmp_dc_concat
                                                    #print("l-1310", dfsv_chr)

                                                    #tmp_dcp_concat = pd.concat([dfsv_chr_paralog, True_exon.set_axis(dfsv_chr_paralog.columns).to_frame().T]).sort_index()
                                                    #dfsv_chr_paralog = tmp_dcp_concat
                                                    #print(dfsv_chr_paralog)

                                                    #tmp_dc_concat = pd.concat([dfsv_chr, True_exon.set_axis(dfsv_chr.columns).to_frame().T]).sort_index()
                                                    #dfsv_chr = tmp_dc_concat

                                                    #tmp_dfsv_concat = pd.concat([dfsv, True_exon.set_axis(dfsv.columns).to_frame().T]).sort_index()
                                                    #dfsv = tmp_dfsv_concat
                                                    #print(dfsv)


                                    else: # グループ内エキソンが複数個ある場合
                                        #non_exon_flag = 0
                                        egb_check_box = [0] * len(eca_group_box[item])
                                        if eca_group_box[item][-1] < from_list_paralog[0]: # TE疑惑エキソングループが先頭エキソンである場合
                                            eca_group_box[item].sort(reverse=True)
                                            for item2 in range(len(eca_group_box[item])):
                                                #if non_exon_flag == 0: # 0でない場合は記述していないが、処理終了
                                                if item2 == 0: # グループ内末尾エキソンの場合
                                                    closest_TE_plus1 = min(from_list_paralog, key=lambda x: abs(x - (eca_group_box[item][0] - 100000000 + 1)))
                                                else: # グループ内末尾エキソンではない場合は、番号の1つ大きいグループ内エキソンか、グループ直後のエキソンが該当する
                                                    if egb_check_box[item2] == 1:
                                                        closest_TE_plus1 = eca_group_box[item][item2 - 1]
                                                    else:
                                                        closest_TE_plus1 = min(from_list_paralog, key=lambda x: abs(x - (eca_group_box[item][0] - 100000000 + 1)))

                                                if dfsv_chr_paralog.loc[dfsv_chr_paralog.index.min(), "inversion"] == "+": # 左端が開放なので、locus_lenだけ手前からということにする
                                                    exon_candidate_area_start = dfsv_chr_paralog[abs(dfsv_chr_paralog["from"]) == closest_TE_plus1]["start"].iloc[0] - locus_len
                                                    exon_candidate_area_end   = dfsv_chr_paralog[abs(dfsv_chr_paralog["from"]) == closest_TE_plus1]["start"].iloc[0]
                                                    # ---- ここから追加（def なし, 隣接パラログで探索窓クランプ in-place） ----
                                                    cur_name = name_unique[n_u_m]
                                                    _parts = cur_name.split("_")
                                                    _is_paralog = (len(_parts) >= 2 and _parts[-2] == "paralog" and _parts[-1].isdigit())
                                                    if _is_paralog:
                                                        _group_prefix = "_".join(_parts[:-1])
                                                        _blocks = (dfsv_chr[dfsv_chr["name"].str.startswith(_group_prefix + "_")]
                                                                   .groupby("name", as_index=False)
                                                                   .agg(block_start=("start", "min"), block_end=("end", "max"))
                                                                   .sort_values("block_start")
                                                                   .reset_index(drop=True))
                                                        _idx_list = _blocks.index[_blocks["name"] == cur_name].tolist()
                                                        if len(_idx_list) > 0:
                                                            _i = int(_idx_list[0])
                                                            if _i == 0:
                                                                _left_limit = int(_blocks.loc[_i, "block_start"] - locus_len)
                                                            else:
                                                                _left_limit = int(_blocks.loc[_i - 1, "block_end"])
                                                            if _i == len(_blocks) - 1:
                                                                _right_limit = int(_blocks.loc[_i, "block_end"] + locus_len)
                                                            else:
                                                                _right_limit = int(_blocks.loc[_i + 1, "block_start"])
                                                            _eca_start = max(exon_candidate_area_start, _left_limit)
                                                            _eca_end   = min(exon_candidate_area_end, _right_limit)
                                                            if _eca_start < _eca_end:
                                                                exon_candidate_area_start = _eca_start
                                                                exon_candidate_area_end   = _eca_end
                                                    # ---- ここまで追加 ----

                                                    dfsv_chr_paralog_eca = TE_dfsv[(TE_dfsv["subject acc.ver"] == chr_ndarray[v]) & (TE_dfsv["inversion"] == "+") & (TE_dfsv["start"] > exon_candidate_area_start) & (TE_dfsv["end"] < exon_candidate_area_end) & (TE_dfsv["from"] == eca_group_box[item][item2])]
                                                else: # 右端が開放なので、locus_lenだけ奥まで、ということにする
                                                    exon_candidate_area_start = dfsv_chr_paralog[abs(dfsv_chr_paralog["from"]) == closest_TE_plus1]["end"].iloc[0]
                                                    exon_candidate_area_end   = dfsv_chr_paralog[abs(dfsv_chr_paralog["from"]) == closest_TE_plus1]["end"].iloc[0] + locus_len
                                                    # ---- ここから追加（def なし, 隣接パラログで探索窓クランプ in-place） ----
                                                    cur_name = name_unique[n_u_m]
                                                    _parts = cur_name.split("_")
                                                    _is_paralog = (len(_parts) >= 2 and _parts[-2] == "paralog" and _parts[-1].isdigit())
                                                    if _is_paralog:
                                                        _group_prefix = "_".join(_parts[:-1])
                                                        _blocks = (dfsv_chr[dfsv_chr["name"].str.startswith(_group_prefix + "_")]
                                                                   .groupby("name", as_index=False)
                                                                   .agg(block_start=("start", "min"), block_end=("end", "max"))
                                                                   .sort_values("block_start")
                                                                   .reset_index(drop=True))
                                                        _idx_list = _blocks.index[_blocks["name"] == cur_name].tolist()
                                                        if len(_idx_list) > 0:
                                                            _i = int(_idx_list[0])
                                                            if _i == 0:
                                                                _left_limit = int(_blocks.loc[_i, "block_start"] - locus_len)
                                                            else:
                                                                _left_limit = int(_blocks.loc[_i - 1, "block_end"])
                                                            if _i == len(_blocks) - 1:
                                                                _right_limit = int(_blocks.loc[_i, "block_end"] + locus_len)
                                                            else:
                                                                _right_limit = int(_blocks.loc[_i + 1, "block_start"])
                                                            _eca_start = max(exon_candidate_area_start, _left_limit)
                                                            _eca_end   = min(exon_candidate_area_end, _right_limit)
                                                            if _eca_start < _eca_end:
                                                                exon_candidate_area_start = _eca_start
                                                                exon_candidate_area_end   = _eca_end
                                                    # ---- ここまで追加 ----

                                                    dfsv_chr_paralog_eca = TE_dfsv[(TE_dfsv["subject acc.ver"] == chr_ndarray[v]) & (TE_dfsv["inversion"] == "-") & (TE_dfsv["start"] > exon_candidate_area_start) & (TE_dfsv["end"] < exon_candidate_area_end) & (TE_dfsv["from"] == eca_group_box[item][item2])]
                                                if dfsv_chr_paralog_eca.empty == False: # 記述していないが、空だったら処理を終える
                                                    #non_exon_flag += 1
                                                    egb_check_box[item2] += 1
                                                    #print("dfsv_chr_paralog_eca is ......")
                                                    #print(dfsv_chr_paralog_eca)
                                                    #print("dfsv_chr_paralog is ......")
                                                    #print(dfsv_chr_paralog)
                                                    #max_bit_score = dfsv_chr_paralog["bit score"].max().index()
                                                    max_bit_score_idx = dfsv_chr_paralog_eca["bit score"].idxmax()
                                                    print(max_bit_score_idx)
                                                    mbsi_check = list(dfsv_chr_paralog_eca["bit score"][dfsv_chr_paralog_eca["bit score"] == dfsv_chr_paralog_eca["bit score"].max()].index)
                                                    #if len(mbsi_check) > 1:
                                                        #max_length_idx =
                                                    #print(dfsv_chr_paralog_eca["bit score"].dtype)
                                                    #print(type(max_bit_score))
                                                    #True_exon = dfsv_chr_paralog_eca[dfsv_chr_paralog_eca["bit score"] == max_bit_score]
                                                    True_exon = dfsv_chr_paralog_eca.loc[max_bit_score_idx]
                                                    if True_exon.name not in Te_check_box:
                                                        Te_check_box.append(True_exon.name)
                                                        if "_single_" in name_unique[n_u_m]:
                                                            dfsv_chr_paralog["from"] = (-1) * dfsv_chr_paralog["from"]
                                                            #dfsv_chr_paralog[dfsv_chr_paralog["from"]>100000000]["from"] = (-1) * dfsv_chr_paralog["from"]
                                                        True_exon["name"] = name_unique[n_u_m]
                                                        #print(True_exon)

                                                        tmp_dcp_concat = pd.concat([dfsv_chr_paralog, True_exon.set_axis(dfsv_chr_paralog.columns).to_frame().T]).sort_index()
                                                        dfsv_chr_paralog = tmp_dcp_concat
                                                        #print("l-1330", dfsv_chr_paralog)
                                                        #print("l-1331", dfsv_chr)
                                                        dfsv_chr.update(dfsv_chr_paralog)
                                                        dfsv_chr = pd.concat([dfsv_chr, dfsv_chr_paralog[~dfsv_chr_paralog.index.isin(dfsv_chr.index)]], axis=0)

                                                        #tmp_dc_concat = pd.concat([dfsv_chr, dfsv_chr_paralog]).sort_index()
                                                        #dfsv_chr = tmp_dc_concat
                                                        #print("l-1337", dfsv_chr)

                                                        #tmp_dfsv_concat = pd.concat([dfsv, dfsv_chr]).sort_index()
                                                        #dfsv = tmp_dfsv_concat
                                                        #print("l-1389", dfsv)


                                        elif eca_group_box[item][0] > from_list_paralog[-1]: # TE疑惑エキソングループが最終エキソンである場合
                                            for item2 in range(len(eca_group_box[item])):
                                                #if non_exon_flag != len(eca_group_box[item]): # 0でない場合は記述していないが、処理終了
                                                if item2 == 0: # グループ内先頭エキソンの場合
                                                    closest_TE_minus1 = min(from_list_paralog, key=lambda x: abs(x - (eca_group_box[item][0] - 100000000 - 1)))
                                                    #print("l-1350", item2, closest_TE_minus1, eca_group_box[item][item2])
                                                else: # グループ内先頭エキソンではない場合は、番号の1つ小さいグループ内エキソンか、グループ直前のエキソンが該当する
                                                    #print(eca_group_box[item])
                                                    #print(eca_group_box[item][item2])
                                                    #print(eca_group_box[item][item2 - 1])
                                                    if egb_check_box[item2] == 1:
                                                        closest_TE_minus1 = eca_group_box[item][item2 - 1]
                                                    else:
                                                        closest_TE_minus1 = min(from_list_paralog, key=lambda x: abs(x - (eca_group_box[item][0] - 100000000 - 1)))
                                                #print("l-1343 closest_TE_minus1 is ......")
                                                #print(closest_TE_minus1)
                                                #print("l-1352 dfsv_chr_paralog is ......")
                                                #print(dfsv_chr_paralog)
                                                if dfsv_chr_paralog.loc[dfsv_chr_paralog.index.min(), "inversion"] == "+": # 右端が開放なので、locus_lenだけ奥までということにする
                                                    exon_candidate_area_start = dfsv_chr_paralog[abs(dfsv_chr_paralog["from"]) == closest_TE_minus1]["end"].iloc[0]
                                                    exon_candidate_area_end   = dfsv_chr_paralog[abs(dfsv_chr_paralog["from"]) == closest_TE_minus1]["end"].iloc[0] + locus_len
                                                    # ---- ここから追加（def なし, 隣接パラログで探索窓クランプ in-place） ----
                                                    cur_name = name_unique[n_u_m]
                                                    _parts = cur_name.split("_")
                                                    _is_paralog = (len(_parts) >= 2 and _parts[-2] == "paralog" and _parts[-1].isdigit())
                                                    if _is_paralog:
                                                        _group_prefix = "_".join(_parts[:-1])
                                                        _blocks = (dfsv_chr[dfsv_chr["name"].str.startswith(_group_prefix + "_")]
                                                                   .groupby("name", as_index=False)
                                                                   .agg(block_start=("start", "min"), block_end=("end", "max"))
                                                                   .sort_values("block_start")
                                                                   .reset_index(drop=True))
                                                        _idx_list = _blocks.index[_blocks["name"] == cur_name].tolist()
                                                        if len(_idx_list) > 0:
                                                            _i = int(_idx_list[0])
                                                            if _i == 0:
                                                                _left_limit = int(_blocks.loc[_i, "block_start"] - locus_len)
                                                            else:
                                                                _left_limit = int(_blocks.loc[_i - 1, "block_end"])
                                                            if _i == len(_blocks) - 1:
                                                                _right_limit = int(_blocks.loc[_i, "block_end"] + locus_len)
                                                            else:
                                                                _right_limit = int(_blocks.loc[_i + 1, "block_start"])
                                                            _eca_start = max(exon_candidate_area_start, _left_limit)
                                                            _eca_end   = min(exon_candidate_area_end, _right_limit)
                                                            if _eca_start < _eca_end:
                                                                exon_candidate_area_start = _eca_start
                                                                exon_candidate_area_end   = _eca_end
                                                    # ---- ここまで追加 ----

                                                    print("e_c_a_start")
                                                    print(exon_candidate_area_start)
                                                    print(exon_candidate_area_end)
                                                    #dfsv_chr_paralog_eca = TE_dfsv[(TE_dfsv["subject acc.ver"] == chr_ndarray[v]) & (TE_dfsv["inversion"] == "+") & (TE_dfsv["start"] > exon_candidate_area_start) & (TE_dfsv["end"] < exon_candidate_area_end)]
                                                    dfsv_chr_paralog_eca = TE_dfsv[(TE_dfsv["subject acc.ver"] == chr_ndarray[v]) & (TE_dfsv["inversion"] == "+") & (TE_dfsv["start"] > exon_candidate_area_start) & (TE_dfsv["end"] < exon_candidate_area_end) & (TE_dfsv["from"] == eca_group_box[item][item2])]
                                                else: # 左端が開放なので、locus_lenだけ手前から、ということにする
                                                    exon_candidate_area_start = dfsv_chr_paralog[abs(dfsv_chr_paralog["from"]) == closest_TE_minus1]["start"].iloc[0] - locus_len
                                                    exon_candidate_area_end   = dfsv_chr_paralog[abs(dfsv_chr_paralog["from"]) == closest_TE_minus1]["start"].iloc[0]
                                                    # ---- ここから追加（def なし, 隣接パラログで探索窓クランプ in-place） ----
                                                    cur_name = name_unique[n_u_m]
                                                    _parts = cur_name.split("_")
                                                    _is_paralog = (len(_parts) >= 2 and _parts[-2] == "paralog" and _parts[-1].isdigit())
                                                    if _is_paralog:
                                                        _group_prefix = "_".join(_parts[:-1])
                                                        _blocks = (dfsv_chr[dfsv_chr["name"].str.startswith(_group_prefix + "_")]
                                                                   .groupby("name", as_index=False)
                                                                   .agg(block_start=("start", "min"), block_end=("end", "max"))
                                                                   .sort_values("block_start")
                                                                   .reset_index(drop=True))
                                                        _idx_list = _blocks.index[_blocks["name"] == cur_name].tolist()
                                                        if len(_idx_list) > 0:
                                                            _i = int(_idx_list[0])
                                                            if _i == 0:
                                                                _left_limit = int(_blocks.loc[_i, "block_start"] - locus_len)
                                                            else:
                                                                _left_limit = int(_blocks.loc[_i - 1, "block_end"])
                                                            if _i == len(_blocks) - 1:
                                                                _right_limit = int(_blocks.loc[_i, "block_end"] + locus_len)
                                                            else:
                                                                _right_limit = int(_blocks.loc[_i + 1, "block_start"])
                                                            _eca_start = max(exon_candidate_area_start, _left_limit)
                                                            _eca_end   = min(exon_candidate_area_end, _right_limit)
                                                            if _eca_start < _eca_end:
                                                                exon_candidate_area_start = _eca_start
                                                                exon_candidate_area_end   = _eca_end
                                                    # ---- ここまで追加 ----

                                                    #dfsv_chr_paralog_eca = TE_dfsv[(TE_dfsv["subject acc.ver"] == chr_ndarray[v]) & (TE_dfsv["inversion"] == "-") & (TE_dfsv["start"] > exon_candidate_area_start) & (TE_dfsv["end"] < exon_candidate_area_end)]
                                                    dfsv_chr_paralog_eca = TE_dfsv[(TE_dfsv["subject acc.ver"] == chr_ndarray[v]) & (TE_dfsv["inversion"] == "-") & (TE_dfsv["start"] > exon_candidate_area_start) & (TE_dfsv["end"] < exon_candidate_area_end) & (TE_dfsv["from"] == eca_group_box[item][item2])]
                                                #print("dfsv_chr_paralog_eca is ......")
                                                #print(dfsv_chr_paralog_eca)
                                                if dfsv_chr_paralog_eca.empty == False: # 記述していないが、空だったら処理を終える
                                                    #non_exon_flag += 1
                                                    egb_check_box[item2] += 1
                                                    #print("l-1382 egb_check_box", egb_check_box)
                                                    #print("dfsv_chr_paralog_eca is ......")
                                                    #print(dfsv_chr_paralog_eca)
                                                    #print("dfsv_chr_paralog is ......")
                                                    #print(dfsv_chr_paralog)
                                                    #max_bit_score = dfsv_chr_paralog["bit score"].max().index()
                                                    max_bit_score_idx = dfsv_chr_paralog_eca["bit score"].idxmax()
                                                    print(max_bit_score_idx)
                                                    mbsi_check = list(dfsv_chr_paralog_eca["bit score"][dfsv_chr_paralog_eca["bit score"] == dfsv_chr_paralog_eca["bit score"].max()].index)
                                                    #if len(mbsi_check) > 1:
                                                        #max_length_idx =
                                                    #print(dfsv_chr_paralog_eca["bit score"].dtype)
                                                    #print(type(max_bit_score))
                                                    #True_exon = dfsv_chr_paralog_eca[dfsv_chr_paralog_eca["bit score"] == max_bit_score]
                                                    True_exon = dfsv_chr_paralog_eca.loc[max_bit_score_idx]
                                                    if True_exon.name not in Te_check_box:
                                                        Te_check_box.append(True_exon.name)
                                                        if "_single_" in name_unique[n_u_m]:
                                                            #dfsv_chr_paralog[dfsv_chr_paralog["from"]>100000000]["from"] = (-1) * dfsv_chr_paralog["from"]
                                                            dfsv_chr_paralog["from"] = (-1) * dfsv_chr_paralog["from"]
                                                        True_exon["name"] = name_unique[n_u_m]
                                                        #print("l-1378", True_exon)

                                                        tmp_dcp_concat = pd.concat([dfsv_chr_paralog, True_exon.set_axis(dfsv_chr_paralog.columns).to_frame().T]).sort_index()
                                                        dfsv_chr_paralog = tmp_dcp_concat
                                                        #print("l-1382", dfsv_chr_paralog)
                                                        #print("l-1383", dfsv_chr)
                                                        dfsv_chr.update(dfsv_chr_paralog)
                                                        dfsv_chr = pd.concat([dfsv_chr, dfsv_chr_paralog[~dfsv_chr_paralog.index.isin(dfsv_chr.index)]], axis=0)

                                                        #tmp_dc_concat = pd.concat([dfsv_chr, dfsv_chr_paralog]).sort_index()
                                                        #dfsv_chr = tmp_dc_concat
                                                        #print("l-1385", dfsv_chr)

                                                        #tmp_dfsv_concat = pd.concat([dfsv, dfsv_chr]).sort_index()
                                                        #dfsv = tmp_dfsv_concat
                                                        #print("l-1389", dfsv)



                                        else: # TE疑惑エキソングループが中間エキソンである場合
                                            for item2 in range(len(eca_group_box[item])):
                                                #if non_exon_flag == 0: # 0でない場合は記述していないが、処理終了
                                                if item2 == 0: # グループ内先頭エキソンの場合
                                                    closest_TE_minus1 = min(from_list_paralog, key=lambda x: abs(x - (eca_group_box[item][0] - 100000000 - 1)))
                                                else: # グループ内先頭エキソンではない場合は、番号の1つ小さいグループ内エキソンか、グループ直前のエキソンが該当する
                                                    #print(eca_group_box[item])
                                                    #print(eca_group_box[item][item2])
                                                    #print(eca_group_box[item][item2 - 1])
                                                    if egb_check_box[item2] == 1:
                                                        closest_TE_minus1 = eca_group_box[item][item2 - 1]
                                                    else:
                                                        closest_TE_minus1 = min(from_list_paralog, key=lambda x: abs(x - (eca_group_box[item][0] - 100000000 - 1)))
                                                    closest_TE_minus1 = eca_group_box[item][item2 - 1]

                                                closest_TE_plus1  = min(from_list_paralog, key=lambda x: abs(x - (eca_group_box[item][-1] - 100000000 + 1)))

                                                print(dfsv_chr_paralog.loc[dfsv_chr_paralog.index.min(), "inversion"])
                                                if dfsv_chr_paralog.loc[dfsv_chr_paralog.index.min(), "inversion"] == "+":
                                                    exon_candidate_area_start = dfsv_chr_paralog[abs(dfsv_chr_paralog["from"]) == closest_TE_minus1]["end"].iloc[0]
                                                    exon_candidate_area_end   = dfsv_chr_paralog[abs(dfsv_chr_paralog["from"]) == closest_TE_plus1]["start"].iloc[0]
                                                    # ---- ここから追加（def なし, 隣接パラログで探索窓クランプ in-place） ----
                                                    cur_name = name_unique[n_u_m]
                                                    _parts = cur_name.split("_")
                                                    _is_paralog = (len(_parts) >= 2 and _parts[-2] == "paralog" and _parts[-1].isdigit())
                                                    if _is_paralog:
                                                        _group_prefix = "_".join(_parts[:-1])
                                                        _blocks = (dfsv_chr[dfsv_chr["name"].str.startswith(_group_prefix + "_")]
                                                                   .groupby("name", as_index=False)
                                                                   .agg(block_start=("start", "min"), block_end=("end", "max"))
                                                                   .sort_values("block_start")
                                                                   .reset_index(drop=True))
                                                        _idx_list = _blocks.index[_blocks["name"] == cur_name].tolist()
                                                        if len(_idx_list) > 0:
                                                            _i = int(_idx_list[0])
                                                            if _i == 0:
                                                                _left_limit = int(_blocks.loc[_i, "block_start"] - locus_len)
                                                            else:
                                                                _left_limit = int(_blocks.loc[_i - 1, "block_end"])
                                                            if _i == len(_blocks) - 1:
                                                                _right_limit = int(_blocks.loc[_i, "block_end"] + locus_len)
                                                            else:
                                                                _right_limit = int(_blocks.loc[_i + 1, "block_start"])
                                                            _eca_start = max(exon_candidate_area_start, _left_limit)
                                                            _eca_end   = min(exon_candidate_area_end, _right_limit)
                                                            if _eca_start < _eca_end:
                                                                exon_candidate_area_start = _eca_start
                                                                exon_candidate_area_end   = _eca_end
                                                    # ---- ここまで追加 ----

                                                    dfsv_chr_paralog_eca = TE_dfsv[(TE_dfsv["subject acc.ver"] == chr_ndarray[v]) & (TE_dfsv["inversion"] == "+") & (TE_dfsv["start"] > exon_candidate_area_start) & (TE_dfsv["end"] < exon_candidate_area_end) & (TE_dfsv["from"] == eca_group_box[item][item2])]
                                                else:
                                                    exon_candidate_area_start = dfsv_chr_paralog[abs(dfsv_chr_paralog["from"]) == closest_TE_plus1]["end"].iloc[0]
                                                    exon_candidate_area_end   = dfsv_chr_paralog[abs(dfsv_chr_paralog["from"]) == closest_TE_minus1]["start"].iloc[0]
                                                    # ---- ここから追加（def なし, 隣接パラログで探索窓クランプ in-place） ----
                                                    cur_name = name_unique[n_u_m]
                                                    _parts = cur_name.split("_")
                                                    _is_paralog = (len(_parts) >= 2 and _parts[-2] == "paralog" and _parts[-1].isdigit())
                                                    if _is_paralog:
                                                        _group_prefix = "_".join(_parts[:-1])
                                                        _blocks = (dfsv_chr[dfsv_chr["name"].str.startswith(_group_prefix + "_")]
                                                                   .groupby("name", as_index=False)
                                                                   .agg(block_start=("start", "min"), block_end=("end", "max"))
                                                                   .sort_values("block_start")
                                                                   .reset_index(drop=True))
                                                        _idx_list = _blocks.index[_blocks["name"] == cur_name].tolist()
                                                        if len(_idx_list) > 0:
                                                            _i = int(_idx_list[0])
                                                            if _i == 0:
                                                                _left_limit = int(_blocks.loc[_i, "block_start"] - locus_len)
                                                            else:
                                                                _left_limit = int(_blocks.loc[_i - 1, "block_end"])
                                                            if _i == len(_blocks) - 1:
                                                                _right_limit = int(_blocks.loc[_i, "block_end"] + locus_len)
                                                            else:
                                                                _right_limit = int(_blocks.loc[_i + 1, "block_start"])
                                                            _eca_start = max(exon_candidate_area_start, _left_limit)
                                                            _eca_end   = min(exon_candidate_area_end, _right_limit)
                                                            if _eca_start < _eca_end:
                                                                exon_candidate_area_start = _eca_start
                                                                exon_candidate_area_end   = _eca_end
                                                    # ---- ここまで追加 ----

                                                    dfsv_chr_paralog_eca = TE_dfsv[(TE_dfsv["subject acc.ver"] == chr_ndarray[v]) & (TE_dfsv["inversion"] == "-") & (TE_dfsv["start"] > exon_candidate_area_start) & (TE_dfsv["end"] < exon_candidate_area_end) & (TE_dfsv["from"] == eca_group_box[item][item2])]
                                                if dfsv_chr_paralog_eca.empty == False: # 記述していないが、空だったら処理を終える
                                                    #non_exon_flag += 1
                                                    egb_check_box[item2] += 1
                                                    #print("dfsv_chr_paralog_eca is ......")
                                                    #print(dfsv_chr_paralog_eca)
                                                    #print("dfsv_chr_paralog is ......")
                                                    #print(dfsv_chr_paralog)
                                                    #max_bit_score = dfsv_chr_paralog["bit score"].max().index()
                                                    max_bit_score_idx = dfsv_chr_paralog_eca["bit score"].idxmax()
                                                    print(max_bit_score_idx)
                                                    mbsi_check = list(dfsv_chr_paralog_eca["bit score"][dfsv_chr_paralog_eca["bit score"] == dfsv_chr_paralog_eca["bit score"].max()].index)
                                                    #if len(mbsi_check) > 1:
                                                        #max_length_idx =
                                                    #print(dfsv_chr_paralog_eca["bit score"].dtype)
                                                    #print(type(max_bit_score))
                                                    #True_exon = dfsv_chr_paralog_eca[dfsv_chr_paralog_eca["bit score"] == max_bit_score]
                                                    True_exon = dfsv_chr_paralog_eca.loc[max_bit_score_idx]
                                                    if True_exon.name not in Te_check_box:
                                                        Te_check_box.append(True_exon.name)
                                                        if "_single_" in name_unique[n_u_m]:
                                                            dfsv_chr_paralog["from"] = (-1) * dfsv_chr_paralog["from"]
                                                            #dfsv_chr_paralog[dfsv_chr_paralog["from"]>100000000]["from"] = (-1) * dfsv_chr_paralog["from"]
                                                        True_exon["name"] = name_unique[n_u_m]
                                                        #print(True_exon)

                                                        tmp_dcp_concat = pd.concat([dfsv_chr_paralog, True_exon.set_axis(dfsv_chr_paralog.columns).to_frame().T]).sort_index()
                                                        dfsv_chr_paralog = tmp_dcp_concat
                                                        #print("l-1451", dfsv_chr_paralog)
                                                        #print("l-1452", dfsv_chr)
                                                        dfsv_chr.update(dfsv_chr_paralog)
                                                        dfsv_chr = pd.concat([dfsv_chr, dfsv_chr_paralog[~dfsv_chr_paralog.index.isin(dfsv_chr.index)]], axis=0)

                                                        #tmp_dc_concat = pd.concat([dfsv_chr, dfsv_chr_paralog]).sort_index()
                                                        #dfsv_chr = tmp_dc_concat
                                                        #print("l-1458", dfsv_chr)

                                                        #tmp_dfsv_concat = pd.concat([dfsv, dfsv_chr]).sort_index()
                                                        #dfsv = tmp_dfsv_concat
                                                        #print("l-1389", dfsv)

                        else:
                            dfsv_chr_paralog_eca = pd.DataFrame()
                            dfsv_chr.set_index("idx", inplace=True)
                            dfsv_chr.update(dfsv_chr_paralog)


                        if "_single_" in name_unique[n_u_m] and dfsv_chr_paralog_eca.empty == False:
                            print("@@@@@@@@@@@@@@  single  @@@@@@@@@@@")
                            print("paralog_num")
                            print(paralog_num)
                            #if paralog_num != 1: # マルチエキソンなパラログがない場合は、paralog_num=1（初期設定値）をそのまま使うのでスキップ
                            paralog_num += 1
                            pre_name = name_unique[n_u_m]
                            name_unique[n_u_m] = name_unique[n_u_m].split("_single_")[0] + "_paralog_" + str(paralog_num)
                            dfsv_chr_paralog["name"] = name_unique[n_u_m]
                            #print(pre_name)
                            #print(dfsv_chr_paralog["name"])
                            #print(dfsv_chr[dfsv_chr["name"] == pre_name]["name"])
                            dfsv_chr_part = dfsv_chr[dfsv_chr["name"] == pre_name]
                            dfsv_chr_part["name"] = name_unique[n_u_m]
                            dfsv_chr[dfsv_chr["name"] == pre_name] = dfsv_chr_part
                            #dfsv_chr[dfsv_chr["name"] == pre_name]["name"] = "^^^^^^^"#name_unique[n_u_m]
                            #print(name_unique[n_u_m])
                            #print(dfsv_chr_paralog)
                            #print(dfsv_chr)
                            #print("l-1132 paralog_num")
                            #print(paralog_num)
                            #if paralog_num == 1:
                            #    paralog_num += 1

                        if dfsv_chr_paralog.loc[dfsv_chr_paralog.index.min(), "inversion"] == "+":
                            true_rank = 1
                            for index, item in dfsv_chr_paralog.iterrows():
                                dfsv_chr_paralog.loc[index, "rank"] = true_rank
                                true_rank += 1
                        else:
                            true_rank = dfsv_chr_paralog.shape[0]
                            for index, item in dfsv_chr_paralog.iterrows():
                                dfsv_chr_paralog.loc[index, "rank"] = true_rank
                                true_rank += -1


                        #dfsv_chr.loc[dfsv_chr[dfsv_chr["name"] == name_unique[n_u_m]].index] = dfsv_chr_paralog
                        mask = (dfsv_chr["name"] == name_unique[n_u_m])
                        dfsv_chr = pd.concat([dfsv_chr.loc[~mask], dfsv_chr_paralog], axis=0)
                        print("^^^^^^^^^^^^^^done^^^^^^^^^^^^^^")
                        print(name_unique[n_u_m])
                        #print(dfsv_chr)
                        dfsv_chr_ri = dfsv_chr.reset_index()
                        dfsv_chr = dfsv_chr_ri
                        #print(dfsv_chr)
                        #if n_u_m != name_unique_num - 1:
                        dfsv_chr.rename(columns={"index":"idx"}, inplace=True)
                        print(dfsv_chr)









                else:
                    print("no TEs") # NOTCH2など


                #print("dfsv_chr is ...")
                #print(dfsv_chr)




                #df_ances.rename(columns={"index":"index1"}, inplace=True)
                #print(dfsv_chr)
                #print("l-1542", dfsv)
                dfsv_chr.set_index("idx", inplace=True)
                ## 元のデータ型を保持する
                original_dtypes = dfsv.dtypes.copy()
                #print("original_dtypes", original_dtypes)
                # ---- ここから追加：採用済みTEの 'from' を 90,000,000 減算（dfsv_chr 側, def なし） ----
                try:
                    _mask_chr_te = (dfsv_chr["from"] >= 100000000) & (dfsv_chr["name"].str.contains("_paralog_"))
                except Exception:
                    _mask_chr_te = pd.Series(False, index=dfsv_chr.index)
                if _mask_chr_te.any():
                    dfsv_chr.loc[_mask_chr_te, "from"] = dfsv_chr.loc[_mask_chr_te, "from"] - 90000000
                # ---- ここまで追加 ----
                dfsv.update(dfsv_chr)
                #print("l-1575", dfsv)
                for col, dtype in original_dtypes.items():
                    if dtype == 'int64':  # 整数型の列だけを対象
                        dfsv[col] = dfsv[col].astype('int64')
                dfsv["rank"] = pd.to_numeric(dfsv["rank"], errors='raise').astype('int64')


                #print("l-1579", dfsv)
                dfsv = pd.concat([dfsv, dfsv_chr[~dfsv_chr.index.isin(dfsv.index)]], axis=0)
                #print("l-1544", dfsv)
                dfsv = pd.concat([dfsv, dfsv_chr[~dfsv_chr.index.isin(dfsv.index)]], axis=0).sort_index()
                #print("l-1545", dfsv)
                for col, dtype in original_dtypes.items():
                    if dtype == 'int64':  # 整数型の列だけを対象
                        dfsv[col] = dfsv[col].astype('int64')
                #dfsv = dfsv.astype(int)
                #print("l-1548", dfsv)

    int_cols = ["start", "end", "ances", "repeat", "rank", "from", "q_whole_start", "q_whole_end",
                "alignment length", "mismatches", "gap opens", "q_head", "q_tail", "s_head", "s_tail"]
    dfsv[int_cols] = dfsv[int_cols].astype('Int64')


    elapsed_l1796 = time.time() - start_time
    #t_l518 = time.time()
    #t_l178_518 = t_l518 - t_l178
    dt_l1796 = datetime.datetime.now()
    print(f"l518-1796, 4GTF_startから: {elapsed_l1796:.2f} 秒, 現在: {dt_l1796}")
    #print(f"l518-1796, done: {t_l178_518:.2f} 秒, 4GTF_startから: {elapsed_l518:.2f} 秒, 現在: {dt_l518}")


    ##### dfsv.to_csv("test4_{}.tsv".format(p_args_br), sep="\t", index=False, mode="w")
    print("dfsv is ......")
    print(dfsv)
    #dups_df = dfsv[dfsv["from"] > 0]
    #dups_df.to_csv("dups_{}.tsv".format(p_args_br), sep="\t", index=False, mode="w")

#    if MODE_SV:
#        copy_sum    = dfsv[dfsv["from"] > 0]["name"].nunique()
#    else:
#        copy_sum    = 1 + dfsv[dfsv["from"] > 0]["name"].nunique()

    if TE_dfsv.empty == False:
        TE_copy_sum = len(TE_dfsv)
    else:
        TE_copy_sum = 0
    #copy_sum = dups_df["name"].nunique()

    if MODE_SV:
        if query_type == "multi-exon": # クエリがマルチエキソンの場合
            SD_copy_sum     = dfsv[dfsv["name"].str.contains("paralog")]["name"].nunique() + 1 # paralog だけだと、itselfを数え落とす、のを防ぐ
            retro_copy_sum  = dfsv[dfsv["name"].str.contains("retro")]["name"].nunique()
            single_copy_sum = dfsv[dfsv["name"].str.contains("single")]["name"].nunique()
        else: # クエリがシングルエキソンの場合は、レトロ型か、単なるシングルエキソンかを判定した上でコピー数を計数する
            SD_copy_sum = dfsv[dfsv["name"].str.contains("paralog")]["name"].nunique()
            if SD_copy_sum == 0: # SD型コピーがなければ、クエリは必ず単なるシングルエキソンであるので。
                dfsv["name"] = dfsv["name"].str.replace("retro", "single")
                retro_copy_sum  = dfsv[dfsv["name"].str.contains("retro")]["name"].nunique() # retroというタグ付けは誤りだったということになるので、付け替える
                single_copy_sum = dfsv[dfsv["name"].str.contains("single")]["name"].nunique() + 1
            else: # SD型コピーがあれば、クエリがシングルエキソンであってもレトロ型ということになるので。
                unique_names = dfsv[dfsv["name"].str.contains("paralog")]["name"].unique()
                print("l1936", unique_names)
                retro_copy_sum  = dfsv[dfsv["name"].str.contains("retro")]["name"].nunique() + 1
                single_copy_sum = dfsv[dfsv["name"].str.contains("single")]["name"].nunique()

    else:
        if query_type == "multi-exon": # クエリがマルチエキソンの場合
            SD_copy_sum     = dfsv[dfsv["name"].str.contains("paralog")]["name"].nunique()
            retro_copy_sum  = dfsv[dfsv["name"].str.contains("retro")]["name"].nunique()
            single_copy_sum = dfsv[dfsv["name"].str.contains("single")]["name"].nunique()

        else: # クエリがシングルエキソンの場合は、レトロ型か、単なるシングルエキソンかを判定した上でコピー数を計数する
            SD_copy_sum = dfsv[dfsv["name"].str.contains("paralog")]["name"].nunique()
            if df_native_align.loc[0, "from"] != 10000: # SD型コピーがなければ、クエリは必ず単なるシングルエキソンであるので。
                dfsv["name"] = dfsv["name"].str.replace("retro", "single")
                retro_copy_sum  = dfsv[dfsv["name"].str.contains("retro")]["name"].nunique() # retroというタグ付けは誤りだったということになるので、付け替える
                single_copy_sum = dfsv[dfsv["name"].str.contains("single")]["name"].nunique()
            else: # SD型コピーがあれば、クエリがシングルエキソンであってもレトロ型ということになるので。
                unique_names = dfsv[dfsv["name"].str.contains("paralog")]["name"].unique()
                #print("l1936", unique_names)
                retro_copy_sum  = dfsv[dfsv["name"].str.contains("retro")]["name"].nunique()
                single_copy_sum = dfsv[dfsv["name"].str.contains("single")]["name"].nunique()

    copy_sum = SD_copy_sum + retro_copy_sum

    #if MODE_SV == False:
    #    copy_sum += 1
    if copy_sum == 0 and single_copy_sum != 0:
        copy_sum = 1
    #copy_sum = dups_df["rank"].value_counts()[1]
    #print("copy_sum is ...")
    #print(copy_sum)

    #print('dups_df["rank"].value_counts() is ...')
    #tmp_r = dups_df["rank"].value_counts()
    #print(dups_df["rank"].value_counts())

    #print('dups_df["from"].value_counts() is ...')
    #tmp_f = dups_df["from"].value_counts()
    #print(tmp_f)

    #if 10000 in dups_df["from"].unique():
    #    retro_sum = dups_df["from"].value_counts()[10000]
    #else:
    #    retro_sum = 0
    #it_tan = copy_sum - retro_sum
    #others_num = len(dfsv) - len(dups_df)

    pd.DataFrame([{"query_transcript_ID": sss.iloc[0,0], "query_gene_name": s6split.iloc[0,1], "query_type": query_type,
                   f"{p_args_sl}_total_CN": copy_sum,
                   f"{p_args_sl}_SD_CN"   : SD_copy_sum,
                   f"{p_args_sl}_retro_CN": retro_copy_sum,
                   f"({p_args_sl}_single)": single_copy_sum,
                   f"({p_args_sl}_TE)"    : TE_copy_sum,
                   "repeat?": "No"}]).to_csv("statistics_{}_{}.tsv".format(p_args_br, s6split.iloc[0,1]), sep="\t", index=False, mode="w")

    #pd.DataFrame([{"query_transcript_ID": sss.iloc[0,0],
    #               "query_gene_name": s6split.iloc[0,1],
    #               "copy_sum": copy_sum, "itself+tandem": it_tan,
    #               "retro": retro_sum, "others(single)": others_num,
    #               "repeat?": "No"}]).to_csv("statistics_{}.tsv".format(p_args_br), sep="\t", index=False, mode="w")
    if __MODE_SV:
        if (query_type == "single-exon") and (retro_copy_sum > 0):
            dfsv.loc[dfsv["ances"] == 1, "from"] = 10000
        dfsv[dfsv["ances"] == 1].to_csv("native_locus_exon_alignments_{}_{}.tsv".format(p_args_br, s6split.iloc[0,1]), sep="\t", index=False, mode="w")

    if dfsv["from"].isin([100000000]).any():
        TE_flag = "Yes"
    if TE_flag == "Yes":
        dfsv.loc[dfsv["name"].str.contains("single"), "from"] = -1
        dfsv[dfsv["name"].str.contains("single")]["from"] = -1
    dfsv_single = dfsv[dfsv["name"].str.contains("single")]
    #dups_df = dfsv[dfsv["from"] > 0]
    dups_df = dfsv[
        ((dfsv["from"] > 0) & (dfsv["from"] < 100000000)) |
        ((dfsv["ances"] == 1) & (dfsv["from"] >= 100000000))
        ]
    dups_df.to_csv("dups_{}_{}.tsv".format(p_args_br, s6split.iloc[0,1]), sep="\t", index=False, mode="w")
    if TE_flag == "Yes":
        dfsv = dups_df

    # 重複配列が落ちている領域をlocusとして括る
    dups_list = dups_df["name"].unique().tolist()
    #if len(dups_list) <= 100:
    cols_dups = ["chr", "start", "end", "inversion", "transcript_ID", "gene_name", "area_name", "numbers_of_exons", "member_exons"]
    dups_area_concat = pd.DataFrame(index=[], columns=cols_dups)

    for dups_name in dups_list:
        df_locus = dups_df[dups_df["name"] == dups_name]
        df_locus = df_locus.reset_index(drop=True)
        num_of_ex = df_locus["rank"].max()
        #print("num_of_ex is ...")
        #print(num_of_ex)
        #print("df_locus is ...")
        #print(df_locus)
        if df_locus.loc[0, "from"] == 10000:
            dn_split = dups_name.split("_")
            dn3_split = dn_split[3].split("+")
            #print("dn3_split is ...")
            #print(dn3_split)
            #dn3_split.remove('')
            #mex_list = dn3_split
            mex_list = list(filter(None, dn3_split))
            #print("mex_list is ...")
            #print(mex_list)
            mex = '_'.join(map(str, mex_list))
        elif df_locus.loc[0, "from"] == 5000:
            mex = 'from_single_exon_query'
        else:
            mex_list = []
            for index,item in df_locus.iterrows():
                mex_list.append(df_locus.loc[index, "from"])
            mex = '_'.join(map(str, mex_list))


        dups_area_record = pd.Series([df_locus.loc[0, "subject acc.ver"],
                                      df_locus.loc[0, "start"],
                                      df_locus.loc[num_of_ex-1, "end"],
                                      df_locus.loc[0, "inversion"],
                                      sss.iloc[0,0], # transcript_ID
                                      s6split.iloc[0,1],# gene_name
                                      df_locus.loc[0, "name"],
                                      num_of_ex,
                                      mex], index=dups_area_concat.columns)
        dups_area_concat = dups_area_concat.append(dups_area_record, ignore_index=True)

        dups_area_concat.to_csv("dups_area_{}_{}.tsv".format(p_args_br, s6split.iloc[0,1]), sep="\t", index=False, mode="w")



elapsed_l1922 = time.time() - start_time
    #t_l518 = time.time()
    #t_l178_518 = t_l518 - t_l178
dt_l1922 = datetime.datetime.now()
print(f"l1796-1922, 4GTF_startから: {elapsed_l1922:.2f} 秒, 現在: {dt_l1922}")
    #print(f"l518-1796, done: {t_l178_518:.2f} 秒, 4GTF_startから: {elapsed_l518:.2f} 秒, 現在: {dt_l518}")


# GTFファイルのために並べ替えなど
dropdf = dfsv.drop(["s_head", "s_tail", "% identity", "repeat",
                    "alignment length", "mismatches", "gap opens",
                    "q_head", "q_tail", "evalue", "bit score",
                    "ances", "q_whole_start", "q_whole_end"], axis=1)

dropdf.loc[:,"subject acc.ver"] = dropdf.loc[:,"subject acc.ver"].astype(str)
#dropdf.loc[:,"subject acc.ver"] = "chr" + dropdf.loc[:,"subject acc.ver"].astype(str)

dropdf.insert(6, "source", ".")
dropdf.insert(7, "feature", "exon")
dropdf.insert(8, "score", ".")
dropdf.insert(9, "frame", ".")
dropdf.insert(10, "attribute", ".")



single_dropdf = dfsv_single.drop(["s_head", "s_tail", "% identity", "repeat",
                    "alignment length", "mismatches", "gap opens",
                    "q_head", "q_tail", "evalue", "bit score",
                    "ances", "q_whole_start", "q_whole_end"], axis=1)

single_dropdf.loc[:,"subject acc.ver"] = single_dropdf.loc[:,"subject acc.ver"].astype(str)
#dropdf.loc[:,"subject acc.ver"] = "chr" + dropdf.loc[:,"subject acc.ver"].astype(str)

single_dropdf.insert(6, "source", ".")
single_dropdf.insert(7, "feature", "exon")
single_dropdf.insert(8, "score", ".")
single_dropdf.insert(9, "frame", ".")
single_dropdf.insert(10, "attribute", ".")


# gene行・transcript行の作成・追加
ddf_name_unique = dropdf["name"].unique()
ddf_name_unique_num = dropdf["name"].nunique()
#print(ddf_name_unique)
#print(ddf_name_unique_num)
for d_n_u in range(ddf_name_unique_num):
    ddf_locus = dropdf[dropdf["name"] == ddf_name_unique[d_n_u]]
    ddf_locus = ddf_locus.sort_values(by = ["subject acc.ver", "start"], ascending = [True, True]).reset_index(drop = True)
    ddf_locus_gene = ddf_locus.iloc[[0]].copy()
    ddf_locus_gene.loc[0, "end"] = ddf_locus.loc[ddf_locus.index[-1], "end"]
    ddf_locus_gene.loc[0, "feature"] = "gene"
    ddf_locus_transcript = ddf_locus_gene.iloc[[0]].copy()
    ddf_locus_transcript.loc[0, "feature"] = "transcript"
    #print(ddf_locus)
    #print(ddf_locus_gene)
    #print(ddf_locus_transcript)
    matching_row_index = dropdf[(dropdf["subject acc.ver"] == ddf_locus_gene["subject acc.ver"].iloc[0]) &
                                (dropdf["start"] == ddf_locus_gene["start"].iloc[0])].index[0]
    dropdf_top    = dropdf.iloc[:matching_row_index]
    dropdf_bottom = dropdf.iloc[matching_row_index:]
    dropdf = pd.concat([dropdf_top, ddf_locus_gene, ddf_locus_transcript, dropdf_bottom]).reset_index(drop=True)


#print(dropdf)
elapsed_l1966 = time.time() - start_time
    #t_l518 = time.time()
    #t_l178_518 = t_l518 - t_l178
dt_l1966 = datetime.datetime.now()
print(f"l1922-1966, 4GTF_startから: {elapsed_l1966:.2f} 秒, 現在: {dt_l1966}")
    #print(f"l518-1796, done: {t_l178_518:.2f} 秒, 4GTF_startから: {elapsed_l518:.2f} 秒, 現在: {dt_l518}")

#全体シートに表を追加していく

#dropdf.to_csv("whole_arr_{}.csv".format(p_args_br), index=False, mode = "w")
gene_version = "1"
if query_type == "single-exon":
    q_type = "qSE"
else:
    q_type = "qME"

dropdf.loc[dropdf["feature"] == "gene", "attribute"] = (
'gene_id "' + dropdf.loc[dropdf["feature"] == "gene", "name"] + '"; '
'gene_version "' + gene_version + '"; '
'gene_name "' + 'BLAST_result_from_' + s6split.iloc[0 ,1] + '_' + sss.iloc[0, 0] + '_' + q_type + '"; '
'gene_source "' + 'unknown' + '"; '
'gene_biotype "' + 'unknown"'
)

dropdf.loc[dropdf["feature"] == "transcript", "attribute"] = (
'gene_id "' + dropdf.loc[dropdf["feature"] == "transcript", "name"] + '"; '
'gene_version "' + gene_version + '"; '
'gene_name "' + 'BLAST_result_from_' + s6split.iloc[0 ,1] + '_' + sss.iloc[0, 0] + '_' + q_type + '"; '
'gene_source "' + 'unknown' + '"; '
'gene_biotype "' + 'unknown' + '"; '
'transcript_id "' + dropdf.loc[dropdf["feature"] == "transcript", "name"] + '_transcript"; '
'transcript_version "' + gene_version + '"; '
'transcript_name "' + 'BLAST_result_from_' + s6split.iloc[0 ,1] + '_' + sss.iloc[0, 0] + '_' + q_type + '_transcript"; '
'transcript_source "' + 'unknown' + '"; '
'transcript_biotype "' + 'unknown"'
)

dropdf.loc[dropdf["feature"] == "exon", "attribute"] = (
'gene_id "' + dropdf.loc[dropdf["feature"] == "exon", "name"] + '"; '
'gene_version "' + gene_version + '"; '
'gene_name "' + 'BLAST_result_from_' + s6split.iloc[0 ,1] + '_' + sss.iloc[0, 0] + '_' + q_type + '"; '
'gene_source "' + 'unknown' + '"; '
'gene_biotype "' + 'unknown' + '"; '
'transcript_id "' + dropdf.loc[dropdf["feature"] == "exon", "name"] + '_transcript"; '
'transcript_version "' + gene_version + '"; '
'transcript_name "' + 'BLAST_result_from_' + s6split.iloc[0 ,1] + '_' + sss.iloc[0, 0] + '_' + q_type + '_transcript"; '
'transcript_source "' + 'unknown' + '"; '
'transcript_biotype "' + 'unknown"' + '"; '
'exon_number "' + dropdf["rank"].astype(str) + '"; '
'exon_id "' + dropdf["name"] + '_' + dropdf["rank"].astype(str) + '"; '
'exon_version "' + gene_version + '"'
)


print("l-1880 dropdf", dropdf)
#print(dropdf["attribute"])

elapsed_l2023 = time.time() - start_time
    #t_l518 = time.time()
    #t_l178_518 = t_l518 - t_l178
dt_l2023 = datetime.datetime.now()
print(f"l1966-2023, 4GTF_startから: {elapsed_l2023:.2f} 秒, 現在: {dt_l2023}")
    #print(f"l518-1796, done: {t_l178_518:.2f} 秒, 4GTF_startから: {elapsed_l518:.2f} 秒, 現在: {dt_l518}")


#dropdf["attribute"] = 'gene_id "' + dropdf["name"] + \
#                  '"; transcript_id "' + dropdf["name"] + \
#                  '"; exon_number "' + dropdf["rank"].astype("str") + \
#                  '"; exon_id "' + dropdf["name"] + '_' + dropdf["rank"].astype("str") + \
#                  '"; gene_name "' + s6split.iloc[:,1] + '_family"'
#for index,item in dropdf.iterrows():
#    s = dropdf.loc[index,"query acc.ver"].split("|")
#    dropdf.loc[index,"attribute"] = dropdf.loc[index, "attribute"].format(s[0],s[3],dropdf.loc[index,"name"])
dr2df = dropdf.drop("query acc.ver", axis=1)
rdf = dropdf.reindex(columns=["subject acc.ver", "source", "feature",
                            "start", "end", "score", "inversion",
                            "frame", "attribute"])

rdf = rdf.reset_index(drop=True)


single_dr2df = single_dropdf.drop("query acc.ver", axis=1)
srdf = single_dropdf.reindex(columns=["subject acc.ver", "source", "feature",
                            "start", "end", "score", "inversion",
                            "frame", "attribute"])

srdf = srdf.reset_index(drop=True)


rdf.to_csv("GTF_from_BLASTresult_{}_{}.gtf".format(p_args_br, s6split.iloc[0,1]), sep="\t", header=False, index=False, mode = "w", quoting=csv.QUOTE_NONE)
srdf.to_csv("single_GTF_from_BLASTresult_{}_{}.gtf".format(p_args_br, s6split.iloc[0,1]), sep="\t", header=False, index=False, mode = "w", quoting=csv.QUOTE_NONE)

# GTFファイルの作成
#import pathlib
#p_file = pathlib.Path("whole_arr_{}.tsv".format(p_args_br))
#p_file.with_suffix(".gtf")


#forGTF(p_args_br)

# if __name__ == "__main__":
#    p = Pool(mp.cpu_count())
#    p.map(forGTF, tl_list[:-1])
#    p.close()


#statistic_files = glob.glob("statistics_*")

##for a in statistic_files:
##    print(a)

#sf_data_list = []

#for file in statistic_files:
#    sf_data_list.append(pd.read_csv(file))

#sf_cat_df = pd.concat(sf_data_list, axis=0, sort=True)

#sf_cat_df.to_csv("total_stat_2HS.csv",index=False)



#dups_area_files = glob.glob("dups_area_*")

##for a in statistic_files:
##   print(a)

#da_data_list = []

#for file in dups_area_files:
#    da_data_list.append(pd.read_csv(file))

#da_cat_df = pd.concat(da_data_list, axis=0, sort=True)

#da_cat_df.to_csv("total_dups_area_2HS.csv",index=False)
