
import numpy as np
import pandas as pd
import gmpy2
from gmpy2 import comb
from itertools import chain
import json
import math 
from multiprocessing import Pool
import sys

# =====================================================================
#
# エンリッチメント解析のために必要なデータの読込
#
# -----------------------------

args = sys.argv

# 除外すべきpathwayのセットの読込                                                                                       
excluded_pathway = set(pd.read_csv('./data/excluded_pathways_ko2.tsv', sep='\t',header=None)[0])

#KEGGデータを選択したときのpathwayデータの読み込み
if args[1] == 'KEGG':
    # pathwayに対応するproteinのデータの読込
    pathway = pd.read_csv('./data/ko210706.tsv',sep='\t')
    
    # 使用するpathwayのみ抽出
    pathway_use = set(pathway['pathway'])-excluded_pathway
    pathway = pathway[pathway['pathway'].isin(pathway_use)]

    # pathwayとproteinを変換するためにdict
    pathway2protein = pathway.groupby('pathway')['protein'].apply(list)
    pathway2protein = dict(zip(pathway2protein.index, pathway2protein))

    protein2pathway = pathway.groupby('protein')['pathway'].apply(list)
    protein2pathway = dict(zip(protein2pathway.index, protein2pathway))

    # hsaからkoに変換するデータの読込
    hsa_ko = pd.read_csv('./data/hsa_ko191007.list', sep='\t', header=None)
    hsa2ko_dict = hsa_ko.groupby(0)[1].apply(list).to_dict()

    # すべての遺伝子の情報の読み込み
    All_genes = pd.read_csv('./data/ko210706.tsv', sep='\t')
    All_genes = set(All_genes[~All_genes['pathway'].isin(excluded_pathway)]['protein'])

    #pathwayに含まれるタンパク質データ    
    df_pathway_in_hsa = pd.read_table('./data/pathway_hsa_list.tsv',sep='\t')

    if args[2] == 'English':
        #pathway詳細データ　[pathway, name_ja, name_en, symptom, KEGG URL]                                                          
        df_pathway = pd.read_table('./data/all_pathway_ID_name.tsv',sep='\t')

    if args[2] == 'Japanese':
        #pathway詳細データ　[pathway, name_ja, name_en, symptom, KEGG URL]
        df_pathway = pd.read_table('./data/all_pathway_ID_name_ja.tsv',sep='\t')
        

#UMLSデータを選択したときのdiseaseデータの読み込み
if args[1] == 'UMLS':
    # 除外すべきpathwayのセットの読込
    excluded_pathway = set(pd.read_csv('./data/excluded_pathways_ko2.tsv', sep='\t',header=None)[0])

    pathway = pd.read_csv('./data/Disease_protein_ootani.tsv',sep='\t')

    # diseaseとproteinを変換するためにdict
    pathway2protein = pathway.groupby('pathway')['protein'].apply(list)
    pathway2protein = dict(zip(pathway2protein.index, pathway2protein))
    protein2pathway = pathway.groupby('protein')['pathway'].apply(list)
    protein2pathway = dict(zip(protein2pathway.index, protein2pathway))

    All_genes = pd.read_csv('./data/Disease_protein_ootani.tsv', sep='\t')
    All_genes = set(All_genes[~All_genes['pathway'].isin(excluded_pathway)]['protein'])

    #diseaseに含まれるタンパク質データ                                                                                      
    df_pathway_in_hsa = pd.read_table('./data/Disease_protein_ootani_hsa.tsv',sep='\t')

    if args[2] == 'English':
        #pathway詳細データ　[pathway, name_ja, name_en, symptom, KEGG URL]                                                          
        df_pathway = pd.read_table('./data/all_disease_ID_name2.tsv',sep='\t')

    if args[2] == 'Japanese':
        #pathway詳細データ　[pathway, name_ja, name_en, symptom, KEGG URL]
        df_pathway = pd.read_table('./data/all_disease_ID_name_ja2.tsv',sep='\t')

#Englishを選択したときの漢方薬、生薬、構成化合物等のデータ読み込み
if args[2] == 'English':
    if args[3] == 'Kampo':
        #生薬-タンパク質スコア {生薬名:{protein1:score1, protein2:score2}}
        infile = './data/kampo_protein_score.json'
        
        #既知の生薬-タンパク質 {生薬名:[protein1,protein2]}
        infile2 = './data/positive_kampo_protein.json'

        #漢方薬-タンパク質スコアの読み込み
        with open(infile) as wf:
            kampo_protein = json.load(wf)

        #既知の漢方薬-タンパク質の読み込み
        with open(infile2) as f:
            known_kampo_protein = json.load(f)

    if args[3] == 'Custom':
        #生薬-タンパク質スコア {生薬名:{protein1:score1, protein2:score2}}
        infile = './data/crude_protein_score_english.json'

        #既知の生薬-タンパク質 {生薬名:[protein1,protein2]}
        infile2 = './data/positive_crude_protein_english.json'

        #生薬-タンパク質スコアの読み込み
        with open(infile) as wf:
            crude_protein = json.load(wf)

        #既知の生薬-タンパク質の読み込み
        with open(infile2) as f:
            known_crude_protein = json.load(f)
        
    #生薬の詳細データ　[crude, kampo, compound, kegg_link, id] 
    df_crude = pd.read_table('./data/all_crude_kampo_compound_english5.tsv',sep='\t')
    
    #化合物詳細データ　[2797化合物, 6レコードでInChI情報なし, link, crude, kampo, Name] 
    df_compound = pd.read_table('./data/all_compound_information.tsv',sep='\t')

    #漢方薬に含まれる生薬データ
    df_kampo_crude = pd.read_table('./data/example_194_kampo_en.tsv',sep='\t')

    #新規適応疾患(正規化のスコア0.5以上)
    input_disease = pd.read_table('./data/Differential_from_kp_matrix5_for_heatmap_normalized_and_scaled_converted.txt',sep='\t')

    #漢方の詳細データ　[KNApSAcK ID, Multi-herbal medicine name, based_link, link, Multi-herbal medicine name japanese]     
    df_kampo = pd.read_table('./data/Kampo_KNApSAcK_ID2.tsv',sep='\t')


#Japaneseを選択したときの漢方薬、生薬、構成化合物等のデータ読み込み
if args[2] == 'Japanese':
    if args[3] == 'Kampo':
        #生薬-タンパク質スコア {生薬名:{protein1:score1, protein2:score2}}                                                  
        infile = './data/kampo_protein_score_ja.json'
	
        #既知の生薬-タンパク質 {生薬名:[protein1,protein2]}                                                                 
        infile2 = './data/positive_kampo_protein_ja.json'

        #漢方薬-タンパク質スコアの読み込み
        with open(infile) as wf:
            kampo_protein = json.load(wf)

        #既知の漢方薬-タンパク質の読み込み
        with open(infile2) as f:
            known_kampo_protein = json.load(f)

    if args[3] == 'Custom':
        #生薬-タンパク質スコア {生薬名:{protein1:score1, protein2:score2}}                                                        
        infile = './data/crude_protein_score_ja.json'

        #既知の生薬-タンパク質 {生薬名:[protein1,protein2]}                                                                       
        infile2 = './data/positive_crude_protein_ja.json'

        #生薬-タンパク質スコアの読み込み                                                                                          
        with open(infile) as wf:
            crude_protein = json.load(wf)

        #既知の生薬-タンパク質の読み込み                                                                                          
        with open(infile2) as f:
            known_crude_protein = json.load(f)

    #生薬の詳細データ　[crude, kampo, compound, kegg_link, id]
    df_crude = pd.read_table('./data/all_crude_kampo_compound4.tsv',sep='\t')

    #化合物詳細データ　[2797化合物, 6レコードでInChI情報なし, link, crude, kampo, Name]                                         
    df_compound = pd.read_table('./data/all_compound_information.tsv',sep='\t')

    #漢方薬に含まれる生薬データ
    df_kampo_crude = pd.read_table('./data/example_194_kampo.tsv',sep='\t')

    #新規適応疾患(正規化のスコア0.5以上)
    input_disease = pd.read_table('./data/Differential_from_kp_matrix5_for_heatmap_normalized_and_scaled_converted_ja.txt',sep='\t')

    #漢方の詳細データ　[KNApSAcK ID, Multi-herbal medicine name, based_link, link, Multi-herbal medicine name japanese] 
    df_kampo = pd.read_table('./data/Kampo_KNApSAcK_ID2_ja.tsv',sep='\t')

#hsaとsymbol対応リスト+タンパク質詳細データ [hsa, symbol, description, NCBI_gene_ID, NCBI_URL]
df_protein = pd.read_table('./data/hsa_symbol_description_geneID_URL.txt',sep='\t')



# =========================================================================

# dictを用いたIDの変換
def convert_by_dict_list(id, conv_dict_list):
    try:
        x = conv_dict_list[id]
    except:
        x = []
    return x


# hsaかkoへのID変換(set型 -> set型)
def hsa2ko(hsa_set):
    return set(chain.from_iterable(list(map(lambda x: convert_by_dict_list(x, hsa2ko_dict),hsa_set))))


# FDR補正をする関数(pd.DataFrame -> pd.DataFrame)
def FDR(p):
    # ---------------------------------------------------------------------------------
    # FDR補正
    # 
    p = p.sort_values(by=1, ascending=False)
    b = p.iloc[0,1]  # 一番大きなp値(0番目のp値)を関数p_to_qの前の値として設定．
                     # bは一つ前の関数p_to_qの出力値

    M = len(p.iloc[:,1]) # P値を計算したパスウェイの数
    def p_to_q(index):             # 1番目から最後まで順に変換 range(1,M)
        nonlocal b                 # 一つ外側の関数の変数を用いる
        p_value = p.iloc[index, 1]
        b = min((M/(M-index))*p_value, b)
        return b

    q = list(map(p_to_q, range(1,M)))   # FDR補正により最も大きいp値以外のp値を補正
    q.insert(0, p.iloc[0,1])            # 最も大きいp値を追加(最も大きいものは補正しない)
    p['q'] = q
    # -----------------------------------------------------------
    # データの整形
    p = p.sort_values(by=0)                    # pathwayのID順にsort
    p = p[[0,'q']]                             # pathway IDと補正したp値(q値)のみにする．
                                               # すなわち，補正前のp値の削除
    p.columns = ['KEGG disease pathway ID','Disease score']  #列の名前を変更

    return p



# ====================================================================
# パスウェイエンリッチメント解析でFDR補正済みp値を計算する関数
# (set, set -> pd.DataFrame)
def enrichiment_analysis(Selected_genes, All_genes):
    # フィッシャーの正確確率検定を用いて値を求める関数
    # l: 対象とする全遺伝子の数
    # r: 化合物等の添加等で変換した遺伝子の数
    # k: パスウェイに含まれる遺伝子の数
    # z: 対象とする全遺伝子と化合物等の添加等で変化した遺伝子の集合の積集合の数
    def p_value(l, r, k, z):
        p = sum(map(lambda i: comb(k, i)*comb(l-k, r-i), range(z, min(r, k)+1)))/comb(l, r)
        return p
    
    # pathwayIDに対してp値を求める関数
    def p_value_pathway(pathway):
        l = len(All_genes)
        r = len(Selected_genes&All_genes)
        k = len(set(pathway2protein[pathway])&All_genes)
        z = len(Selected_genes&set(pathway2protein[pathway])&All_genes)
        return p_value(l, r, k, z)
    
    # 全てのパスウェイにおけるp値を計算
    p = pd.DataFrame(list(map(lambda x:(x, p_value_pathway(x)), pathway2protein.keys())))
    return FDR(p)

#=======================================================================
#パスウェイエンリッチメント解析への関数

def kampo2genes(x):
    return hsa2ko(set(input_data_group.get_group(x)['hsa']))

def enrichiment_analysis_kampo(x):
    return enrichiment_analysis(kampo2genes(x), All_genes)

def kampo2genes2(x):
    return set(input_data_group.get_group(x)['hsa'])

def enrichiment_analysis_kampo2(x):
    return enrichiment_analysis(kampo2genes2(x), All_genes)

#=========================================================================


# ====================main=============================
input_crude_dict = dict()

if args[3] == 'Kampo':
    selected_kampo = args[4]

    df_crude_id_selected = pd.DataFrame(columns = ['Name','Kampo','Compound','KEGG drug URL','KEGG drug ID'])
    df_kampo_crude = df_kampo_crude[df_kampo_crude['kampo'] == selected_kampo]
    crude_list = df_kampo_crude.iat[0,1].split(',')
    for i in crude_list:
        df = df_crude[df_crude['Name'] == i]
        df_crude_id_selected = pd.concat([df_crude_id_selected, df], axis=0)

if args[3] == 'Custom':
    #入力した生薬と比率を辞書型
    ratio_sum = 0
    for i in range(len(args)):
        if i > 3:
            if i % 2 == 0:
                ratio_sum += float(args[i+1])
        
    #入力した比率を正規化
    #入力した生薬の情報を抽出
    df_crude_id_selected = pd.DataFrame(columns = ['Name','Kampo','Compound','KEGG drug URL','KEGG drug ID'])
    for i in range(len(args)):
        if i > 3:
            if i % 2 == 0:
                input_crude_dict[args[i]] = float(args[i+1]) / ratio_sum
    
                df = df_crude[df_crude['Name'] == args[i]]
                df_crude_id_selected = pd.concat([df_crude_id_selected, df], axis=0)

#生薬情報を出力
crude_list = df_crude_id_selected.to_numpy().tolist()
df_crude_id_selected = df_crude_id_selected.reindex(columns=['Name','KEGG drug ID','KEGG drug URL'])
df_crude_id_selected['Name'] = df_crude_id_selected['Name'].str.replace('_', ' ')
df_crude_id_selected = df_crude_id_selected[~df_crude_id_selected.duplicated()]
df_crude_id_selected = df_crude_id_selected.sort_values('KEGG drug URL')
df_crude_id_selected.to_csv('./output/inputCrudeDrugs.csv',sep=',',index=False)


#入力した生薬を含む漢方薬の情報を抽出
df_crude_contain_kampo = pd.DataFrame(columns = ['KNApSAcK ID', 'Name','KNApSAcK URL', 'Crude drug'])
for crude_kampo in crude_list:
    contain_kampo = str(crude_kampo[1]).split(',')
    for kampo in contain_kampo:
        df = df_kampo[df_kampo['Name'] == kampo]
        df['Crude drug'] = crude_kampo[0]
        df_crude_contain_kampo = pd.concat([df_crude_contain_kampo, df], axis=0)
        

#生薬を含む漢方薬の情報を出力
df_crude_contain_kampo = df_crude_contain_kampo.reindex(columns=['Crude drug', 'Name', 'KNApSAcK ID', 'KNApSAcK URL'])
df_crude_contain_kampo['Crude drug'] = df_crude_contain_kampo['Crude drug'].str.replace('_', ' ')
df_crude_contain_kampo = df_crude_contain_kampo.fillna(0)
df_crude_contain_kampo = df_crude_contain_kampo[~df_crude_contain_kampo.duplicated()]
df_crude_contain_kampo.to_csv('./output/multiHerbalMedicines_in_theCrudeDrug.csv',sep=',',index=False)



#入力した生薬に含まれる化合物の情報を抽出
df_crude_constituent_compound = pd.DataFrame(columns = ['Crude drug', 'Name', 'KNApSAcK compounds ID', 'KNApSAcK compounds URL'])
df_constituent_compound = pd.DataFrame(columns = ['Name', 'KNApSAcK compounds ID', 'KNApSAcK compounds URL'])
for crude_compound in crude_list:
    constituent_compound = str(crude_compound[2]).split(',')
    for compound in constituent_compound:
        df = df_compound[df_compound['KNApSAcK compounds ID'] == compound]
        df_constituent_compound = pd.concat([df_constituent_compound, df], axis=0)
        df['Crude drug'] = crude_compound[0]
        df_crude_constituent_compound = pd.concat([df_crude_constituent_compound, df], axis=0)

#生薬に含まれる化合物の情報を出力
df_crude_constituent_compound['Crude drug'] = df_crude_constituent_compound['Crude drug'].str.replace('_', ' ')
df_crude_constituent_compound = df_crude_constituent_compound[~df_crude_constituent_compound.duplicated()]
df_crude_constituent_compound.to_csv('./output/constituentCompounds_in_theCrudeDrugs.csv',sep=',',index=False)

#入力された生薬に含まれる化合物全てを出力
df_constituent_compound = df_constituent_compound[~df_constituent_compound.duplicated()]
df_constituent_compound.to_csv('./output/constituentCompounds.csv',sep=',',index=False)



if args[3] == 'Kampo':
    #予測されたタンパク質に*をつけたリスト作成                                                                                    
    new_medicine = list()
    for medicine,protein_dict in kampo_protein.items():
        if selected_kampo == medicine:
            for protein, score in protein_dict.items():
                new_medicine.append([protein, score, '*'])

    df = pd.DataFrame(new_medicine,columns=['hsa','Interaction score','known'])
    df_s = df.sort_values('Interaction score', ascending=False)


if args[3] == 'Custom':
    #入力された生薬の生薬-タンパク質スコアを抽出
    #生薬-タンパク質スコアに比率をかける
    new_herb_dict = dict()
    for herb, ratio in input_crude_dict.items():
        p_s = new_herb_dict.setdefault(herb,dict())
        for crude, protein_dict in crude_protein.items():
            if herb == crude:
                for protein,score in protein_dict.items():
                    p_s[protein] = score * ratio

    #重複しているタンパク質スコアを足す
    new_multi_herb_protein = dict()
    for crude,protein_dict in new_herb_dict.items():
        c_s = new_multi_herb_protein.setdefault('multi_herbal_medicine',dict())
        for protein, score in protein_dict.items():
            if protein not in c_s:
                c_s[protein] = score
            else:
                c_s[protein] = c_s[protein] + score

    #予測されたタンパク質に*をつけたリスト作成                                                                                    
    new_medicine = list()
    for medicine,protein_dict in new_multi_herb_protein.items():
        for protein, score in protein_dict.items():
            new_medicine.append([protein, score, '*'])
                
    #new_medicineのリストをデータフレームに変換
    #スコアの正規化                                                                                                               
    #閾値0.77以上のタンパク質を抽出
    df = pd.DataFrame(new_medicine,columns=['hsa','Interaction score','known'])
    df_s = df.sort_values('Interaction score', ascending=False)
    max = df_s['Interaction score'].max()
    df_s['Interaction score'] = df_s['Interaction score'] / max

df_input = df_s[df_s['Interaction score'].astype('float') > 0.77]



#known_protein_list  スコア0.77以上のタンパク質に含まれない既知のタンパク質リスト
#known_protein_list2 スコア0.77以上のタンパク質に含まれる既知のタンパク質リスト
#known_protein_list3 全ての既知のタンパク質リスト
known_protein_list = list()
known_protein_list2 = list()
known_protein_list3 = list()
if args[3] == 'Kampo':
    for known_kampo, known_protein in known_kampo_protein.items():
        if selected_kampo == known_kampo:
            for protein in known_protein:
                if [protein, 1, ''] not in known_protein_list:
                    if protein not in df_input['hsa'].to_list():
                        known_protein_list.append([protein, 1, ''])
                    else:
                        known_protein_list2.append(protein)
                if protein not in known_protein_list3:
                    known_protein_list3.append(protein)

if args[3] == 'Custom':
    for input_crude, input_ratio in input_crude_dict.items():
        for known_crude, known_protein in known_crude_protein.items():
            if input_crude == known_crude:
                for protein in known_protein:
                    if [protein, 1, ''] not in known_protein_list:
                        if protein not in df_input['hsa'].to_list():
                            known_protein_list.append([protein, 1, ''])
                        else:
                            known_protein_list2.append(protein)
                    if protein not in known_protein_list3:
                        known_protein_list3.append(protein)

    
#スコア0.77以上の既知のタンパク質の＊を消す
for known_under_077 in known_protein_list2:
    for i  in range(len(df_input)):
        if known_under_077 == df_input.iat[i,0]:
            df_input.iat[i,2] = ''
                    
#スコア0.77以上のタンパク質に含まれない既知のタンパク質リストをデータフレームに変換
#スコア0.77以上のタンパク質と結合
df_known_protein = pd.DataFrame(known_protein_list, columns=['hsa','Interaction score','known'])
df_input = pd.concat([df_input, df_known_protein], axis=0)

for i in range(len(df_s)):
    if df_s.iat[i,0] in known_protein_list3:
        df_s.iat[i,1] = 1


#標的のタンパク質に詳細情報を結合
df_input2 = df_input
df_input = pd.merge(df_input, df_protein)

#未知のタンパク質symbolの末尾に*をつける
for i in range(len(df_input)):
    if df_input.iat[i,2] == '*':
        df_input.iat[i,3] = df_input.iat[i,3] + '*'

#標的タンパク質情報を出力
df_input = df_input.reindex(columns=['Name', 'Description','Interaction score', 'NCBI gene ID','NCBI URL'])               
df_input = df_input.sort_values('Name')
df_input = df_input[~df_input.duplicated()]
df_input.to_csv('./output/targetProtein.csv',sep=',',index=False)


#パスウェイ解析の標的とするタンパク質データ
if args[3] == 'Kampo':
    df_input2['medicine'] = selected_kampo
if args[3] == 'Custom':
    df_input2['medicine'] = 'multi_herbal_medicine'
    
input_data = df_input2
input_data_group = input_data.groupby('medicine')
kampo_list = list(input_data_group.groups.keys())


# すべての漢方で計算(マルチスレット計算)
#p-valueを-log10をとり、Disease scoreに変換
if args[1] == 'KEGG':
    kampo_result = dict(zip(kampo_list, list(map(enrichiment_analysis_kampo, kampo_list))))
if args[1] == 'UMLS':
    kampo_result = dict(zip(kampo_list, list(map(enrichiment_analysis_kampo2, kampo_list))))
    
for kampo,pathway_matrix in kampo_result.items():
    for pathway_num in range(len(pathway_matrix)):
        pathway_matrix.iat[pathway_num,1] = -1 * math.log10(pathway_matrix.iat[pathway_num,1])

        
# ファイルの書き出しする関数
def write_data(data, file_name):
    #Disease scoreが高いもの10を抽出
    pathway_include_ko = list()
    data2 = data.sort_values('Disease score',ascending=False)
    data2 = data2[data2['Disease score'] > 1.3]

    
    #pathwayに含まれるタンパク質情報を結合
    df_pathway_in_protein = pd.DataFrame(columns=['pathway','hsa'])
    for i in range(len(data2)):
        df = df_pathway_in_hsa[df_pathway_in_hsa['pathway'] == data2.iat[i,0]]
        df_pathway_in_protein = pd.concat([df_pathway_in_protein, df], axis=0)

        
    #標的タンパク質スコア情報と結合
    df = pd.merge(df_pathway_in_protein, df_s)
    df = pd.merge(df,df_protein)
    df = df.sort_values('pathway')

    #エンリッチされたタンパク質数を数える
    top10_list = data2['KEGG disease pathway ID'].to_list()
    top10_count_list = list()
    for top10 in top10_list:
        df_bool = df[df['pathway'] == top10]
        df_bool2 = (df_bool['Interaction score'].astype('float') >= 0.77)
        top10_count_list.append([top10,df_bool2.sum()])

    if args[1] == 'KEGG':
        #Disease scoreがtop10のパスウェイ情報を出力
        df_count = pd.DataFrame(top10_count_list, columns=['KEGG disease pathway ID','Count'])
        data3 = pd.merge(data2, df_pathway)
        data3 = pd.merge(data3, df_count)
        data3 = data3.reindex(columns=['Name', 'Count', 'Disease score', 'KEGG disease pathway ID', 'KEGG disease pathway map URL'])

    if args[1] == 'UMLS':
        data2 = data2.rename(columns={'KEGG disease pathway ID':'UMLS ID'})
        df_count = pd.DataFrame(top10_count_list, columns=['UMLS ID','Count'])
        data3 = pd.merge(data2, df_pathway)
        data3 = pd.merge(data3, df_count)
        data3 = data3.reindex(columns=['Name', 'Count', 'Disease score', 'UMLS ID', 'UMLS ID URL'])

        
    if args[3] == 'Kampo':
        input_disease2 = input_disease.reindex(columns=['Pathway',selected_kampo])
        input_disease3 = input_disease2[input_disease[selected_kampo].astype('float') > 0.5]
        new_disease_list = list(input_disease3['Pathway'])

        for i in range(len(data3)):
            if data3.iat[i,3] in new_disease_list:
                data3.iat[i,0] = data3.iat[i,0] + '*'

    data3 = data3[~data3.duplicated()]
    data3.to_csv('./output/applicableDisease.csv',sep=',',index=False)

    df['color'] = ''
    for i in range(len(df)):
        if float(df.iat[i,2]) >= 0.77:
            if str(df.iat[i,1]) not in known_protein_list3:
                df.iat[i,4] = df.iat[i,4] + '*'
        else:
            df.iat[i,8] = 'gray'    

            
    #pathwayに含まれるタンパク質情報を出力
    if args[1] == 'KEGG':
        df = df.reindex(columns=['pathway', 'Name', 'Description', 'Interaction score', 'NCBI gene ID', 'NCBI URL', 'color'])
    if args[1] == 'UMLS':
        df = df.rename(columns={'pathway':'UMLS ID'})
        df = df.reindex(columns=['UMLS ID', 'Name', 'Description', 'Interaction score', 'NCBI gene ID', 'NCBI URL', 'color'])
    df = df[~df.duplicated()]
    df.to_csv('./output/protein_list_in_count_applicableDisease.csv',sep=',',index=False)


# ファイルの書き出しの実行
list(map(lambda x: write_data(kampo_result[x], x), kampo_list))
