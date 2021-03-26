import os
from Bio import SeqIO
from collections import Counter
import pandas as pd
import matplotlib.pyplot as plt


'''
暂时储存说明：
s = 物种列表
r = 残基数量
p = 混杂的各个基因产物列表
'''

s = []
r = []
p = []


#读取已经筛选的gbk文件
path = r'D:/1/sdp'
files_list = os.listdir(path)
f_num = 0 #NRPS-based计数初始化
for filename in files_list:
    filepath = 'D:/1/sdp/' + filename #获得绝对路径
    rec = SeqIO.read(filepath, "genbank")


    #统计物种
    taxonomy = rec.annotations["taxonomy"]
    species = taxonomy[-1]
    s.append(species)


    #统计NRPS based
    features = rec.features

    for feature in features:
        f = feature.qualifiers
        try:
            us = str(f['sec_met']).upper() #字符串+大写化
            if 'NRPS' in us:
                f_num += 1
                pass
            else:
                continue
        except:
            continue


    #统计残基
    for i in range(len(features)):
        try:
            a = features[i]
            b = a.qualifiers
            c = b["translation"]
            r = r + c
        except:
            continue




'''
spresult = Counter(s)
df1 = pd.DataFrame.from_dict(spresult, orient='index').reset_index()
df1a = df1.rename(columns={'index':'species', 0:'num'})
df1a.to_excel('spnum.xlsx',index=None)
'''
allAA = "".join(r)
ctresult = Counter(allAA) #储存计数字典
AA = list(ctresult.keys()) #从dictkeys转化为列表储存
num = list(ctresult.values())
df2 = pd.DataFrame({'AA': AA,'num': num})
df2.to_excel('rsdnum.xlsx',index=None)

print('NRPS-based数量为：' + str(f_num) + '个\n')




#product

#可视化
