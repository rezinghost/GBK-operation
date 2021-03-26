import os
from Bio import SeqIO
import shutil
import glob
import xlrd

#所有创建的新目录不要和原数据在一个路径下
#全部铁载体GBK → 筛选特定种属
#（否定）读取每个GBK文件逐一选出

objName = 'Enterobacter' #特定的种属名

drlist = [] #空表放符合条件文件名
a = 0 #计数

path = "D:/1/test data" #存放原数据的位置
files_list = os.listdir(path)
targetPath = r'D:/1/siderophore'  # 目标文件夹
if not os.path.exists(targetPath):
    os.makedirs(targetPath)  # 创建目标文件夹

for drname in files_list:
    filepath = 'D:/1/test data/' + drname
    listnames = glob.glob(filepath + "/*.xls")
    listname = ''.join(listnames)
    print('==============loading==============' + '\n' + listname)
    try:
        workbk = xlrd.open_workbook(listname)
        workst = workbk.sheet_by_index(0)
        for row in range(workst.nrows):
            for col in range(workst.ncols):
                cell_value = str(workst.cell_value(row, col)).strip()
                if cell_value == "siderophore":
                    fileindex = row
                    Fname = workst.cell_value(fileindex, 0)
                    sourcePath = r'D:/1/test data/' + drname

                    objFileName = Fname + '.' + 'cluster' + str(fileindex).zfill(3) + '.gbk'
                    print(objFileName)
                    for i in os.listdir(sourcePath):
                        if i == objFileName:
                            a += 1
                            str_list = [drname,str(a)] #“文件名 铁载体num”
                            sta = ' '.join(str_list)
                            drlist.append(sta)
                            shutil.copyfile(sourcePath + '/' + i, targetPath + '/' + i)

        a = 0 #重置计数
                    #print(fileindex, Fname)
    except:
        print(drname + '\n' + '■■■■■■■■FAIL■■■■■■■■')
        continue

#文件名列表写入
f = open("list.txt", "w+")
for line in drlist:
    f.write(line + '\n')
f.close()


#读取已经筛选的gbk文件
sdppath = r'D:/1/siderophore/' #已筛选为铁载体的GBK文件路径
objpath = r'D:/1/' + objName
if not os.path.exists(objpath):
    os.makedirs(objpath)  #创建目标文件夹

files_list2 = os.listdir(sdppath)
b = 0 #计数
for filename in files_list2:
    filepath = sdppath + filename #获得绝对路径
    rec = SeqIO.read(filepath, "genbank")
    taxonomy = rec.annotations["taxonomy"]
    genus = taxonomy[-2] #种:-1 属:-2 ...以此类推
    if genus == objName:
        shutil.copyfile(filepath, objpath + '/' + filename)
        print(filename + "拷贝成功")
        b += 1
print('共得到 [%s] 的GBK文件 %d 个.' %(objName, b))
