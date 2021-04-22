from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
import os
import re


def getFileNames(path):
    '''
    gets the file names for all gbk files in the specific folder
    '''
    files = []
    for file in os.listdir(os.getcwd()+path):
        if file.endswith(".gbk"):
            files.append(os.getcwd()+path+file)
    return files

def findDomains(filename):
    for rec in SeqIO.parse(filename, 'genbank'):
        results = []
        seq = rec.seq
        id = rec.id
        for feature in rec.features:
            if feature.type == "aSDomain":
                f = feature.qualifiers

                domain_type = f['domain'][0]
                domain_seq = feature.extract(seq)
                domain_len = len(feature)
                region =  (feature.location.start , feature.location.end)# get a landscape of domain location


                Cdomain = 'Condensation'
                Adomain = 'AMP-binding'
                Tdomain = 'PCP'


                if  domain_type == Adomain:
                    substrate_sp = f['specificity']
                    results.append((id, domain_type, domain_len, domain_seq, substrate_sp))
                else:
                    results.append((id, domain_type, domain_len, domain_seq))

        return results


def main(folder, path):

    #make results dir
    if not os.path.exists(folder):
        os.makedirs(folder)

    for fileName in getFileNames(path):
        domains = findDomains(fileName)

        fpath,fname=os.path.split(fileName)

        savefilename = folder + fname.split('.gbk')[0] + '_domains.fasta'
        saveFiledm = open(savefilename, 'w')
        for index in range(len(domains)):
            saveFiledm.write('> %s_domain%d_%s_length:%d\n' % (domains[index][0], index + 1, domains[index][1],domains[index][2])) #id, domain_num, type, length
            if len(domains[index]) == 5:
                saveFiledm.write(str(domains[index][4]) + '\n') #specificity
            saveFiledm.write(str(domains[index][3]) + '\n') #sequence
            saveFiledm.write('\n')

        saveFiledm.close()


main('./domains/', '\\test\\')



