import pandas as pd
import re
import argparse
from multiprocessing import Pool
import numpy as np
parser = argparse.ArgumentParser(description='creating dataframes where rows indicate CG dinucleotide pos and columns indicate numbers of met and unmet reads at indicated positions')
parser.add_argument('-i', action = 'store', dest = 'input_file', required = True, type=str, help='a sorted tsv file output by nanopolish')
parser.add_argument('-o',action='store',dest='output_file',required=True,type=str,help='destination file for dataframe')
args = parser.parse_args()
input_file=args.input_file
out=args.output_file
#hard threshold of -2 for unmet and +2 for met
def parse_cg(row):
    cgposlist=[]
    for i,cg in enumerate(re.finditer('CG',row.iloc[10])):
        if i>0:
            cgposlist.append(int(row.iloc[2])+cg.start()-first_cg_pos)
        else:
            first_cg_pos=cg.start()
    return(cgposlist)
def met_array_construct(chunk):
    ###############
    #pos met unmet#
    ###############
    array={}
    for i in range(len(chunk)):
        pos,ratio,num_motifs=chunk.iloc[i,2],chunk.iloc[i,5],chunk.iloc[i,9]
        if num_motifs==1:
            if pos in array.keys():
                if ratio<-2:
                    if 'unmet' in array[pos].keys():
                        array[pos]['unmet']+=1
                    else:
                        array[pos]['unmet']=1
                elif ratio > 2:
                    if 'met'  in array[pos].keys():
                        array[pos]['met']+=1
                    else:
                        array[pos]['met']=1
            else:
                if ratio <-2:
                    array[pos]={}
                    array[pos]['unmet']=1
                elif ratio >2:
                    array[pos]={}
                    array[pos]['met']=1
        else:
            pos_list=parse_cg(chunk.iloc[i])
            for _pos_ in pos_list:
                if _pos_ in array.keys():
                    if ratio<-2*num_motifs:
                        if 'unmet' in array[_pos_].keys():
                            array[_pos_]['unmet']+=1
                        else:
                            array[_pos_]['unmet']=1
                    elif ratio > 2*num_motifs:
                        if 'met'  in array[_pos_].keys():
                            array[_pos_]['met']+=1
                        else:
                            array[_pos_]['met']=1
                else:
                    if ratio <-2*num_motifs:
                        array[_pos_]={}
                        array[_pos_]['unmet']=1
                    elif ratio >2*num_motifs:
                        array[_pos_]={}
                        array[_pos_]['met']=1
    return(pd.DataFrame(array))


if __name__ == '__main__':
    pool=Pool(10)
    dflist=pool.map(met_array_construct,[chunk for chunk in  pd.read_table(input_file,chunksize=100000,header=None)])
    _=pd.DataFrame(columns=['met','unmet'])
    for i,df in enumerate(dflist):
        df=df.replace(np.nan,0)
        newdf=pd.concat((_,df.T.dropna()))
        newdfc=newdf.copy()
        newdfc.columns=['met','total']
        newdfc.iloc[:,0]=newdf.iloc[:,0].astype('int')
        if i==0:
            merge=newdfc
        else:
            merge=pd.concat((merge,newdfc))
    merge=merge.groupby(merge.index).aggregate({'met':np.sum,'total':np.sum})
    #smooth out edges, in case of 0 met, add 1 to met and total
    indx=np.where(merge.iloc[:,0]==0)[0]
    merge.iloc[indx,0]=merge.iloc[indx,0]+1
    merge.iloc[:,1]=(merge.iloc[:,0]+merge.iloc[:,1]).astype('int')
    merge.to_csv(out,sep='\t',header=None)

