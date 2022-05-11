import pandas as pd
import argparse
import functools
import re
parser = argparse.ArgumentParser(description='make merged dataframe in the format of pos, unmet in p, total in p, unmet in c, total in c. The script can also merge multiple probands or controls.')
parser.add_argument('--i', action = 'store', dest = 'input', required = True, type=str,nargs='+', help='input file list for proband and sibling')
parser.add_argument('-o',action='store',dest='output',required=True,type=str,help='output filename')
parser.add_argument('--col',action='store',dest='columns',required=True,type=str,nargs='+',help='prefix for column names i.e. [12220_p1, 12220_s1] in which case the column names will be [12220_p1_met, 12220_p1_total,\
                   12220_s1_met,12220_s1_total]')
args  = parser.parse_args()
infile=list(args.input)
outfile=str(args.output)
col=list(args.columns)
rp=re.compile(r'.*p[0-9].*')
#rp=re.compile(r'.*p[0-9].*')
rc=re.compile(r'.*s[0-9].*')
pinfile=[i for i in infile if rp.match(i) ]
cinfile=[i for i in infile if rc.match(i) ]
#prefixlist=['11918/11918_p1_','12220/12220_p1_','11918/11918_s1_','12220/12220_s1_']  
columns=[]
for _ in col:
    for __ in ['_met','_total']:
        columns.append(_+__)
pmet_list=[]
cmet_list=[]
for i, _pinfile in enumerate(pinfile):
    pmet_list.append(pd.read_table(_pinfile,header=None))
for j,_cinfile in enumerate(cinfile):
    cmet_list.append(pd.read_table(_cinfile,header=None))
df=functools.reduce(lambda x, y: pd.merge(x[x.iloc[:,2]>10],y[y.iloc[:,2]>10],on=[0]).drop_duplicates(),[i for i in pmet_list]+[i for i in cmet_list])
df.columns=['pos']+columns
df.to_csv(outfile,index=None,sep='\t')

