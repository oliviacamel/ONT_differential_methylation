import numpy as np
from scipy import stats
from multiprocessing import Pool
import pandas as pd
import argparse
parser = argparse.ArgumentParser(description='mcmc approach to call significantly differential cpg dinucleotides')
parser.add_argument('-i', action = 'store', dest = 'inputfile', required = True, type=str, help='file containing df, grouped proband unmet,n and grouped control unmet,n')
parser.add_argument('-o',action='store',dest='output',required=True,type=str,help='output containing pos, unmet_ratio of proband over control, mean_theta1, mean_theta2, pval')
parser.add_argument('-p', action = 'store', dest = 'n_prob', required = True, type=str, help='number of probands')
parser.add_argument('-c',action='store',dest='n_ctrl',required=True,type=str,help='number of control')
args = parser.parse_args()
infile=str(args.inputfile)
outfile=str(args.output)
n_p=int(args.n_prob)
n_c=int(args.n_ctrl)

def dbeta(theta, a, b):
	return stats.beta.logpdf(theta, a, b)
def dbinom(x, n,p):
	return stats.binom.logpmf(x,n,p)


def like(x1,n1 ,x2, n2,para):
    [theta1,theta2]=para
    return np.sum([dbinom(x1, n1, theta1),dbinom(x2, n2, theta2)])

def prior(para):

    [theta1, theta2] = para
    return np.sum([dbeta(theta1, 0.5, 0.5), dbeta(theta2, 0.5, 0.5)])

def posterior(x1, n1,x2, n2,para):
    [theta1, theta2]=para
    return like(x1,n1,x2,n2, [theta1, theta2]) + prior(para)

def computeHDI(chain, interval = .95):
	'''
	Compute 95% highest density interval (HDI)
	'''
	# sort chain using the first axis which is the chain
	chain.sort()
	# how many samples did you generate?
	nSample = chain.size    
	# how many samples must go in the HDI?
	nSampleCred = int(np.ceil(nSample * interval))
	# number of intervals to be compared
	nCI = nSample - nSampleCred
	# width of every proposed interval
	width = np.array([chain[i+nSampleCred] - chain[i] for  i in range(nCI)])
	# index of lower bound of shortest interval (which is the HDI) 
	best  = width.argmin()
	# put it in a dictionary
	#HDI   = {'Lower': chain[best], 'Upper': chain[best + nSampleCred], 'Width': width.min()}
	HDI = [chain[best], chain[best + nSampleCred]]
	return HDI
    

def beta_bayes( id, s1, s2, niter = 18000, nburn_in = 1000):

    #	s to be 2d arrays of [[k1,n1],[k2,n2]...]
    theta1 = []	#means sampled by MCMC for s1
    theta2 = []	#means sampled by MCMC for s2
     	   # run MCMC (Metropolis-Hastings)
    
    #smoothing in the case of p=0
    x1,n1 ,x2, n2 = np.sum(s1[:,0]),np.sum(s1[:,1]), np.sum(s2[:,0]),np.sum(s2[:,1])
    if x1==0:
        x1=1
        n1+=1
    if x2==0:
        x2=1
        n2+=1
    if x1==n1:
        n1+=1
    if x2==n2:
        n2+=1
    _theta1,_theta2=x1/n1,x2/n2
    parameters = np.array([_theta1,_theta2])	
    #generate jitter from normal distribution with mu=0, std=0.1
    _drawlist=[np.random.normal(0,0.1,niter*20),np.random.normal(0,0.1,niter*20)]
    drawlist=[(parameters[0] + _drawlist[0]) [((parameters[0] + _drawlist[0]) > 0)&((parameters[0] + _drawlist[0])<1)], (parameters[1] + _drawlist[1]) [((parameters[1] + _drawlist[1]) > 0)&((parameters[1] + _drawlist[1]) < 1)]]
    for iteration in np.arange(1,niter):
        candidates=[drawlist[0][iteration], drawlist[1][iteration]]
        ratio = np.exp(posterior(x1,n1 ,x2, n2, candidates) - posterior(x1,n1 ,x2, n2, parameters))
        if np.random.uniform() < ratio:
            parameters = candidates
        if iteration < nburn_in:
            continue
        theta1.append(parameters[0])
        theta2.append(parameters[1])
    	
    # calculate estimated means			
    theta1 = np.array(theta1)
    theta2 = np.array(theta2)
    est_theta1 = theta1.mean()	#estimated mu1
    est_theta2 = theta2.mean()	#estimated mu2
    # calculate probability
    diff = (theta1 - theta2)
    diff_median = np.median(diff)
    if diff_median < 0:
        prob = np.mean(diff > 0)
    elif diff_median > 0:
        prob = np.mean(diff < 0)     	
     # calculate HDI
    #diff_HDI_h, diff_HDI_l = computeHDI(diff)      	
    # CpG_ID, mean of group1, mean of group2, diff of mean, 95%HDI_low, HDI_high, probability
    results= [ id,  _theta1/_theta2,est_theta1, est_theta2, prob]
    
    #results=[id,0,0,0,0]
    return(results)
if __name__=='__main__':
    df=pd.read_table(infile)
    def dffilter(df,n_p,n_c):
        #drop duplicates
        df=df.loc[ df.iloc[:,1:].drop_duplicates().index].reset_index(drop=True)
        #do differential calling on CpGs with only 2 fold differences.
        indxlist=[]
        for i in range(n_p*2+1,df.shape[1],2):
            for j in range(1,n_p*2+1,2): 
                indx=np.where((df.iloc[:,j]/df.iloc[:,j+1]/(df.iloc[:,i]/df.iloc[:,i+1])<2)&(df.iloc[:,j]/df.iloc[:,j+1]/(df.iloc[:,i]/df.iloc[:,i+1])>0.5))[0]
                indxlist.append(indx)
        if len(indxlist)==1:
            target=indxlist[0]
        else:    
            for i,indx in enumerate(indxlist):
                if i==0:
                    target=np.intersect1d(indxlist[i],indxlist[i+1])
                else:
                    target=np.intersect1d(target,indxlist[i])
        df=df.loc[~df.index.isin(target)].reset_index(drop=True)
        return(df)
    filtereddf=dffilter(df=df,n_p=n_p,n_c=n_c)
    n=0
    for k in range(0,filtereddf.index.values[-1]-int(filtereddf.shape[0]/250), int(filtereddf.shape[0]/250))  :
        results_list=[]
        with Pool() as pool:
            multiple_results = [pool.apply_async(beta_bayes, (filtereddf.values[i,0],filtereddf.iloc[i,1:2*n_p+1].values.reshape([n_p,2]) , filtereddf.iloc[i,2*n_p+1:].values.reshape([n_c,2])))  for i in range(k,k+int(filtereddf.shape[0]/250))]
            for res in multiple_results:
                res.wait()
                results_list.append(res.get())
        pool.close()
        pool.join()
        if n==0:
            results_list=pd.DataFrame(results_list)
            results_list.columns=['pos','ratio','theta1','theta2','p_raw']
            results_list.to_csv(outfile,mode='a',index=None,sep='\t')
        else:
            pd.DataFrame(results_list).to_csv(outfile,mode='a',index=None,sep='\t',header=None)
        n+=1
    results_list=[]
    with Pool() as pool:
        multiple_results = [pool.apply_async(beta_bayes, (filtereddf.values[i,0],filtereddf.iloc[i,1:2*n_p+1].values.reshape([n_p,2]) , filtereddf.iloc[i,2*n_p+1:].values.reshape([n_c,2])))  for i in range(k+int(filtereddf.shape[0]/250),filtereddf.index.values[-1] + 1 )]
        for res in multiple_results:
            res.wait()
            results_list.append(res.get())
    pool.close()
    pool.join()
    pd.DataFrame(results_list).to_csv(outfile,mode='a',index=None,sep='\t',header=None)
    
    
    
