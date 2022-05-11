import glob
import re
indx=['s1','p1']
chr=[i for i in range(1,23)]
familyid=['14455','12220','11918','12456']
rule all:
    input:
        expand("{familyid}/{familyid}_chr{chr}.nanopolish.CG.dinuc.differential.significance.tsv",familyid=familyid,chr=chr),  
        expand("{familyid}/{familyid}_{indx}_chr{chr}.nanopolish.CG.dinuc.met.perc.tsv",familyid=familyid,indx=indx,chr=chr)

checkpoint get_dist:
    input:
        rawtsv='{familyid}/{familyid}_{indx}_chr{chr}_sorted.tsv'
    output:
        outputdist='{familyid}/{familyid}_{indx}_chr{chr}.nanopolish.CG.dinuc.met.perc.tsv'
    resources:
        mem = 3,
        hrs = 120,
        threads = 10
    threads: 10
    shell:
        "module load miniconda/4.9.2;"
        "python get_percentage_met_at_CG_dinucleotide_position.py  -i {input.rawtsv} -o  {output.outputdist}"
        
def get_prob_files(wildcards):
    with checkpoints.get_dist.get(**wildcards).output[0].open() as f:
        l=glob.glob(str(wildcards.familyid)+'/'+str(wildcards.familyid)+'_'+str(wildcards.indx)+'_chr'+str(wildcards.chr)+'.nanopolish.CG.dinuc.met.perc.tsv')
        r=re.compile(r'.*p[0-9].*')
        newl=[i for i in l if r.match(i)]
        return(newl)

def get_sib_files(wildcards):
    with checkpoints.get_dist.get(**wildcards).output[0].open() as f:
        l=glob.glob(str(wildcards.familyid)+'/'+str(wildcards.familyid)+'_'+str(wildcards.indx)+'_chr'+str(wildcards.chr)+'.nanopolish.CG.dinuc.met.perc.tsv')
        r=re.compile(r'.*s[0-9].*')
        newl=[i for i in l if r.match(i)]
        return(newl)
rule merge:
    input: expand("{{familyid}}/{{familyid}}_{indx}_chr{{chr}}.nanopolish.CG.dinuc.met.perc.tsv" , indx=indx)
    output: "{familyid}/{familyid}_chr{chr}.combined.nanopolish.CG.dinuc.met.perc.tsv",
    resources:
        mem = 3,
        hrs = 120,
        threads = 1
    params:
        pcol="{familyid}_p1",
        scol="{familyid}_s2",
    threads: 1
    shell:
        "module load miniconda/4.9.2;"
        "python merge_proband_sib_CpG.py  --i {input}  -o {output} --col {params.pcol} {params.scol}" 
        
rule MCMC_bayes:
    input: "{familyid}/{familyid}_chr{chr}.combined.nanopolish.CG.dinuc.met.perc.tsv"
    output: "{familyid}/{familyid}_chr{chr}.nanopolish.CG.dinuc.differential.significance.tsv"
    resources:
        mem=3,    
        hrs=400,
        threads = 15
    threads : 15
    shell:
        "module load miniconda/4.9.2;"
        "python MCMC_bayes_for_differential_CpG_met.py  -i {input} -o {output} -p 1 -c 1" 
        

