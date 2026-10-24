import glob
import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
root=input("enter the path to root directory").strip().strip('"')
cutoff=float(input("enter the minimum coverage"))
csv=glob.glob("**/dates.csv",root_dir=root,recursive=True,include_hidden=True)
tsv=glob.glob("**/nc.tsv",root_dir=root,recursive=True,include_hidden=True)
path1=os.path.join(root,tsv[0])
path2=os.path.join(root,csv[0])
t1=pd.read_csv(path1,sep="\t",usecols=["coverage","index","seqName","clade","aaDeletions","aaInsertions","aaSubstitutions","proposedSubclade","subclade","short-clade","qc.overallStatus","legacy-clade","glycosylation","qc.overallScore","privateAaMutations.unlabeledSubstitutions"])
t2=pd.read_csv(path2,sep=",",usecols=["acc_no","dates"])
print(t1)
print(t2)
t2["dates"]=pd.to_datetime(t2["dates"],errors="coerce")
df=pd.merge(t1,t2,left_on="seqName",right_on="acc_no",how="left")
df=df.dropna(subset=["dates"])
df=df[(df["qc.overallStatus"].isin(["good","mediocre"]))&(df["coverage"]>=cutoff)]
target=input("enter the names of target  protein,seperated by space").split()
L=[]
for name,date,subs,pvtmut in zip(df["seqName"],df["dates"],df["aaSubstitutions"],df["privateAaMutations.unlabeledSubstitutions"]):
    if pd.isna(subs):
        continue
    value=subs.split(",")
    for mut in value:
        gene=mut.split(":")[0]
        if gene in target[:]:
            L.append({"seqName":name,"dates":date,"mutation":mut,"privateAaMutations":pvtmut})
            print(mut)        
df2=pd.DataFrame(L)            
if df2.empty:
    print("enter the correct name of spike protein")
else:
    df2=df2.dropna(subset=["dates"])
    print(df2)
df2=df2.drop_duplicates(["seqName","mutation"])    
df2["period"]=df2["dates"].dt.to_period("M")
df["period"]=df["dates"].dt.to_period("M")
counts=df2.groupby(["period","mutation"]).size().unstack(fill_value=0)
print(counts)
denum=df.groupby(["period"])["seqName"].nunique()
freq=counts.div(denum,axis=0)
print(freq)
freq2=freq.reset_index().melt(id_vars="period",value_name="frequency",var_name="mutation")
freq2["period"]=freq2["period"].astype(str)
lineplot=sns.lineplot(data=freq2,x="period",y="frequency",hue="mutation",marker="o")
plt.legend(loc="upper left",title="mutation")
plt.savefig("lineplot.png")
plt.show()
heatmap=sns.heatmap(freq.T,annot=True,vmin=0,vmax=1,cmap="coolwarm")
plt.savefig("heatmap.png")
plt.show()