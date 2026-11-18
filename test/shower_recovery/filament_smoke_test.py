"""Small display-inspired toy study; run the actual C++ splitter, no tuning.
Artifact hits have truth label -1 and are excluded from ownership metrics.
"""
import argparse
import json
import subprocess
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CASES = [('contact','Small outward branch contact'),('streak','Fixed-wire connecting streak'),
         ('unequal','Unequal showers + gaps'),('offset','Branch contact; vertex shifted 2 cm'),
         ('single','Single shower with forward branches')]


def generate(kind, seed):
    rng=np.random.default_rng(seed)
    rows=[]
    def add(x,y,label,scale=1.):
        for w,t in zip(x,y):
            rows.append([round(w/.3)*.3,t+rng.normal(0,.045),rng.lognormal(0,.4)*scale,label])
    def shower(label, end):
        slope=0.16 if label else -0.16
        x=np.arange(5.2,end,.24)
        y=slope*x+.12*np.sin(x/2.8)
        add(x,y,label)
        for start,fan in [(12.,.07),(18.,-.06),(23.,.12),(28.,-.09)]:
            if start>=end-3: continue
            x=np.arange(start,end-rng.uniform(0.,3.),.26)
            y=slope*x+.12*np.sin(x/2.8)+fan*(x-start)+.22*np.sin((x-start)/2.)
            # Narrow strands with occasional missing hits.
            keep=rng.random(len(x))>.07
            add(x[keep],y[keep],label,.7)
    shower(0,42.)
    if kind!='single': shower(1,28. if kind=='unequal' else 40.)
    if kind in ('contact','offset','unequal'):
        x=np.arange(7.,17.,.22)
        # A forward-going offshoot bends toward the lower shower locally.
        y=.16*x+.12*np.sin(x/2.8)-3.45*np.exp(-((x-11.)/2.4)**2)
        add(x,y,1,.7)
    if kind=='streak':
        y=np.arange(-3.8,4.2,.22)
        add(np.full(len(y),15.),y,-1,.45)
    data=np.array(rows)
    if kind=='unequal':
        data=data[~((data[:,0]>23.7)&(data[:,0]<24.3))]
    return data,np.array([0.,2.]) if kind=='offset' else np.zeros(2)


def run(driver,data,vertex):
    inp=data[:,:3].copy();inp[:,:2]-=vertex
    text=str(len(inp))+'\n'+'\n'.join(' '.join(map(str,r)) for r in inp)+'\n'
    out=subprocess.run([str(driver)],input=text,text=True,capture_output=True,check=True).stdout.splitlines()
    accepted=bool(int(out[0].split()[0]));labels=np.array(list(map(int,out[1:])))
    assert len(labels)==len(data)
    result={'accepted':accepted,'charge_accuracy':None,'min_shower_recall':None}
    if accepted:
        assert set(labels)=={0,1}
        valid=data[:,3]>=0
        truth=data[:,3].astype(int)
        if len(set(truth[valid]))==2:
            score=lambda z:data[valid&(z==truth),2].sum()
            if score(1-labels)>score(labels):labels=1-labels
            result['charge_accuracy']=float(score(labels)/data[valid,2].sum())
            result['min_shower_recall']=float(min((labels[truth==i]==i).mean() for i in (0,1)))
    return result,labels


def main():
    p=argparse.ArgumentParser();p.add_argument('--driver',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    summaries=[];examples=[];trials=[]
    for kind,title in CASES:
        results=[]
        for seed in (310,311,312):
            data,vertex=generate(kind,seed)
            result,labels=run(args.driver,data,vertex);results.append(result)
            trials.append(dict(case=kind,seed=seed,**result))
            if seed==310:
                examples.append((kind,title,data,vertex,result,labels))
                np.savetxt(args.output/f'{kind}.csv',np.column_stack((data,labels)),delimiter=',',
                  header='wire_cm,drift_cm,charge,truth,reco',comments='')
        accepted=[r for r in results if r['accepted']]
        scores=[r['charge_accuracy'] for r in accepted if r['charge_accuracy'] is not None]
        summaries.append({'case':kind,'accepted':len(accepted),'trials':3,
          'mean_charge_accuracy':float(np.mean(scores)) if scores else None,
          'min_shower_recall':min((r['min_shower_recall'] for r in accepted if r['min_shower_recall'] is not None),default=None)})
    fig,axes=plt.subplots(2,3,figsize=(14,7.5),sharex=True,sharey=True)
    palette=np.array(['#1677b8','#e88129'])
    for ax,(kind,title,data,vertex,result,labels) in zip(axes.flat,examples):
        true=data[:,3].astype(int);physical=true>=0
        colors=np.full(len(data),'#777777',dtype=object)
        if result['accepted']:colors[physical]=palette[labels[physical]]
        ax.scatter(data[:,0],data[:,1],c=colors,s=np.clip(9*data[:,2],3,20),linewidths=0)
        wrong=physical&(labels!=true) if result['accepted'] and kind!='single' else np.zeros(len(data),bool)
        ax.scatter(data[wrong,0],data[wrong,1],s=23,facecolors='none',edgecolors='#ce2363',linewidths=.8)
        ax.scatter(0,0,c='black',marker='*',s=95,zorder=4)
        if np.any(vertex):ax.scatter(*vertex,c='#ce2363',marker='x',s=60,zorder=5)
        status='No split' if not result['accepted'] else ('FALSE SPLIT' if kind=='single' else f"Correct charge: {result['charge_accuracy']:.1%}")
        ax.set_title(title+'\n'+status,fontsize=10)
        ax.set_xlim(-2,44);ax.set_ylim(-12,13);ax.set_aspect('equal');ax.grid(alpha=.15)
        ax.set_xlabel('Wire coordinate [cm]');ax.set_ylabel('Drift coordinate [cm]')
    axes.flat[-1].axis('off')
    axes.flat[-1].text(.05,.85,'Actual C++ splitter; defaults unchanged\n\nBlue/orange: output partitions\nPink rings: incorrect ownership\nGray: artifact or unsplit hits\nStar: true vertex; x: shifted vertex\n\nFirst fixed seed shown for each case.\n15 toy tests total; no detector simulation.\nNo cross-plane or 3D reconstruction.',transform=axes.flat[-1].transAxes,va='top',fontsize=11)
    fig.suptitle('Display-inspired filament tests',fontsize=15)
    fig.tight_layout();fig.savefig(args.output/'filament_tests.png',dpi=150);plt.close(fig)
    (args.output/'summary.json').write_text(json.dumps(summaries,indent=2)+'\n')
    (args.output/'trials.json').write_text(json.dumps(trials,indent=2)+'\n')
    print(json.dumps(summaries,indent=2))

if __name__=='__main__':main()
