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

CASES = [('outward','Branching showers; no added crossing branch'),
         ('toward','Branch grows toward other shower'),
         ('unequal','Unequal branching showers + contact'),
         ('offset','Branching contact; vertex shifted 2 cm'),
         ('single','Single shower with several branches')]


def generate(kind, seed):
    rng=np.random.default_rng(seed)
    rows=[]
    axes=[-0.20,0.20]
    def strand(start, angle, length, label, charge, jitter=.07):
        # Persistent, slowly changing direction; all segments progress outward.
        # No deliberately returning/bending-back bridge.
        pos=np.array(start,dtype=float)
        points=[pos.copy()]
        heading=angle
        for step in range(int(length/.25)):
            heading += rng.normal(0.,.022)
            heading = np.clip(heading,-.85,.85)
            pos = pos + .25*np.array([np.cos(heading),np.sin(heading)])
            points.append(pos.copy())
        points=np.array(points)
        for i,(w,t) in enumerate(points):
            # Sparse local gaps and small transverse variation.
            if i>3 and rng.random()<.055: continue
            rows.append([round(w/.3)*.3,t+rng.normal(0.,jitter),
                         charge*rng.lognormal(0.,.5),label])
        return points
    nshowers=1 if kind=='single' else 2
    trunks=[]
    for label in range(nshowers):
        length=31. if kind=='unequal' and label==1 else 44.
        trunk=strand((0.,0.),axes[label],length,label,1.4,.04)
        trunks.append(trunk)
        # Both showers have multiple daughter strands and secondary branching.
        # Start positions, lengths, and directions change across seeds.
        for j,dist in enumerate((7.,12.,17.,22.,28.,34.)):
            if dist>length-4: continue
            start=trunk[min(int(dist/.25),len(trunk)-1)]
            fan=rng.uniform(.07,.25)*(-1 if (j+label)%2 else 1)
            branch=strand(start,axes[label]+fan,rng.uniform(6.,13.),label,.75)
            if j%2==0:
                k=len(branch)//2
                strand(branch[k],axes[label]+fan+rng.uniform(-.18,.18),
                       rng.uniform(3.,7.),label,.4)
    if kind in ('toward','unequal','offset'):
        # Upper-shower daughter continues forward/downward into the lower shower.
        # Its trajectory does not turn back toward its parent axis.
        start=trunks[1][int(10./.25)]
        strand(start,-.40,17.,1,.8)
    data=np.array(rows)
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
    fig,axes=plt.subplots(5,2,figsize=(12,14),layout='constrained')
    palette=np.array(['#1677b8','#e88129'])
    for pair,(kind,title,data,vertex,result,labels) in zip(axes,examples):
        true=data[:,3].astype(int)
        for col,ax in enumerate(pair):
            colors=palette[true] if col==0 else (palette[labels] if result['accepted'] else '#777777')
            ax.scatter(data[:,0],data[:,1],c=colors,s=np.clip(9*data[:,2],3,20),linewidths=0)
            if col==1 and result['accepted'] and kind!='single':
                wrong=labels!=true
                ax.scatter(data[wrong,0],data[wrong,1],s=23,facecolors='none',edgecolors='#ce2363',linewidths=.8)
            ax.scatter(0,0,c='black',marker='*',s=95,zorder=4)
            if np.any(vertex):ax.scatter(*vertex,c='#ce2363',marker='x',s=60,zorder=5)
            status='No split' if not result['accepted'] else ('FALSE SPLIT' if kind=='single' else f"Correct charge: {result['charge_accuracy']:.1%}")
            ax.set_title(title+'\n'+('Synthetic truth' if col==0 else status),fontsize=10)
            ax.set_xlim(-2,54);ax.set_ylim(-15,15);ax.set_aspect('equal');ax.grid(alpha=.15)
            ax.set_xlabel('Wire coordinate [cm]');ax.set_ylabel('Drift coordinate [cm]')
    fig.suptitle('Shared-vertex showers with multiple branches: truth (left), splitter output (right)\n'
                 'Pink rings: wrong ownership; star: true vertex; x: assumed vertex. First seed shown.\n'
                 'Actual C++ splitter, unchanged defaults; geometric toys, no detector or 3D simulation.',fontsize=12)
    fig.savefig(args.output/'branching_tests.png',dpi=150);plt.close(fig)
    (args.output/'summary.json').write_text(json.dumps(summaries,indent=2)+'\n')
    (args.output/'trials.json').write_text(json.dumps(trials,indent=2)+'\n')
    print(json.dumps(summaries,indent=2))

if __name__=='__main__':main()
