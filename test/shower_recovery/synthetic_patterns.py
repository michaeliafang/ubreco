"""Toy geometry stress study; NOT a detector simulation or physics efficiency estimate.
Run with the compiled vertex_split_driver and an output directory.
"""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CASES = [
    ('bridge', 'Two widening showers + thin bridge', 2),
    ('unequal', 'Unequal showers + missing-hit bands', 2),
    ('crossing', 'Branch crosses into other shower', 2),
    ('overlap', 'Two showers with narrower opening', 2),
    ('offset', 'Crossing branch; vertex off by 2 cm', 2),
    ('narrow', 'Single narrow shower', 1),
    ('broad', 'Single broad shower', 1),
    ('late', 'Single shower + late branch', 1),
    ('early', 'Single shower + early branch', 1),
]


def generate(kind, seed):
    rng = np.random.default_rng(seed)
    records = []
    def add(xy, label, branch=False):
        q = rng.lognormal(0., .55, len(xy))
        records.extend((x, y, z, label, int(branch)) for (x, y), z in zip(xy, q))
    def shower(deg, length, n, width=.075, label=0):
        a = np.deg2rad(deg)
        u = np.array([np.cos(a), np.sin(a)])
        v = np.array([-u[1], u[0]])
        # Continuous narrow backbone plus a widening, unevenly charged cloud.
        r = np.arange(5.5, length, .3)
        add(r[:, None]*u + rng.normal(0., .045, (len(r), 1))*v, label)
        r = rng.uniform(6., length, n)
        transverse = rng.normal(0., .12+width*(r-5.))
        add(r[:, None]*u + transverse[:, None]*v, label)
    def branch(start, end, label=0):
        start, end = np.array(start), np.array(end)
        n = int(np.linalg.norm(end-start)/.22)+1
        t = np.linspace(0., 1., n)
        xy = start + t[:, None]*(end-start)
        add(xy+rng.normal(0., .035, xy.shape), label, True)
        for _ in range(2):
            add(xy+rng.normal(0., .18, xy.shape), label, True)
    vertex = np.zeros(2)
    if kind in ('bridge', 'unequal', 'crossing', 'offset', 'overlap'):
        angle = 7. if kind == 'overlap' else 19.
        shower(-angle, 42., 650, label=0)
        shower(angle, 28. if kind == 'unequal' else 39.,
               280 if kind == 'unequal' else 600, label=1)
        y = 20.*np.tan(np.deg2rad(angle))
        branch((20., -y), (23., 23.*np.tan(np.deg2rad(angle))))
        if kind in ('crossing','offset'):
            branch((23., 23.*np.tan(np.deg2rad(angle))), (30., 17.))
        if kind == 'offset': vertex = np.array([0., 2.])
    else:
        shower(0., 42., 850, width=.035 if kind == 'narrow' else .11)
        if kind == 'late': branch((24., 0.), (39., 16.))
        if kind == 'early': branch((9., 0.), (32., 19.))
    data = np.array(records)
    # Approximate wire sampling; no electronics response, RMS, or induction model.
    data[:, 0] = np.round(data[:, 0]/.3)*.3
    if kind == 'unequal':
        keep = ~(((data[:,0] > 12.) & (data[:,0] < 12.6)) |
                 ((data[:,0] > 30.) & (data[:,0] < 30.6)))
        data = data[keep]
    return data, vertex


def components(xy, radius=.8):
    # Describe connectivity rather than pretend every supplied toy is one
    # connected component of the original (different) proximity algorithm.
    parent = list(range(len(xy)))
    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    cells = {}
    for i, point in enumerate(xy):
        cell = tuple(np.floor(point/radius).astype(int))
        for dx in (-1,0,1):
            for dy in (-1,0,1):
                for j in cells.get((cell[0]+dx,cell[1]+dy), ()):
                    if np.sum((point-xy[j])**2) <= radius**2:
                        parent[root(i)] = root(j)
        cells.setdefault(cell, []).append(i)
    ids = np.array([root(i) for i in range(len(xy))])
    _, counts = np.unique(ids, return_counts=True)
    return len(counts), int(max(counts))/len(xy)


def evaluate(binary, data, vertex, true_count):
    translated = data[:,:3].copy()
    translated[:,:2] -= vertex
    payload = str(len(data))+'\n'+'\n'.join(' '.join(map(str,row)) for row in translated)+'\n'
    proc = subprocess.run([str(binary)], input=payload, capture_output=True, text=True, check=True)
    lines = proc.stdout.splitlines()
    valid, cut, opening = lines[0].split()
    labels = np.array([int(s) for s in lines[1:]])
    assert len(labels) == len(data)
    valid = bool(int(valid))
    result = dict(accepted=valid, cut_fraction=float(cut) if valid else None,
                  opening_degrees=float(opening), hit_accuracy=None,
                  charge_accuracy=None, branch_accuracy=None, min_shower_recall=None)
    if valid:
        assert set(labels) == {0,1}
        if true_count == 2:
            truth = data[:,3].astype(int)
            # Pick permutation by correct charge, then use it for every metric.
            if data[labels == truth,2].sum() < data[1-labels == truth,2].sum(): labels = 1-labels
            correct = labels == truth
            result['hit_accuracy'] = float(correct.mean())
            result['charge_accuracy'] = float(data[correct,2].sum()/data[:,2].sum())
            branch = data[:,4] == 1
            result['branch_accuracy'] = float(correct[branch].mean())
            result['min_shower_recall'] = float(min(correct[truth==t].mean() for t in (0,1)))
    return result, labels


def plot_examples(examples, out):
    for page, subset in enumerate((examples[:5], examples[5:])):
        if not subset:
            continue
        fig, axes = plt.subplots(len(subset), 2, figsize=(12, 3.35*len(subset)), squeeze=False)
        for row, (name,data,vertex,result,labels,true_count) in enumerate(subset):
            truth = data[:,3].astype(int)
            for col, colors in enumerate((truth, labels)):
                ax=axes[row,col]
                palette=np.array(['#1677b8','#e88129'])
                c=palette[colors] if np.all(colors>=0) else '#76818c'
                ax.scatter(data[:,0],data[:,1],s=np.clip(data[:,2]*5,2,22),c=c,alpha=.65,linewidths=0)
                ax.scatter([0],[0],marker='*',s=100,c='black',label='True vertex',zorder=5)
                if np.any(vertex): ax.scatter(*vertex,marker='x',s=65,c='#bd145b',label='Used vertex',zorder=6)
                if col==1 and result['accepted'] and true_count==2:
                    wrong=labels!=truth
                    ax.scatter(data[wrong,0],data[wrong,1],s=15,facecolors='none',edgecolors='#c01858',linewidths=.6)
                ax.set_aspect('equal',adjustable='box'); ax.set_xlim(-2,46); ax.set_ylim(-20,23)
                ax.set_xlabel('Wire coordinate [cm]'); ax.set_ylabel('Drift coordinate [cm]')
                ax.grid(alpha=.15)
            axes[row,0].set_title(name+' — assigned truth',fontsize=11)
            title='No split' if not result['accepted'] else ('FALSE SPLIT (one true shower)' if true_count==1 else
                f"Split: {result['charge_accuracy']:.1%} charge assigned correctly")
            axes[row,1].set_title(title,fontsize=11)
        fig.suptitle('Synthetic stress tests — actual C++ splitter, unchanged default settings\n'
                     'Pink rings: wrong shower assignment. Marker area: synthetic charge. Not detector simulation.',fontsize=12)
        fig.tight_layout(rect=(0,0,1,.96))
        fig.savefig(out/f'synthetic_examples_{page+1}.png',dpi=150)
        plt.close(fig)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--driver',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--seeds',type=int,default=20)
    args=parser.parse_args(); args.output.mkdir(parents=True,exist_ok=True)
    rows=[]; examples=[]; summaries=[]
    for kind,name,ntruth in CASES:
        case_rows=[]
        for seed in range(args.seeds):
            data,vertex=generate(kind,1000+seed)
            result,labels=evaluate(args.driver,data,vertex,ntruth)
            nc,largest=components(data[:,:2])
            result.update(case=kind,seed=1000+seed,n_hits=len(data),components=nc,largest_component_fraction=largest)
            rows.append(result); case_rows.append(result)
            if seed==0:
                examples.append((name,data,vertex,result,labels,ntruth))
                np.savetxt(args.output/f'{kind}_example.csv',np.column_stack((data,labels)),delimiter=',',
                           header='w_cm,t_cm,charge,truth_shower,is_branch,reco_partition',comments='')
        accepted=[r for r in case_rows if r['accepted']]
        summary=dict(case=kind,name=name,true_showers=ntruth,trials=args.seeds,accepted=len(accepted),
                     mean_largest_component_fraction=float(np.mean([r['largest_component_fraction'] for r in case_rows])))
        if ntruth==2 and accepted:
            for key in ('charge_accuracy','hit_accuracy','branch_accuracy','min_shower_recall'):
                summary['mean_'+key]=float(np.mean([r[key] for r in accepted]))
        summaries.append(summary)
    with (args.output/'results.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    (args.output/'summary.json').write_text(json.dumps(summaries,indent=2)+'\n')
    plot_examples(examples,args.output)
    print(json.dumps(summaries,indent=2))

if __name__=='__main__': main()
