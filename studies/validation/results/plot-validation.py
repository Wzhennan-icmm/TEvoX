import json,pathlib,re,collections
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=pathlib.Path('/workspace/tevox-v05-validation-20261006')
out=pathlib.Path('/workspace/tevox-review-20261003/early-development/v0.5-source/studies/validation/results')
poly=json.loads((root/'polyploid/summary.json').read_text())['runs']
scale=json.loads((root/'scaling/summary.json').read_text())['runs']
values=collections.defaultdict(list)
for r in poly:
 path=next(s for s in r['command'] if s.endswith('manifest.tsv'))
 copies,missing,seed=map(int,re.search(r'copies(\d+)_missing(\d+)_seed(\d+)',path).groups())
 values[(r['condition'],missing,copies)].append(r['metrics']['b_cubed_f1'])
fig,ax=plt.subplots(1,3,figsize=(12,3.7),layout='constrained')
for condition,color,name in [('with_synteny','#0072B2','With prior'),('without_synteny','#D55E00','Without prior')]:
 for missing,style in [(0,'-'),(1,'--')]:
  y=[sum(values[(condition,missing,c)])/len(values[(condition,missing,c)]) for c in [2,4,6]]
  ax[0].plot([2,4,6],y,style,marker='o',color=color,label=name+(' / missing annotation' if missing else ' / complete'))
ax[0].set(ylim=(0,1.025),xticks=[2,4,6],xlabel='Retained homeologous copies',ylabel='B-cubed F1',title='Controlled copy-context ablation')
ax[0].legend(fontsize=6.8,loc='lower left',frameon=False)
x=[r['counts']['tes'] for r in scale]
for panel,key,factor,label,title,color in [(ax[1],'wall_seconds',1,'Wall time (s)','Sparse synthetic scaling','#0072B2'),(ax[2],'peak_rss_kib',1024**2,'Peak RSS (GiB)','Core + gzip resource use','#009E73')]:
 panel.loglog(x,[r[key]/factor for r in scale],marker='o',color=color)
 panel.set(xlabel='Input TE annotation fragments',ylabel=label,title=title)
for panel in ax:
 panel.grid(alpha=.18);panel.spines['top'].set_visible(False);panel.spines['right'].set_visible(False)
fig.suptitle('TEvoX frozen integration build: synthetic validation only',fontsize=12)
fig.savefig(out/'validation-summary.png',dpi=180)
fig.savefig(out/'validation-summary.pdf')
print('saved',out/'validation-summary.png')
