# Fortran DDS on Athabasca: figures (best KGE so far, best parameter sets, run time vs heightCanopyTop)
# and the run-time analysis in PROGRESS.md (Nov 10-12).
# Usage: python analyze_fortran_dds.py <folder with one subfolder per run: results/evaluations.csv, dds_run.log> [figure folder]
import re,json,numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
import sys
D=sys.argv[1] if len(sys.argv)>1 else '.'; O=sys.argv[2] if len(sys.argv)>2 else '.'
BLUE,ORANGE,INK,MUTED,GRID='#2a78d6','#eb6834','#0b0b0b','#52514e','#e4e3df'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.edgecolor':MUTED,'axes.labelcolor':INK,
 'xtick.color':MUTED,'ytick.color':MUTED,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,
 'grid.color':GRID,'grid.linewidth':0.6})
lab={'full':'Original s42','seed43':'s43','seed44':'s44','seed45':'s45','seed43b':'43b','seed44b':'44b','seed45b':'45b',
     'defstart_s42':'d42','defstart_s43':'d43','defstart_s44':'d44'}
# Fig 1 curves from logs
fig,ax=plt.subplots(figsize=(6.5,3.6))
ends=[]
for d,l in lab.items():
    v=[float(m.group(1)) for m in re.finditer(r'objective= *([-0-9.Ee+]+)',open(f'{D}/{d}/dds_run.log').read())]
    b=np.maximum.accumulate(v); x=np.arange(1,len(b)+1)
    dflt=d.startswith('def'); c=BLUE if dflt else ORANGE
    ls='--' if d.endswith('b') else '-'
    ax.plot(x,b,color=c,lw=1.6 if dflt else 1.2,ls=ls)
    ends.append((b[-1],x[-1],l,c))
ax.axhline(0.8925,color=MUTED,lw=1,ls=':'); ax.text(150,0.8945,'async DDS v2 (0.8925)',color=MUTED,fontsize=8)
ends.sort()
last=-1
for y,x,l,c in ends:
    yy=max(y,last+0.0085); last=yy
    ax.text(5180,yy,f'{l} {y:.3f}',va='center',fontsize=7.5,color=INK)
ax.set_xlim(0,5150); ax.set_ylim(0.65,0.92); ax.set_xlabel('Model runs completed'); ax.set_ylabel('Best KGE so far')
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([],[],color=BLUE,lw=1.6,label='Default start'),Line2D([],[],color=ORANGE,lw=1.2,label='Random start'),
  Line2D([],[],color=ORANGE,lw=1.2,ls='--',label='Random start, rerun')],loc='lower right',frameon=False,fontsize=8)
fig.tight_layout(); fig.subplots_adjust(right=0.84); fig.savefig(f'{O}/fig1_curves.png',dpi=200)
# Fig 2 run time vs canopy
runs=['full','seed43b','seed44b','seed45b','defstart_s42','defstart_s43','defstart_s44']
a=pd.concat([pd.read_csv(f'{D}/{d}/results/evaluations.csv').assign(job=d) for d in runs])
fig,ax=plt.subplots(figsize=(6.5,3.2))
fast=a.job.isin(['full','defstart_s42','defstart_s43','defstart_s44'])
ax.scatter(a[fast].heightCanopyTop,a[fast].seconds/60,s=4,color=BLUE,alpha=.25,lw=0,label='Fast jobs (original, d42–d44)')
ax.scatter(a[~fast].heightCanopyTop,a[~fast].seconds/60,s=4,color=ORANGE,alpha=.25,lw=0,label='Slow jobs (43b–45b)')
ax.set_xscale('log'); ax.set_ylim(5,40); ax.set_xlabel('heightCanopyTop (m, log scale)'); ax.set_ylabel('Run time (min)')
ax.axvline(20,color=MUTED,lw=1,ls=':'); ax.text(21,37,'20 m (default)',fontsize=8,color=MUTED)
lg=ax.legend(loc='upper left',frameon=False,fontsize=8,markerscale=3)
for h in lg.legend_handles: h.set_alpha(1)
fig.tight_layout(); fig.savefig(f'{O}/fig2_runtime.png',dpi=200)
# Fig 3 best params normalized
P=list(a.columns[6:21]); logp=['k_soil','aquiferBaseflowRate','qSurfScale','heightCanopyTop','routingGammaScale']
tr=lambda p,v: np.log10(v) if p in logp else v
lo={p:tr(p,a[p]).min() for p in P}; hi={p:tr(p,a[p]).max() for p in P}
best={d:a[a.job==d].loc[a[a.job==d].kge.idxmax()] for d in runs}
order=['qSurfScale','routingGammaScale','routingGammaShape','frozenPrecipMultip','aquiferBaseflowRate','aquiferBaseflowExp','heightCanopyTop',
 'theta_sat','k_soil','vGn_n','Fcapil','summerLAI','heightCanopyBottom','windReductionParam','tempCritRain']
fig,ax=plt.subplots(figsize=(6.5,4.6))
yl=np.arange(len(order))[::-1]
mk={'full':'v','seed43b':'s','seed44b':'^','seed45b':'D'}
for d in runs:
    b=best[d]; xs=[(tr(p,b[p])-lo[p])/(hi[p]-lo[p]) for p in order]
    dflt=d.startswith('def')
    ax.scatter(xs,yl+(0 if dflt else 0),s=30 if dflt else 34,color=BLUE if dflt else ORANGE,marker='o' if dflt else mk[d],
       edgecolor='white',lw=0.8,zorder=3 if dflt else 4,alpha=0.9,
       label=(f'{lab[d]} ({b.kge:.3f})'))
ax.set_yticks(yl); ax.set_yticklabels([p+(' (log)' if p in logp else '') for p in order],fontsize=8)
ax.set_xlim(-0.03,1.03); ax.set_xlabel('Position within sampled range (0 = lower bound, 1 = upper bound)')
ax.grid(axis='y',color=GRID); ax.axhspan(yl[6]-0.5,yl[0]+0.5,color='#f3f2ee',zorder=0)
ax.legend(loc='upper center',bbox_to_anchor=(0.5,-0.13),ncol=4,frameon=False,fontsize=7.5)
fig.tight_layout(); fig.savefig(f'{O}/fig3_params.png',dpi=200)
# --- run-time analysis ---
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_val_score
X=np.column_stack([tr(p,a[p]) for p in P]); y=np.log(a.seconds.values)
print('CV R2, params -> log run time: %.3f'%cross_val_score(RandomForestRegressor(200,n_jobs=-1,random_state=0),X,y,cv=4).mean())
rf=RandomForestRegressor(300,n_jobs=-1,random_state=0).fit(X,y)
for i,p in sorted(zip(rf.feature_importances_,P),reverse=True)[:4]: print('  %-22s %.3f'%(p,i))
for d in runs:
    g=a[a.job==d]; lo_,hi_=g[g.heightCanopyTop<20],g[g.heightCanopyTop>=20]
    print(f'{d:13s} <20m n={len(lo_):5d} med={lo_.seconds.median()/60:5.1f}  >=20m n={len(hi_):5d} med={hi_.seconds.median()/60:5.1f}')
