#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Author: Anonymous
"""Generate the explanatory figure. Usage: python3 make_figure.py OUTDIR.
Requires numpy and matplotlib; the curves are illustrations, not proofs."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from legacy_checks import HOUSE_EDGES, house_formula

plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                     'figure.facecolor':'white','savefig.facecolor':'white'})
fig,axs=plt.subplots(2,2,figsize=(11.2,6.3))
ax=axs[0,0]
pos={0:(0,1),1:(1.4,1),2:(1.4,0),3:(0,0),4:(.7,1.85)}
for i,j in HOUSE_EDGES:ax.plot([pos[i][0],pos[j][0]],[pos[i][1],pos[j][1]],color='#B9ADB5',lw=2.3,zorder=1)
for i,(x,y) in pos.items():
 color='#70365F' if i==4 else '#167C80' if i in (0,1) else '#B48328'
 ax.scatter([x],[y],s=430,color=color,zorder=2)
 ax.text(x,y,str(i),ha='center',va='center',color='white',weight='bold',fontsize=12)
ax.text(1.85,1.65,'Roof: vertex 4\nDegree 2',color='#70365F',va='center')
ax.text(1.85,.75,'Vertices 0, 1: degree 3\nVertices 2, 3: degree 2',va='center')
ax.set_xlim(-.35,3.4);ax.set_ylim(-.25,2.15);ax.set_aspect('equal');ax.axis('off')
ax.set_title('(a) The house graph',loc='left',fontweight='bold')
for ax,e,upper,col,letter in [(axs[0,1],1/3,1,'#70365F','b'),(axs[1,0],.5,10,'#B48328','c'),(axs[1,1],.85,60,'#167C80','d')]:
 mus=np.linspace(0,upper,600)
 ps=np.array([float(house_formula(mu,e)) for mu in mus]);standard=float(house_formula(1,e))
 ax.plot(mus,1e4*(ps-standard),color=col,lw=2.2)
 ax.axhline(0,color='#A497A0',lw=.8,ls='--')
 ax.axvline(1,color='#A497A0',lw=.8,ls=':')
 ax.set_xlabel('Backtracking weight, μ')
 ax.set_ylabel('10⁴ × [π_μ(4) − π_1(4)]')
 label='1/3' if e==1/3 else '1/2' if e==.5 else '17/20'
 ax.set_title(f'({letter}) Damping ε = {label}',loc='left',fontweight='bold')
 ax.grid(axis='y',color='#EEE8ED',lw=.7)
fig.tight_layout(pad=1.7)
root=Path(sys.argv[1]); root.mkdir(parents=True, exist_ok=True)
fig.savefig(root/'house_response.pdf',bbox_inches='tight',metadata={'CreationDate':None,'ModDate':None})
fig.savefig(root/'house_response.png',dpi=180,bbox_inches='tight')
