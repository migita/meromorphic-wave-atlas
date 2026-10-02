"""Exact affine-Weierstrass construction for the even C4 equation.

v'''' + A*v*v'' + B*(v')² + C*v³ + a2*v'' + a0*v + a00 = 0.
The displayed sector is v=alpha*wp+beta; no multiple-pole/genus-two
completeness claim is made. All equations below are derived directly.
"""
from pathlib import Path
import json
import numpy as np
import sympy as s
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from render import PAPER,PANEL,INK,MUTED,GRID,CMAP

HERE=Path(__file__).resolve().parent
X,g2,g3,A,B,C,alpha,beta,a2,a0,a00=s.symbols('X g2 g3 A B C alpha beta a2 a0 a00')
R=s.Rational
PRESETS=[('Lax core',(10,5,10),(-2,.25)),('Sawada–Kotera',(15,0,15),(-3,.25)),
         ('Kaup–Kupershmidt',(10,R(15,2),R(20,3)),(-8,.5)),('Pure cubic',(0,0,-30),(-.6,.6)),
         ('Generic core',(5,30,30),(-3.2,.4))]


def residual():
    cubic=4*X**3-g2*X-g3
    D2=lambda f:s.expand(s.diff(f,X,2)*cubic+s.diff(f,X)*(6*X**2-g2/2))
    v=alpha*X+beta
    return s.Poly(s.expand(D2(D2(v))+A*v*D2(v)+B*alpha**2*cubic+C*v**3+a2*D2(v)+a0*v+a00),X)


def readout(aa,bb,cc,al,bt=None):
    aa,bb,cc,al=map(s.sympify,(aa,bb,cc,al))
    bt=s.cancel(-2*a2/(2*aa+cc*al)) if bt is None else bt
    G2=s.cancel((6*cc*bt**2+2*a0)/(36+(aa+2*bb)*al))
    G3=s.cancel((-al*(aa*bt+a2)*G2+2*cc*bt**3+2*a0*bt+2*a00)/(2*al*(12+bb*al)))
    return bt,G2,G3


def verify():
    res=residual();lead=C*alpha**2+(6*A+4*B)*alpha+120
    expected=[alpha*lead,3*alpha*((2*A+C*alpha)*beta+2*a2),
        alpha*(6*C*beta**2+2*a0-(36+(A+2*B)*alpha)*g2)/2,
        (-alpha*(A*beta+a2)*g2+2*C*beta**3+2*a0*beta+2*a00-2*alpha*(12+B*alpha)*g3)/2]
    assert all(s.expand(u-v)==0 for u,v in zip(res.all_coeffs(),expected))
    bt,G2,G3=readout(A,B,C,alpha)
    for j in (2,1,0):assert s.cancel(res.nth(j).subs({beta:bt,g2:G2,g3:G3}))==0
    checks=[]
    for name,(aa,bb,cc),interval in PRESETS:
        for al in s.solve(lead.subs({A:aa,B:bb,C:cc}),alpha):
            if 2*aa+cc*al==0:
                bt,G2,G3=readout(aa,bb,cc,al,beta)
                equation=res.as_expr().subs({A:aa,B:bb,C:cc,alpha:al,a2:0,g2:G2.subs(a2,0),g3:G3.subs(a2,0)})
                assert s.cancel(equation)==0
                checks.append(dict(core=name,alpha=str(al),condition='a2=0',free='beta',identity='exact zero'))
            else:
                bt,G2,G3=readout(aa,bb,cc,al)
                assert s.cancel(res.as_expr().subs({A:aa,B:bb,C:cc,alpha:al,beta:bt,g2:G2,g3:G3}))==0
                checks.append(dict(core=name,alpha=str(al),identity='exact zero'))
    special=[
      ({A:1,B:1,C:0,alpha:-12,beta:-a2,a0:0,a00:0},'g2,g3 free'),
      ({A:0,B:1,C:-R(4,27),alpha:-18,beta:-R(3,4),a2:1,a0:R(1,4),a00:0,g3:g2/12-R(1,864)},'g2 free'),
      ({A:0,B:1,C:-R(1,2),alpha:-12,beta:-R(1,3),a2:1,a0:0,a00:R(4,27),g2:-R(1,36)},'g3 free'),
    ]
    for be in [s.Integer(0),R(1,2),-R(1,2)]:
        special.append(({A:3,B:2,C:1,alpha:-6,beta:be,a2:0,a0:-1,a00:0,g2:(2-6*be**2)/6},'g3 free; constrained beta'))
    for be in [R(1,2),-R(1,2)]:
        special.append(({A:4,B:1,C:R(4,3),alpha:-6,beta:be,a2:0,a0:-1,a00:0,g3:-be*g2/3+be/54},'g2 free; constrained beta'))
    for sub,note in special:
        assert s.expand(res.as_expr().subs(sub))==0
        checks.append(dict(stratum=note,substitutions={str(k):str(v) for k,v in sub.items()},identity='exact zero'))
    # The Lax one-gap pulse belongs to the same fixed equation as its curve.
    T=s.Symbol('T');v=(1-T*T)/2
    D=lambda f:s.expand((1-T*T)*s.diff(f,T)/2)
    assert s.expand(D(D(D(D(v))))+10*v*D(D(v))+5*D(v)**2+10*v**3-v)==0
    data=dict(equation=__doc__,sector='v=alpha*wp+beta',undivided_conditions=[str(s.factor(c)) for c in expected],
        generic_readout={'beta':str(readout(A,B,C,alpha)[0]),'g2':str((6*C*beta**2+2*a0)/(36+(A+2*B)*alpha)),
        'g3':str((-alpha*(A*beta+a2)*g2+2*C*beta**3+2*a0*beta+2*a00)/(2*alpha*(12+B*alpha)))},
        status='PASS',generic_identity=True,checks=checks,lax_pulse='v=sech(z/2)^2/2, a2=0,a0=-1,a00=0',pulse_identity=True)
    (HERE/'data/c4_formulas.json').write_text(json.dumps(data,indent=2,ensure_ascii=False)+'\n')
    return data


def figure():
    fig,axes=plt.subplots(2,3,figsize=(15,10.4));axes=axes.ravel()
    fig.suptitle('C4 · affine elliptic waves and a free lattice curve',x=.055,y=.96,ha='left',fontsize=24)
    fig.text(.055,.903,r"$v''''+A v v''+B(v')^2+C v^3+a_2v''+a_0v+a_{00}=0,\qquad v=\alpha\wp+\beta$",fontsize=15,color=MUTED)
    records=[]
    for i,(name,core,interval) in enumerate(PRESETS+[('Lax · one fixed equation',(10,5,10),(-.21,.21))]):
        ax=axes[i];aa,bb,cc=core;var=beta if i==5 else a0
        branches=[]
        if i==5:
            branches=[(-2,R(1,2)-15*beta**2,35*beta**3-beta)]
        else:
            for al in s.solve(cc*alpha**2+(6*aa+4*bb)*alpha+120,alpha):
                if 2*aa+cc*al==0:continue
                if name=='Pure cubic' and al<0:continue
                bt,G2,G3=readout(aa,bb,cc,al)
                branches.append((al,G2.subs({a2:1,a00:0}),G3.subs({a2:1,a00:0})))
        norm=Normalize(*interval);allpts=[];branchdata=[]
        for j,(al,G2,G3) in enumerate(branches):
            values=np.linspace(*interval,1801)
            f=s.lambdify(var,(G2,G3),'numpy');x,y=f(values);xy=np.column_stack((np.broadcast_to(x,values.shape),np.broadcast_to(y,values.shape)))
            allpts.append(xy);lc=LineCollection(np.stack((xy[:-1],xy[1:]),axis=1),array=values[:-1],cmap=CMAP,norm=norm,linewidth=2.0,zorder=3)
            if j:lc.set_linestyle((0,(3,1.3)))
            ax.add_collection(lc)
            delta=s.factor(G2**3-27*G3**2);roots=s.Poly(delta,var).intervals(eps=R(1,10**18))
            nodes=[]
            for bounds,_ in roots:
                root=float((bounds[0]+bounds[1])/2)
                if interval[0]<=root<=interval[1]:
                    Xv,Yv=map(float,f(root));ax.plot(Xv,Yv,'D',ms=4.8,color='#d86849',mec=PANEL,mew=.9,zorder=5);nodes.append(root)
            branchdata.append(dict(alpha=str(al),g2=str(G2),g3=str(G3),degenerations=nodes))
        points=np.vstack(allpts);lo=points.min(axis=0);hi=points.max(axis=0);span=hi-lo
        ax.set_xlim(lo[0]-.07*span[0],hi[0]+.07*span[0]);ax.set_ylim(lo[1]-.1*span[1],hi[1]+.1*span[1])
        xmax=ax.get_xlim()[1]
        if xmax>0:
            x=np.linspace(max(ax.get_xlim()[0],0),xmax,601);y=np.sqrt(x**3/27)
            ax.fill_between(x,-y,y,color='#edf3e8',zorder=0)
            for sign in [-1,1]:ax.plot(x,sign*y,color='#929ba3',lw=1,zorder=1)
        if i==5:ax.plot(1/42,0,'o',color=INK,ms=6,mec=PANEL,zorder=5);ax.annotate('α = −6',(1/42,0),xytext=(-8,-17),textcoords='offset points',fontsize=8,ha='right')
        ax.grid(color=GRID,lw=.6);ax.set_axisbelow(True)
        for side in ['top','right']:ax.spines[side].set_visible(False)
        ax.set_title(name+('\na₂ = 0, a₀ = −1, a₀₀ = 0' if i==5 else f'\n(A, B, C) = ({aa}, {bb}, {cc})'),loc='left',fontsize=11)
        ax.set_xlabel('g₂',fontsize=10);ax.set_ylabel('g₃',fontsize=10);ax.tick_params(labelsize=8)
        ax.ticklabel_format(style='sci',scilimits=(-2,2),axis='both');ax.locator_params(nbins=4)
        records.append(dict(name=name,core=list(map(str,core)),parameter=str(var),interval=interval,branches=branchdata))
    fig.subplots_adjust(left=.065,right=.975,top=.805,bottom=.12,wspace=.34,hspace=.58)
    fig.text(.055,.06,'Panels 1–5: a₂ = 1, a₀₀ = 0; a₀ varies. Panel 6: one fixed Lax equation, with β free. Diamonds mark Δ = 0.',fontsize=10,color=MUTED)
    fig.text(.055,.034,'The explorer lets you edit all six equation coefficients, select the sweep coordinate, focus a branch, and inspect its profile. Only the affine-℘ sector is shown.',fontsize=9.5,color=MUTED)
    (HERE/'data/c4_presets.json').write_text(json.dumps(records,indent=2,ensure_ascii=False)+'\n')
    return fig


def main():
    data=verify();fig=figure()
    for ext in ['png','svg','pdf']:fig.savefig(HERE/'figures'/('c4.'+ext),dpi=170)
    plt.close(fig)
    print('PASS:',len(data['checks']),'symbolic family/stratum checks, generic coefficient identity, and the Lax sech² pulse.')


if __name__=='__main__':main()
