"""Kawahara's regular real branch in wavelength--amplitude--mean coordinates.

Equation: v'''' + a v'' - v + v^2/2 = 0, -13/6 < a < 13/6.
The larger g2 root has Delta > 0. Its pole-free real oval is parameterized
by Jacobi sn, so both the profile and its physical observables are explicit.
The infinite-period endpoint is a limit, never a finite-wavelength point.
"""
from pathlib import Path
import csv
import json
import math
import sys

import mpmath as mp
import numpy as np
import sympy as s
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
from scipy.special import ellipj
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import plotly.graph_objects as go
from plotly.subplots import make_subplots

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from render import PAPER,PANEL,INK,MUTED,GRID,CMAP

AC=13/6
NORM=Normalize(-AC,AC)
SMALL_PERIOD=2*math.pi*math.sqrt(6)
PULSE_AMPLITUDE=35/12
PULSE_DEFICIT=70*math.sqrt(6)/9


def state(a=None,*,side=None,distance=None,dps=65):
    """Use extra precision before subtracting near-coincident lattice roots."""
    if distance is not None:
        dps=max(dps,int(-math.log10(float(distance)))+55)
    with mp.workdps(dps):
        critical=mp.mpf(13)/6
        av=mp.mpf(str(a)) if side is None else side*(critical-mp.mpf(str(distance)))
        assert -critical<av<critical
        rad=mp.sqrt(653016**2*av**4-4*662158224*(1457*av**4-28561))
        g2=(653016*av**2+rad)/(2*662158224)
        g3=av*(31*av**2-42588*g2)/4745520
        radius=mp.sqrt(g2/12)
        theta=mp.acos(3*mp.sqrt(3)*g3/g2**mp.mpf('1.5'))/3
        e1=2*radius*mp.cos(theta)
        e2=2*radius*mp.cos(theta-2*mp.pi/3)
        e3=2*radius*mp.cos(theta+2*mp.pi/3)
        G=e1-e3; m=(e2-e3)/G
        K=mp.ellipk(m); E=mp.ellipe(m)
        period=2*K/mp.sqrt(G)
        background=-g2/10-(31*av**2+507)/851760
        def V(X):return -1680*(X**2+av*X/78+background)
        extrema=[V(e3),V(e2)]
        vertex=-av/156
        if e3<vertex<e2:extrema.append(V(vertex))
        low,high=min(extrema),max(extrema)
        average_X=e1-G*E/K
        mean=-1680*(g2/12+av*average_X/78+background)
        rms=mp.sqrt(mean*(2-mean))
        values=dict(a=av,g2=g2,g3=g3,e1=e1,e2=e2,e3=e3,G=G,m=m,
            wavelength=period,amplitude=high-low,mean=mean,minimum=low,maximum=high,rms=rms,
            background=background,distance_to_endpoint=critical-abs(av))
        result={key:float(value) for key,value in values.items()}
        result["a_exact_decimal"]=mp.nstr(av,dps-10)
        result["g2_exact_decimal"]=mp.nstr(g2,dps-10)
        result["g3_exact_decimal"]=mp.nstr(g3,dps-10)
        result["m_exact_decimal"]=mp.nstr(m,dps-10)
        result["dps"]=dps
        if side is not None and float(distance)<1e-7:
            result["a_label"]=f"{'−' if side<0 else ''}13/6 {'+' if side<0 else '−'} {float(distance):.1e}"
        else:result["a_label"]=f"{float(av):.5g}"
        return result


def profile(record,z):
    """Place a deepest trough at zero; reduce to the fundamental real cell."""
    L=record["wavelength"]
    z=(np.asarray(z)+L/2)%L-L/2
    sn,cn,dn,phase=ellipj(np.sqrt(record["G"])*z,record["m"])
    X=record["e3"]+(record["e2"]-record["e3"])*sn**2
    return -1680*(X**2+record["a"]*X/78+record["background"])


def pulse(z):
    return 2-(35/12)/np.cosh(np.asarray(z)/(2*np.sqrt(6)))**4


def build_data():
    records=[state(a=a) for a in np.linspace(-AC+1e-4,AC-1e-4,1001)]
    # Near the pulse the period grows logarithmically while a becomes extremely
    # close to -13/6. Retain the distance and a high-precision decimal in the CSV.
    records.extend(state(side=-1,distance=10.**(-k)) for k in np.linspace(4,80,240))
    records.extend(state(side=1,distance=10.**(-k)) for k in np.linspace(4,20,70))
    # Sorting by double precision a would identify a long sequence of distinct
    # long-period states. The order in the tail is explicitly retained.
    negative=[r for r in records if r["a"]<0]
    negative.sort(key=lambda r:r["distance_to_endpoint"])
    positive=[r for r in records if r["a"]>=0]
    positive.sort(key=lambda r:r["a"] if r["distance_to_endpoint"]>1e-10 else AC)
    near_positive=[r for r in positive if r["distance_to_endpoint"]<=1e-10]
    near_positive.sort(key=lambda r:r["distance_to_endpoint"],reverse=True)
    positive=[r for r in positive if r["distance_to_endpoint"]>1e-10]+near_positive
    return negative+positive


def validate(records):
    a,x,g2,g3=s.symbols("a x g2 g3")
    Q=4*x**3-g2*x-g3
    D2=lambda f:s.expand(s.diff(f,x,2)*Q+s.diff(f,x)*(6*x**2-g2/2))
    V=-1680*(x**2+a*x/78-g2/10-(31*a*a+507)/851760)
    residual=s.expand(D2(D2(V))+a*D2(V)-V+V*V/2)
    line=a*(31*a*a-42588*g2)/s.Integer(4745520)
    condition=662158224*g2**2-653016*a*a*g2+1457*a**4-28561
    for c in s.Poly(residual.subs(g3,line),x).all_coeffs():
        assert s.Poly(c,g2,domain=s.QQ.frac_field(a)).rem(s.Poly(condition,g2,domain=s.QQ.frac_field(a))).is_zero
    T=s.Symbol("T")
    P=2-s.Rational(35,12)*(1-T*T)**2
    D=lambda f:s.expand((1-T*T)*s.diff(f,T)/(2*s.sqrt(6)))
    assert s.expand(D(D(D(D(P))))-s.Rational(13,6)*D(D(P))-P+P*P/2)==0
    quadrature=[]
    for av in [-2.16,-2,-1,0,1,2,2.16]:
        r=state(a=av);L=r["wavelength"]
        avg=quad(lambda zz:float(profile(r,zz)),0,L/2,epsabs=2e-12,epsrel=2e-12)[0]*2/L
        avg2=quad(lambda zz:float(profile(r,zz))**2,0,L/2,epsabs=2e-12,epsrel=2e-12)[0]*2/L
        mean_error=abs(avg-r["mean"]);balance_error=abs(avg2-2*r["mean"])
        z=np.linspace(0,L,30001,endpoint=False);v=profile(r,z)
        amp_error=abs(np.ptp(v)-r["amplitude"])
        assert mean_error<2e-10 and balance_error<2e-10 and amp_error<2e-7
        quadrature.append(dict(a=av,mean_error=mean_error,mean_square_balance_error=balance_error,
                              sampled_extrema_error=amp_error))
    long=state(side=-1,distance=1e-80)
    small=state(side=1,distance=1e-20)
    assert abs(small["wavelength"]-SMALL_PERIOD)<1e-7
    assert small["amplitude"]<1e-8 and abs(small["mean"]-2)<1e-14
    assert abs(long["amplitude"]-PULSE_AMPLITUDE)<1e-12
    mass=(2-long["mean"])*long["wavelength"]
    assert abs(mass-PULSE_DEFICIT)<2e-12
    pulse_error=float(np.max(abs(profile(long,np.linspace(-20,20,2001))-pulse(np.linspace(-20,20,2001)))))
    assert pulse_error<2e-12
    assert all(0<r["mean"]<2+1e-14 and r["amplitude"]>0 for r in records)
    minimum=minimize_scalar(lambda aa:state(a=aa)["wavelength"],bounds=(-2.12,0),method="bounded",options={"xatol":1e-12})
    mean_min=minimize_scalar(lambda aa:state(a=aa)["mean"],bounds=(-2.15,-1),method="bounded",options={"xatol":1e-12})
    report=dict(status="PASS",symbolic_ode_identity=True,exact_pulse_identity=True,quadrature=quadrature,
        shortest_wave=state(a=float(minimum.x)),lowest_mean=state(a=float(mean_min.x)),
        small_amplitude_limit=dict(wavelength=SMALL_PERIOD,amplitude=0,mean=2),
        pulse_limit=dict(wavelength="infinity",amplitude=PULSE_AMPLITUDE,background=2,mean_limit=2,
                         integrated_deficit=PULSE_DEFICIT,profile_check_max_error=pulse_error),
        sampled_period_range=[min(r["wavelength"] for r in records),max(r["wavelength"] for r in records)])
    (HERE/"checks.json").write_text(json.dumps(report,indent=2)+"\n")
    return report


def style(ax):
    ax.grid(color=GRID,lw=.7);ax.set_axisbelow(True)
    for side in ("top","right"):ax.spines[side].set_visible(False)
    ax.tick_params(labelsize=9)


def colored_curve(ax,records,xkey,ykey,lw=2.7):
    xy=np.array([[r[xkey],r[ykey]] for r in records])
    lc=LineCollection(np.stack([xy[:-1],xy[1:]],axis=1),array=np.array([r["a"] for r in records[:-1]]),
                      cmap=CMAP,norm=NORM,linewidth=lw,zorder=3)
    ax.add_collection(lc);ax.autoscale_view();return lc


def add_markers(ax,examples,xkey,ykey,offsets=None):
    for j,r in enumerate(examples):
        x,y=r[xkey],r[ykey]
        ax.scatter([x],[y],s=40,c=[CMAP(NORM(r["a"]))],edgecolor=PANEL,linewidth=1.1,zorder=5)
        dx,dy=offsets[j] if offsets else (5,7)
        ax.annotate(chr(65+j),(x,y),xytext=(dx,dy),textcoords="offset points",fontsize=9,weight="bold",color=INK)


def coordinate_figure(records,examples,report):
    fig=plt.figure(figsize=(15,10.4))
    fig.text(.055,.958,"Kawahara in physical coordinates",fontsize=25,weight="medium")
    fig.text(.055,.912,r"$v''''+a v''-v+v^2/2=0$    ·    smooth real periodic branch    ·    $-13/6<a<13/6$",fontsize=13,color=MUTED)
    gs=fig.add_gridspec(2,3,left=.065,right=.97,bottom=.20,top=.845,hspace=.48,wspace=.34,height_ratios=[1.15,1])
    configs=[("wavelength","amplitude","Wavelength → amplitude","Peak-to-trough amplitude  A"),
             ("wavelength","mean","Wavelength → mean","Period average  M"),
             ("mean","amplitude","Mean → amplitude","Peak-to-trough amplitude  A")]
    for j,(xkey,ykey,title,ylabel) in enumerate(configs):
        ax=fig.add_subplot(gs[0,j]);style(ax);colored_curve(ax,records,xkey,ykey)
        ax.set_title(title,loc="left",fontsize=12,pad=13);ax.set_ylabel(ylabel,fontsize=10)
        if xkey=="wavelength":
            ax.set_xscale("log");ax.set_xlim(11,300)
            ticks=[12,20,40,80,160];ax.set_xticks(ticks,[str(v) for v in ticks]);ax.set_xlabel("Wavelength L  (log scale)",fontsize=10)
            ax.minorticks_off()
        else:
            ax.set_xlim(.42,2.06);ax.set_xlabel("Period average  M",fontsize=10)
        ax.set_ylim((-.12,3.12) if ykey=="amplitude" else (.43,2.12))
        if j==0:
            ax.scatter([SMALL_PERIOD],[0],marker="o",s=55,facecolor=PANEL,edgecolor=INK,zorder=5)
            ax.annotate("small-amplitude limit",(SMALL_PERIOD,0),xytext=(30,.32),arrowprops=dict(arrowstyle="-",color=MUTED),fontsize=8.7,color=MUTED)
            ax.axhline(PULSE_AMPLITUDE,color=MUTED,lw=.8,ls=(0,(4,4)))
            ax.text(38,3.015,r"$A\to35/12$ as $L\to\infty$",fontsize=9,color=MUTED)
            add_markers(ax,examples,xkey,ykey,[(7,-10),(7,-12),(-10,6),(5,-17)])
        elif j==1:
            ax.axhline(2,color=MUTED,lw=.8,ls=(0,(4,4)))
            ax.text(43,2.045,r"$M\to2$ as $L\to\infty$",fontsize=9,color=MUTED)
            add_markers(ax,examples,xkey,ykey,[(7,-4),(7,-4),(7,-5),(6,6)])
        else:
            ax.scatter([2],[0],marker="o",s=55,facecolor=PANEL,edgecolor=INK,zorder=5)
            ax.scatter([2],[PULSE_AMPLITUDE],marker="D",s=55,facecolor=PANEL,edgecolor=INK,zorder=5)
            ax.annotate("pulse limit",(2,PULSE_AMPLITUDE),xytext=(1.55,3.02),fontsize=9,color=MUTED,arrowprops=dict(arrowstyle="-",color=MUTED))
            add_markers(ax,examples,xkey,ykey,[(7,-10),(6,-14),(6,7),(-8,-15)])
    # Three profiles with a common physical z-axis, plus a long-period cell
    # compared directly with the exact solitary pulse.
    z=np.linspace(-24,24,3201)
    for j,r in enumerate(examples[:3]):
        ax=fig.add_subplot(gs[1,j]);style(ax)
        ax.plot(z,profile(r,z),color=CMAP(NORM(r["a"])),lw=2)
        ax.axhline(r["mean"],color=MUTED,lw=.9,ls=(0,(3,3)))
        ax.set_xlim(-24,24);ax.set_ylim(-1.05,2.6)
        half=r["wavelength"]/2
        ax.annotate("",(-half,2.36),(half,2.36),arrowprops=dict(arrowstyle="<->",lw=.8,color=MUTED))
        ax.text(0,2.43,"one full period L",fontsize=8,color=MUTED,ha="center")
        ax.set_title(f"{chr(65+j)}   a = {r['a']:g}   ·   L = {r['wavelength']:.2f}",loc="left",fontsize=10.7,pad=11)
        ax.set_xlabel("Travelling coordinate z",fontsize=10);ax.set_ylabel("Wave profile v(z)",fontsize=10)
        ax.text(.03,.055,f"A = {r['amplitude']:.3f}    M = {r['mean']:.3f}",transform=ax.transAxes,fontsize=9,color=MUTED,
                bbox=dict(facecolor=PANEL,edgecolor="none",alpha=.85,pad=3))
    cax=fig.add_axes([.065,.119,.55,.014])
    cb=fig.colorbar(plt.cm.ScalarMappable(norm=NORM,cmap=CMAP),cax=cax,orientation="horizontal")
    cb.set_label("Dispersion coefficient a",fontsize=10);cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=9,length=3)
    fig.text(.68,.125,"Dots A–D identify the same waves in each projection.\nD is a long-period wave (L = 38.32).\nDashed profile lines show the period averages.",fontsize=9.3,color=MUTED,va="top",linespacing=1.5)
    fig.text(.055,.049,"L is the full repeat distance; A = max(v) − min(v); M is the average over one period. All units are those of the displayed normalized equation.",fontsize=9,color=MUTED)
    fig.text(.055,.025,"The pulse has infinite wavelength and is shown only as a limit. At a = 13/6 the profile becomes constant, so its wavelength is also a branch limit.",fontsize=9,color=MUTED)
    return fig


def pulse_figure(records,report):
    fig=plt.figure(figsize=(13.8,8.5))
    fig.text(.065,.944,"A periodic train opens into a solitary pulse",fontsize=24)
    fig.text(.065,.888,r"$a\to-13/6$:  $L\to\infty$,   $A\to35/12$,   $M\to2$",fontsize=16,color=MUTED)
    ax=fig.add_axes([.075,.27,.57,.53]);style(ax)
    z=np.linspace(-48,48,4801)
    for a,ls in [(-2.16,(0,(4,3))),(-2.16666,(0,(2,2)))]:
        r=state(a=a)
        ax.plot(z,profile(r,z),lw=1.5,ls=ls,label=f"a = {a:g},  L = {r['wavelength']:.2f}",alpha=.75)
    ax.plot(z,pulse(z),color=INK,lw=2.5,label="Exact solitary pulse")
    ax.axhline(2,color=MUTED,lw=.8,ls=(0,(4,4)))
    ax.set_xlim(-48,48);ax.set_ylim(-1.15,2.55)
    ax.set_xlabel("Travelling coordinate z");ax.set_ylabel("v(z)")
    ax.legend(loc="lower left",frameon=True,facecolor=PANEL,edgecolor=GRID,fontsize=9)
    ax2=fig.add_axes([.73,.52,.24,.28]);style(ax2)
    large=[r for r in records if r["a"]<0 and r["wavelength"]>20]
    L=np.array([r["wavelength"] for r in large]);M=np.array([r["mean"] for r in large])
    ax2.plot(1/L,M,color=CMAP(.02),lw=2.5,label="Periodic waves")
    inv=np.linspace(0,.05,150)
    ax2.plot(inv,2-PULSE_DEFICIT*inv,color=INK,ls=(0,(4,3)),lw=1,label="Pulse asymptotic")
    ax2.scatter([0],[2],marker="D",s=35,facecolor=PANEL,edgecolor=INK,zorder=4)
    ax2.set_xlim(-.002,.052);ax2.set_ylim(1,2.08)
    ax2.set_xlabel("Inverse wavelength 1/L",fontsize=10);ax2.set_ylabel("Mean M",fontsize=10)
    ax2.set_title("How the mean returns to 2",loc="left",fontsize=10.5,pad=12)
    ax2.legend(loc="lower left",fontsize=8,frameon=False)
    fig.text(.70,.355,r"$v_*(z)=2-\frac{35}{12}\,\operatorname{sech}^{4}\!\left(\frac{z}{2\sqrt{6}}\right)$",fontsize=16)
    fig.text(.70,.245,r"$M=2-\frac{70\sqrt{6}}{9L}+o(L^{-1})$",fontsize=15)
    fig.text(.065,.125,"The pulse keeps a finite depth while the separation between pulses diverges. Its mean approaches the background because a single depression",fontsize=10,color=MUTED)
    fig.text(.065,.098,"occupies a vanishing fraction of the period. The full ODE identity, mean formula, extrema, and this limiting profile are checked by the accompanying script.",fontsize=10,color=MUTED)
    return fig


def three_dimensional(records,examples):
    fig=plt.figure(figsize=(10,8))
    ax=fig.add_subplot(111,projection="3d")
    xyz=np.array([[r["wavelength"],r["amplitude"],r["mean"]] for r in records if r["wavelength"]<=80])
    colors=np.array([r["a"] for r in records if r["wavelength"]<=80])
    ax.add_collection3d(Line3DCollection(np.stack([xyz[:-1],xyz[1:]],axis=1),array=colors[:-1],cmap=CMAP,norm=NORM,linewidth=2.5))
    ax.set(xlim=(10,80),ylim=(0,3.05),zlim=(.4,2.06),xlabel="Wavelength L",ylabel="Amplitude A",zlabel="Mean M")
    ax.view_init(elev=22,azim=-57)
    for i,r in enumerate(examples):
        if r["wavelength"]<=80:
            ax.scatter(r["wavelength"],r["amplitude"],r["mean"],color=CMAP(NORM(r["a"])),s=38,edgecolor=PANEL)
            ax.text(r["wavelength"],r["amplitude"],r["mean"]+.04,chr(65+i),fontsize=10)
    ax.set_box_aspect((1.7,1,1));fig.suptitle("One branch in (wavelength, amplitude, mean)",fontsize=20,y=.94)
    fig.text(.10,.07,"Color follows a. The long-period tail continues beyond the window toward A = 35/12 and M = 2.",fontsize=10,color=MUTED)
    return fig


def interactive(records,examples):
    # Plain offline Plotly: rotate the three-dimensional curve and use the
    # selector to compare profiles, with all data embedded in this one file.
    fig=make_subplots(rows=1,cols=2,specs=[[{"type":"scene"},{"type":"xy"}]],column_widths=[.60,.40],horizontal_spacing=.09,
        subplot_titles=("Wavelength · amplitude · mean","A real wave on the same spatial scale"))
    data=records[::3]
    cs=[[i/10,matplotlib.colors.to_hex(CMAP(i/10))] for i in range(11)]
    custom=[[r["a_label"],r["distance_to_endpoint"]] for r in data]
    fig.add_trace(go.Scatter3d(x=[r["wavelength"] for r in data],y=[r["amplitude"] for r in data],z=[r["mean"] for r in data],
        mode="lines",line=dict(color=[r["a"] for r in data],colorscale=cs,cmin=-AC,cmax=AC,width=6,colorbar=dict(title="a",x=.55,len=.6)),
        customdata=custom,hovertemplate="a = %{customdata[0]}<br>L = %{x:.5g}<br>A = %{y:.5g}<br>M = %{z:.5g}<extra></extra>",name="Physical branch"),row=1,col=1)
    initial=examples[1]
    fig.add_trace(go.Scatter3d(x=[initial["wavelength"]],y=[initial["amplitude"]],z=[initial["mean"]],mode="markers",marker=dict(size=6,color=INK),name="Selected wave"),row=1,col=1)
    z=np.linspace(-30,30,1601)
    fig.add_trace(go.Scatter(x=z,y=profile(initial,z),mode="lines",line=dict(color=INK,width=2.6),name="v(z)"),row=1,col=2)
    fig.add_trace(go.Scatter(x=[-30,30],y=[initial["mean"]]*2,mode="lines",line=dict(color=MUTED,width=1,dash="dash"),name="Period mean"),row=1,col=2)
    choices=[state(a=a) for a in [2.1666,2.16,2.1,2,1,0,-1,-1.5,-2,-2.16,-2.16666]]
    choices.extend([state(side=-1,distance=1e-10),state(side=-1,distance=1e-30)])
    buttons=[]
    for r in choices:
        buttons.append(dict(label="a = "+r["a_label"],method="update",args=[{
            "x":[[r["wavelength"]],z,[-30,30]],"y":[[r["amplitude"]],profile(r,z),[r["mean"]]*2],"z":[[r["mean"]],None,None]},
            {"title.text":f"Kawahara physical branch · L = {r['wavelength']:.4f}, A = {r['amplitude']:.4f}, M = {r['mean']:.4f}"},[1,2,3]]))
    fig.update_layout(title=dict(text="Kawahara physical branch · select a wave and rotate the curve",x=.035,font=dict(size=22)),height=750,
        paper_bgcolor=PAPER,plot_bgcolor=PANEL,font=dict(family="Arial",color=INK,size=12),showlegend=False,
        margin=dict(t=130,b=75,l=25,r=35),updatemenus=[dict(buttons=buttons,active=5,x=.73,y=1.19,xanchor="left",yanchor="top")],
        scene=dict(xaxis=dict(title="Wavelength L",type="log",range=[math.log10(11),math.log10(300)],tickvals=[12,20,40,80,160]),
                   yaxis=dict(title="Amplitude A",range=[0,3.1]),zaxis=dict(title="Mean M",range=[.4,2.1]),
                   bgcolor=PANEL,aspectmode="manual",aspectratio=dict(x=1.4,y=1,z=1),camera=dict(eye=dict(x=1.6,y=-1.65,z=1.15))))
    fig.update_xaxes(title_text="Travelling coordinate z",range=[-30,30],row=1,col=2)
    fig.update_yaxes(title_text="v(z)",range=[-1.05,2.6],row=1,col=2)
    body=fig.to_html(full_html=False,include_plotlyjs=True,config={"responsive":True,"displaylogo":False},div_id="kawahara-physical")
    html='''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Kawahara in physical coordinates</title><style>body{margin:0;background:#f8f5ed;color:#203344;font:16px/1.55 system-ui,sans-serif}main{max-width:1450px;margin:auto;padding:22px}a{color:#177989}p{max-width:1000px}.note{font-size:14px;color:#687581}</style><main>
<p><a href="../index.html">Lattice atlas</a> · <a href="kawahara_physical_coordinates.png">Static figure</a> · <a href="kawahara_physical.pdf">PDF figures</a> · <a href="kawahara_physical.csv">CSV</a></p>
<p>Exact smooth real waves of <strong>v⁗ + a v″ − v + v²/2 = 0</strong>. L is the full repeat distance, A is the peak-to-trough amplitude, and M is the period average. Changing a changes the equation's dispersion coefficient.</p>
'''+body+'''<p class="note">The wavelength axis is logarithmic. At the small-amplitude endpoint, (L,A,M) tends to (2π√6,0,2). At the pulse endpoint, L tends to infinity, A tends to 35/12, and M tends to 2. Neither endpoint is a nonconstant periodic wave. The mean of an isolated pulse is quoted only as the limit of periodic averages.</p>
<p class="note">Near the pulse, many distinct finite periods require values of a much closer to −13/6 than double precision can distinguish. The CSV retains the distance to the endpoint and high-precision coefficient decimals. This is an existence plot; it does not encode stability.</p></main></html>'''
    (HERE/"kawahara_physical.html").write_text(html)


def main():
    records=build_data();report=validate(records)
    columns=list(records[0])
    with (HERE/"kawahara_physical.csv").open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=columns);writer.writeheader();writer.writerows(records)
    examples=[state(a=a) for a in [2.16,0,-2]]+[state(side=-1,distance=1e-10)]
    figures=[("kawahara_physical_coordinates",coordinate_figure(records,examples,report)),
             ("kawahara_pulse_limit",pulse_figure(records,report)),
             ("kawahara_physical_3d",three_dimensional(records,examples))]
    with PdfPages(HERE/"kawahara_physical.pdf",metadata={"Title":"Kawahara waves in wavelength, amplitude, and mean coordinates"}) as pdf:
        for name,fig in figures:
            for extension in ("png","svg"):fig.savefig(HERE/(name+"."+extension),dpi=170)
            pdf.savefig(fig);plt.close(fig)
    interactive(records,examples)
    print(f"PASS: {len(records)} states; exact ODE and pulse identities; seven independent quadrature checks.")
    print("Shortest sampled branch location:",report["shortest_wave"]["a"],report["shortest_wave"]["wavelength"])
    print("Saved physical-coordinate figures, pulse comparison, 3-D view, offline interactive HTML, CSV, and checks.")


if __name__=="__main__":main()
