"""Build interactive atlas data. Use --figures for legacy publication plates."""
import argparse
import csv
import json
import math
import logging
import re
import textwrap
import warnings
from pathlib import Path

import numpy as np
import sympy as s
import mpmath as mp
from scipy.optimize import linear_sum_assignment
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap,Normalize
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D

from models import HERE,build_models,serializable,t,g2

PAPER="#f8f5ed"; PANEL="#fffdf8"; INK="#203344"; MUTED="#687581"
GRID="#dfe4e3"; CUSP="#929ba3"; INSIDE="#edf3e8"; ORANGE="#d86849"
CMAP=LinearSegmentedColormap.from_list("wave_parameter",["#177989","#4b8f9c","#beaa80","#e59761","#ad3c63"])
plt.rcParams.update({"figure.facecolor":PAPER,"axes.facecolor":PANEL,"text.color":INK,
    "axes.labelcolor":INK,"xtick.color":MUTED,"ytick.color":MUTED,"axes.edgecolor":GRID,
    "font.family":"DejaVu Sans","font.size":10,"axes.titleweight":"normal",
    "savefig.facecolor":PAPER,"svg.fonttype":"none","pdf.fonttype":42})
logging.getLogger("fontTools.ttLib.tables._h_e_a_d").setLevel(logging.ERROR)


def real_scalar(value):
    value=complex(value)
    if not np.isfinite(value) or abs(value.imag)>1e-11+abs(value.real)*1e-7:return None
    return float(value.real)


def real_polynomial_roots(expression,variable=t):
    poly=s.Poly(s.sympify(expression),variable)
    return [float((bounds[0]+bounds[1])/2) for bounds,multiplicity in poly.intervals(eps=s.Rational(1,10**22))]


def is_degenerate(point):
    scale=max(abs(point["g2"])**3,27*point["g3"]**2,1e-35)
    return abs(point["delta"])<=2e-6*scale


def boundary_data(model):
    path=HERE/"data"/"boundaries.json"
    return json.loads(path.read_text()).get(model["slug"],{}) if path.exists() else {}


def special_parameters(model):
    special=[float(x) for x in model.get("specials",[])]
    lo,hi=model["interval"]
    for factor in boundary_data(model).get("cusp_factors",[]):
        special.extend(real_polynomial_roots(factor))
    for br in model["branches"]:
        if "root_variable" in br or model.get("kawahara"):continue
        delta=s.factor(br["g2"]**3-27*br["g3"]**2)
        numerator=s.fraction(delta)[0]
        if not numerator.has(t):continue
        for factor,_ in s.factor_list(numerator,t)[1]:
            poly=s.Poly(factor,t)
            for r in np.roots([float(c) for c in poly.all_coeffs()]):
                r=real_scalar(r)
                if r is not None and lo<=r<=hi:special.append(r)
    return sorted({round(x,13) for x in special if lo<=x<=hi})


def sample(model,count=1801):
    lo,hi=model["interval"]
    specials=special_parameters(model)
    boundaries=boundary_data(model)
    folds=[v for f in boundaries.get("fold_factors",[]) for v in real_polynomial_roots(f) if lo<=v<=hi]
    breaks=model.get("breaks",[])+[v for f in boundaries.get("leading_factors",[]) for v in real_polynomial_roots(f)]
    parameters=np.linspace(lo,hi,count)
    if "zoom_parameter" in model:
        parameters=np.r_[parameters,np.linspace(*model["zoom_parameter"],1001)]
    parameters=np.unique(np.r_[parameters,specials,folds])
    # Symbolic and numerical event locators can name the same parameter a few
    # ulps apart. Keep one frame, so the exceptional-fibre replacement is unique.
    parameters=parameters[np.r_[True,np.diff(parameters)>1e-11]]
    for brk in model.get("breaks",[]):parameters=parameters[abs(parameters-brk)>1e-7]
    tracks=[]
    frames=[dict(t=float(T),points=[]) for T in parameters]
    for br_index,br in enumerate(model["branches"]):
        variable=br.get("root_variable")
        args=(t,variable) if variable is not None else (t,)
        functions={name:s.lambdify(args,br[name],"numpy",cse=True) for name in ("g2","g3","constant")}
        operators={j:s.lambdify(args,val,"numpy",cse=True) for j,val in br["operator"].items()}
        if variable is not None:
            coefficients=s.lambdify(t,s.Poly(br["root_polynomial"],variable).all_coeffs(),"numpy",cse=True)
            accurate=model["slug"]=="n2_p5_quadratic_fifth"
            if accurate:
                mpcoeff=s.lambdify(t,s.Poly(br["root_polynomial"],variable).all_coeffs(),"mpmath",cse=True)
                mpfunctions={name:s.lambdify(args,br[name],"mpmath",cse=True) for name in ("g2","g3","constant")}
                mpoperators={j:s.lambdify(args,val,"mpmath",cse=True) for j,val in br["operator"].items()}
            degree=s.degree(br["root_polynomial"],variable)
            previous=None
            roots_grid=[]
            for T in parameters:
                cf=np.asarray(coefficients(T),dtype=complex)
                roots=np.roots(cf)
                if any(abs(T-event)<1e-12 for event in folds):
                    roots=np.array([complex(z.real) if abs(z.imag)<3e-6*(1+abs(z.real)) else z for z in roots])
                if len(roots)!=degree:
                    roots=np.r_[roots,np.full(int(degree)-len(roots),complex(np.nan))]
                if previous is not None and np.isfinite(roots).all() and np.isfinite(previous).all():
                    _,perm=linear_sum_assignment(abs(previous[:,None]-roots[None,:]))
                    roots=roots[perm]
                else:roots=roots[np.lexsort((roots.imag,roots.real))]
                roots_grid.append(roots);previous=roots
            roots_grid=np.array(roots_grid)
        else:
            degree=1;roots_grid=np.zeros((len(parameters),1));accurate=False
        for root_index in range(int(degree)):
            xy=np.full((len(parameters),2),np.nan)
            meta=[]
            for k,T in enumerate(parameters):
                root=roots_grid[k,root_index]
                if variable is not None and real_scalar(root) is None:continue
                values=(complex(T),root) if variable is not None else (complex(T),)
                high_values=None
                if accurate and not any(abs(T-event)<1e-12 for event in folds):
                    with mp.workdps(55):
                        tt=mp.mpf(float(T));rr=mp.mpf(float(root.real));cf=mpcoeff(tt)
                        derivative=[c*(len(cf)-j-1) for j,c in enumerate(cf[:-1])]
                        for _ in range(12):
                            denominator=mp.polyval(derivative,rr)
                            if not denominator:break
                            correction=mp.polyval(cf,rr)/denominator;rr-=correction
                            if abs(correction)<mp.mpf('1e-48'):break
                        high_values=(tt,rr)
                        G2=float(mpfunctions["g2"](*high_values));G3=float(mpfunctions["g3"](*high_values))
                        coeff={str(j):float(f(*high_values)) for j,f in mpoperators.items()}
                        C=float(mpfunctions["constant"](*high_values))
                with np.errstate(all="ignore"):
                    if high_values is None:
                        G2=real_scalar(functions["g2"](*values));G3=real_scalar(functions["g3"](*values))
                if model.get("kawahara") and abs(abs(T)-13/6)<1e-10:
                    G2=1/432;G3=np.sign(T)/46656
                if G2 is None or G3 is None:continue
                xy[k]=G2,G3
                if high_values is None:
                    coeff={str(j):real_scalar(f(*values)) for j,f in operators.items()}
                    C=real_scalar(functions["constant"](*values))
                point=dict(g2=G2,g3=G3,delta=G2**3-27*G3**2,operator=coeff,constant=C,
                           branch=len(tracks),character=br["character"])
                if br.get("physical_profile"):point["physical"]=True
                if variable is not None:point["root"]=float(root.real)
                frames[k]["points"].append(point)
            tracks.append(dict(xy=xy,parameter=parameters,character=br["character"],label=br.get("label",f"branch {len(tracks)+1}"),
                               roots=roots_grid[:,root_index],breaks=breaks))
    for fibre in boundaries.get("exceptional_fibres",[]):
        T=float(s.sympify(fibre["t"]));R0=float(s.sympify(fibre["a2"]))
        k=int(np.argmin(abs(parameters-T)))
        if abs(parameters[k]-T)>1e-10:continue
        indices=[i for i,tr in enumerate(tracks) if abs(tr["roots"][k]-R0)<1e-5]
        roots=real_polynomial_roots(fibre["g2_polynomial"],g2)
        pts=[]
        f3=s.lambdify(g2,s.sympify(fibre["g3"]),"numpy")
        fc=s.lambdify(g2,s.sympify(fibre["C"]),"numpy")
        for G2 in roots:
            G3=float(f3(G2))
            pts.append(dict(g2=G2,g3=G3,delta=G2**3-27*G3**2,root=R0,character=1,
                operator={"4":1,"3":T,"2":R0,"1":float(s.sympify(fibre["a1"]))},constant=float(fc(G2)),exceptional=True))
        frames[k]["points"]=[pt for pt in frames[k]["points"] if pt["branch"] not in indices]
        if indices and pts:
            old=np.array([tracks[i]["xy"][max(0,k-1)] for i in indices])
            new=np.array([[pt["g2"],pt["g3"]] for pt in pts])
            costs=np.nan_to_num(np.linalg.norm(old[:,None,:]-new[None,:,:],axis=2),nan=1e20)
            aa,bb=linear_sum_assignment(costs)
            for ia,ib in zip(aa,bb):
                idx=indices[ia];pt=pts[ib];pt["branch"]=idx
                tracks[idx]["xy"][k]=pt["g2"],pt["g3"]
                frames[k]["points"].append(pt)
            for ia,idx in enumerate(indices):
                if ia not in aa:
                    ib=int(np.argmin(costs[ia]))
                    tracks[idx]["xy"][k]=pts[ib]["g2"],pts[ib]["g3"]
    # Numerical coordinates of distinct waves can agree at a coalescence.
    # Deduplication is only used for explicitly marked coalescences, never to
    # infer that two arbitrary profiles on the same lattice are one wave.
    counts=[]
    for frame in frames:
        pts=frame["points"]
        if (model.get("kawahara") or model["slug"]=="n3_p4_cubic_positive") and any(abs(frame["t"]-v)<1e-10 for v in specials):
            pts=pts[:1];frame["points"]=pts
        if any(abs(frame["t"]-event)<1e-12 for event in folds):
            unique=[]
            for pt in pts:
                if not any(abs(pt.get("root",1e10)-q.get("root",-1e10))<1e-5 and
                    abs(pt["g2"]-q["g2"])<2e-9 and abs(pt["g3"]-q["g3"])<2e-10 for q in unique):unique.append(pt)
            pts=unique;frame["points"]=pts
        counts.append(len(pts))
    specials=[T for T in specials if any(is_degenerate(pt) for pt in min(frames,key=lambda f:abs(f["t"]-T))["points"])]
    return dict(model=model,tracks=tracks,frames=frames,specials=specials,counts=np.array(counts),parameters=parameters)


def view_limits(data,zoom=False):
    model=data["model"]
    if zoom and "zoom_limits" in model:return np.array(model["zoom_limits"],dtype=float)
    if not zoom and "view_limits" in model:return np.array(model["view_limits"],dtype=float)
    if zoom and "zoom_parameter" in model:
        a,b=model["zoom_parameter"]
        points=np.vstack([tr["xy"][(tr["parameter"]>=a)&(tr["parameter"]<=b)] for tr in data["tracks"]])
    else:points=np.vstack([tr["xy"] for tr in data["tracks"]])
    points=points[np.isfinite(points).all(axis=1)]
    quant=model.get("view_quantiles",(0,1))
    if model["slug"]=="n2_p5_quadratic_fifth":quant=(.035,.965)
    limits=np.quantile(points,quant,axis=0).T
    if zoom and "zoom_parameter" not in model:
        special=[]
        for T in data["specials"]:
            frame=min(data["frames"],key=lambda f:abs(f["t"]-T))
            special.extend((p["g2"],p["g3"]) for p in frame["points"])
        if special:
            center=np.mean(special,axis=0)
        else:center=np.median(points,axis=0)
        span=np.maximum(limits[:,1]-limits[:,0],1e-12)*np.array([.32,.16])
        limits=np.column_stack((center-span/2,center+span/2))
    for j in range(2):
        low,high=limits[j];span=high-low
        if span<max(abs(high),abs(low),1e-12)*1e-8:
            span=max(abs(high)*.75,.02 if j==0 else .002)
            low-=span/2;high+=span/2
        else:
            low-=span*.07;high+=span*.07
        limits[j]=low,high
    return limits


def powers(limits):
    result=[]
    for low,high in limits:
        maximum=max(abs(low),abs(high))
        exponent=int(math.floor(math.log10(maximum)))
        result.append(exponent if abs(exponent)>=2 else 0)
    return result


def plane(ax,data,limits=None,small=False,zoom=False):
    limits=view_limits(data,zoom=zoom) if limits is None else np.array(limits)
    exponents=powers(limits);scale=10.**(-np.array(exponents))
    ax.set_xlim(limits[0]*scale[0]);ax.set_ylim(limits[1]*scale[1])
    ax.grid(color=GRID,lw=.6,alpha=.7,zorder=0)
    ax.set_axisbelow(True)
    for side in ("top","right"):ax.spines[side].set_visible(False)
    if limits[0,1]>0:
        x=np.linspace(max(limits[0,0],0),limits[0,1],800)
        y=np.sqrt(x**3/27)
        ax.fill_between(x*scale[0],-y*scale[1],y*scale[1],color=INSIDE,zorder=0)
        for sign in (-1,1):ax.plot(x*scale[0],sign*y*scale[1],color=CUSP,lw=1.4 if not small else 1,zorder=2)
    model=data["model"];norm=Normalize(*model["interval"])
    for tr in data["tracks"]:
        xy=tr["xy"]*scale
        segments=np.stack((xy[:-1],xy[1:]),axis=1)
        good=np.isfinite(segments).all(axis=(1,2))
        window=(limits[:,1]-limits[:,0])*scale
        good &= np.max(abs(segments[:,1]-segments[:,0])/window,axis=1)<.35
        # Break across the missing normalization chart; do not draw an
        # artificial chord between opposite sides of an asymptote.
        delta=np.diff(tr["parameter"])
        for brk in tr["breaks"]:
            good&=~((tr["parameter"][:-1]<brk)&(tr["parameter"][1:]>brk))
        width=2.5 if not small else 1.6
        lc=LineCollection(segments[good],array=tr["parameter"][:-1][good],cmap=CMAP,norm=norm,linewidth=width,zorder=4)
        if len(model["branches"])>1 and tr["character"]>1:lc.set_linestyle((0,(3,1.1)))
        ax.add_collection(lc)
    for T in data["specials"]:
        frame=min(data["frames"],key=lambda f:abs(f["t"]-T))
        for point in frame["points"]:
            if not is_degenerate(point):continue
            ax.plot(point["g2"]*scale[0],point["g3"]*scale[1],"D",color=ORANGE,mec=PANEL,mew=1.2,
                    ms=6.5 if not small else 4.2,zorder=6)
    fs=10 if not small else 8
    for j,name in enumerate(("g_2","g_3")):
        label=rf"${name}$"+(rf"  ($10^{{{exponents[j]}}}$)" if exponents[j] else "")
        (ax.set_xlabel if j==0 else ax.set_ylabel)(label,fontsize=fs)
    ax.tick_params(labelsize=fs-1)
    ax.locator_params(axis="both",nbins=5 if not small else 4)
    return limits


def plate(data):
    model=data["model"]
    fig=plt.figure(figsize=(13.8,8.5))
    fig.text(.065,.95,model["title"],fontsize=22,weight="medium",va="top")
    fig.text(.945,.951,rf"$(n,p)=({model['n']},{model['p']})$"+f"  ·  q = {model['q']}",ha="right",va="top",fontsize=12,color=MUTED)
    fig.text(.065,.876,"$"+model["equation"]+"$",fontsize=16)
    ax=fig.add_axes([.075,.24,.55,.575]);plane(ax,data)
    ax.text(.02,.965,"REAL LATTICE LOCUS",transform=ax.transAxes,fontsize=8.5,color=MUTED,va="top")
    zoom=fig.add_axes([.725,.51,.235,.305]);plane(zoom,data,small=True,zoom=True)
    zoom.set_title("Detail near the discriminant",loc="left",fontsize=10,color=MUTED,pad=10)
    countax=fig.add_axes([.725,.265,.235,.145])
    counts=data["counts"]
    countax.plot(data["parameters"],counts,color=INK,lw=1.6,drawstyle="steps-mid")
    for T in data["specials"]:countax.axvline(T,color=ORANGE,lw=.8,alpha=.55)
    countax.set_ylim(-.15,max(1,int(counts.max()))+.7)
    countax.set_yticks(range(max(1,int(counts.max()))+1))
    countax.set_xlim(*model["interval"])
    countax.set_xlabel(model["parameter"],fontsize=8)
    countax.tick_params(labelsize=8)
    countax.grid(color=GRID,lw=.6)
    for side in ("top","right"):countax.spines[side].set_visible(False)
    count_title="Real lattice solutions" if model["count_kind"]=="waves" else "Real compatible coefficient choices"
    if model.get("fixed_equation"):count_title="Selected energy: one lattice; energy is free"
    countax.set_title(count_title,loc="left",fontsize=9,pad=9)
    cax=fig.add_axes([.075,.151,.55,.017])
    cb=fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(*model["interval"]),cmap=CMAP),cax=cax,orientation="horizontal")
    cb.outline.set_visible(False);cb.ax.tick_params(labelsize=8,length=3)
    cb.set_label(model["parameter"],fontsize=9,labelpad=3)
    for j,note in enumerate(model["notes"]):fig.text(.065,.062-j*.025,note,fontsize=9.2,color=MUTED)
    fig.text(.725,.135,"Grey cusp: Δ = 0\nPale region: Δ > 0\nDiamonds: real degeneration parameters",
             fontsize=8.5,color=MUTED,linespacing=1.55)
    fig.text(.065,.012,"Lattice of the full twisted period group Γ  ·  Real invariants do not alone imply a bounded real profile.",fontsize=8,color=MUTED)
    return fig


def poster(collection,primary=True):
    selected=[d for d in collection if d["model"]["primary"]==primary]
    cols=4 if primary else 3;rows=4 if primary else 2
    fig,axes=plt.subplots(rows,cols,figsize=(19,19 if primary else 11))
    axes=np.array(axes).ravel()
    for ax,data in zip(axes,selected):
        m=data["model"];plane(ax,data,small=True)
        title=m["title"].split(" · ")[0]
        ax.set_title(f"({m['n']}, {m['p']})   {title}",loc="left",fontsize=10.5,pad=10)
    if primary:
        ax=axes[len(selected)];ax.axis("off")
        table=[["n \\ p"]+[str(p) for p in range(2,7)]]
        for n in range(2,8):table.append([str(n)]+[(str(p//(n-1)) if p%(n-1)==0 else "—") for p in range(2,7)])
        tab=ax.table(cellText=table,loc="center",cellLoc="center",bbox=[0,.15,1,.75])
        for (r,c),cell in tab.get_celld().items():
            cell.set_edgecolor(PAPER);cell.set_linewidth(3)
            cell.set_facecolor(INSIDE if r and c and cell.get_text().get_text()!="—" else "#eaece8")
            cell.get_text().set_color(INK);cell.get_text().set_fontsize(11)
        ax.set_title("All 30 low-order pairs",loc="left",fontsize=11)
        ax.text(0,.04,"Numbers are pole orders q = p/(n − 1).\nDashes: no integral pole balance.",fontsize=10,color=MUTED)
        ax=axes[len(selected)+1];ax.axis("off")
        ax.text(0,.95,"Reading the plane",fontsize=15,va="top",color=INK)
        ax.text(0,.81,"Δ = g₂³ − 27g₃²\n\nGrey cusp: degenerate lattice\nPale interior: Δ > 0\nOff the cusp: elliptic lattice\nOrigin: rational degeneration\n\nColors follow the parameter\nwithin each declared coefficient slice.",fontsize=11,color=MUTED,va="top",linespacing=1.6)
        ax=axes[-1];ax.axis("off")
        ax.text(0,.95,"13 pairs · 19 slices",fontsize=15,va="top")
        ax.text(0,.81,"The full book contains every plate,\ncompanion slices, and enlarged details.\n\nOpen index.html for the slider gallery.\n\nCounts are modulo translation and phase.\nA real lattice need not give a smooth,\nbounded real-valued wave.\n\nAt p = 2 the conservative equation\ncan carry a continuum of waves.",fontsize=10.5,color=MUTED,va="top",linespacing=1.6)
    fig.suptitle("Travelling waves in the plane of lattices" if primary else "Companion slices · ordinary and twisted branches",x=.04,y=.985,ha="left",fontsize=26)
    fig.text(.04,.948,"All pole-compatible (n, p) with 2 ≤ p ≤ 6 and 2 ≤ n ≤ 7" if primary else "Additional normalized coefficient families for the same low-order pairs",fontsize=12,color=MUTED)
    fig.subplots_adjust(left=.052,right=.968,top=.91,bottom=.045,hspace=.48,wspace=.38)
    return fig


def write_csv(data):
    path=HERE/"data"/(data["model"]["slug"]+".csv")
    with path.open("w",newline="") as f:
        writer=csv.writer(f);writer.writerow(["parameter","branch","g2","g3","discriminant","character","operator_ascending","constant"])
        p=data["model"]["p"]
        for frame in data["frames"]:
            for pt in frame["points"]:
                op=[pt["operator"].get(str(j),0) for j in range(p)]+[1]
                writer.writerow([frame["t"],pt["branch"],pt["g2"],pt["g3"],pt["delta"],pt["character"],json.dumps(op),pt["constant"]])


def gallery_data(collection):
    out=[]
    for data in collection:
        m=data["model"]
        indices=np.unique(np.r_[np.linspace(0,len(data["frames"])-1,701,dtype=int),
            [int(np.argmin(abs(data["parameters"]-T))) for T in data["specials"]]]).astype(int)
        limits=view_limits(data)
        record={k:m[k] for k in ("slug","n","p","q","title","parameter","interval","notes","count_kind","primary")}
        record.update(limits=limits.tolist(),fixed_equation=m.get("fixed_equation",False),
                      K=int((-1)**m["p"]*s.rf(m["q"],m["p"])),specials=data["specials"],
                      frames=[data["frames"][i] for i in indices],tracks=[])
        for tr in data["tracks"]:
            points=[]
            for i in indices:
                xy=tr["xy"][i]
                points.append([float(v) for v in xy]+[float(tr["parameter"][i])] if np.isfinite(xy).all() else None)
            record["tracks"].append(dict(points=points,character=tr["character"],breaks=tr["breaks"]))
        out.append(record)
    return out


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--figures',action='store_true',help='Also regenerate the archived static plates and PDF book.')
    args=parser.parse_args()
    models=build_models()
    (HERE/"figures").mkdir(exist_ok=True)
    (HERE/"data"/"models.json").write_text(json.dumps(serializable(models),ensure_ascii=False,indent=2)+"\n")
    collection=[]
    for m in models:
        data=sample(m);collection.append(data);write_csv(data)
        print("SAMPLED",m["slug"],len(data["frames"]),"parameter values",flush=True)
    if args.figures:
        with PdfPages(HERE/"g2g3_atlas.pdf",metadata={"Title":"Travelling waves in the plane of lattices","Author":"Computational companion to the meromorphic travelling-wave papers"}) as pdf:
            for primary,slug in [(True,"atlas_overview"),(False,"companion_overview")]:
                fig=poster(collection,primary)
                fig.savefig(HERE/"figures"/(slug+".png"),dpi=145)
                fig.savefig(HERE/"figures"/(slug+".pdf"))
                pdf.savefig(fig);plt.close(fig)
            for data in collection:
                fig=plate(data);slug=data["model"]["slug"]
                for extension in ("png","svg","pdf"):
                    fig.savefig(HERE/"figures"/(slug+"."+extension),dpi=170)
                pdf.savefig(fig);plt.close(fig)
                print("RENDERED",slug,flush=True)
            from c4 import figure as c4_figure
            fig=c4_figure();pdf.savefig(fig);plt.close(fig)
    payload=gallery_data(collection)
    from profiles import attach_profiles
    attach_profiles(payload)
    (HERE/"data"/"atlas.json").write_text(json.dumps(payload,separators=(",",":"),ensure_ascii=False,allow_nan=False)+"\n")
    # The template needs no server, CDN, package installation, or network.
    from journeys import compile_site,build_journey
    compile_site(payload,build_journey())
    (HERE/"data"/"coverage.json").write_text(json.dumps({"pairs":[[m["n"],m["p"]] for m in models if m["primary"]],
        "excluded_pairs":[[n,p] for p in range(2,7) for n in range(2,8) if p%(n-1)],"slices":len(models),
        "sample_counts":{d["model"]["slug"]:sum(len(f["points"]) for f in d["frames"]) for d in collection}},indent=2)+"\n")
    print("WROTE interactive index.html, atlas data, CSVs and validation fixtures",flush=True)


if __name__=="__main__":main()
