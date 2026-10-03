"""Compile every selected lattice into a small, browser-evaluable profile.

Only polynomial coefficients and elliptic geometry are shipped, not sampled
curves. The browser evaluates Jacobi functions and continues the characters
through zeros. High precision keeps nearly singular lattices distinguishable.
"""
import json
import re
from functools import lru_cache
from pathlib import Path
import mpmath as mp
import sympy as s
from models import HERE,build_models,t,X,Y,g2,g3,read
from verify import power_expression

R=s.Symbol("profile_root")
ARGS=(t,R,g2,g3)


@lru_cache(None)
def evaluator(expression):
    return s.lambdify(ARGS,s.sympify(expression),"mpmath",cse=True)


def value(expression,values):
    v=evaluator(str(expression))(*values)
    if isinstance(v,mp.mpc) and abs(v.imag)<mp.mpf('1e-25'):v=v.real
    return v


def geometry(G2,G3,node=False):
    D=G2**3-27*G3**2
    if abs(G2)+abs(G3)<mp.mpf('1e-40'):
        return dict(type="rational",delta=0,period=None,scale=1)
    if node or abs(D)<mp.mpf('1e-43')*(abs(G2)**3+27*abs(G3)**2):
        e=mp.sign(-G3)*mp.sqrt(max(G2,0)/12)
        if e>0:
            return dict(type="rect",delta=0,node=True,e=[float(e),float(e),float(-2*e)],G=float(3*e),D=float(3*e),m=1,mc=0,K=None,period=None)
        G=-3*e
        return dict(type="rect",delta=0,node=True,e=[float(-2*e),float(e),float(e)],G=float(G),D=0,m=0,mc=1,K=float(mp.pi/2),period=float(mp.pi/mp.sqrt(G)))
    if D>0:
        radius=mp.sqrt(G2/12);theta=mp.acos(max(-1,min(1,3*mp.sqrt(3)*G3/G2**mp.mpf('1.5'))))/3
        roots=[2*radius*mp.cos(theta-j*2*mp.pi/3) for j in (0,1,-1)]
        G=roots[0]-roots[2];gap=roots[1]-roots[2];m=gap/G;K=mp.ellipk(m)
        return dict(type="rect",delta=1,e=list(map(float,roots)),G=float(G),D=float(gap),m=float(m),mc=float(1-m),K=float(K),period=float(2*K/mp.sqrt(G)))
    root=lambda x:mp.sign(x)*abs(x)**(mp.mpf(1)/3)
    rad=mp.sqrt(G3**2/64-G2**3/1728)
    e=root(G3/8+rad)+root(G3/8-rad)
    H=mp.sqrt(3*e**2-G2/4);m=mp.mpf('.5')-3*e/(4*H);K=mp.ellipk(m)
    return dict(type="one",delta=-1,e=float(e),H=float(H),m=float(m),mc=float(1-m),K=float(K),period=float(2*K/mp.sqrt(H)))


def base_at(g,x,regular):
    if g["type"]=="rational":return 1/x**2,-2/x**3
    m=mp.mpf(g["m"]);u=mp.sqrt(g.get("G",g.get("H")))*x
    if g["type"]=="one":u*=2
    if m==1:sn=mp.tanh(u);cn=dn=1/mp.cosh(u)
    elif m==0:sn=mp.sin(u);cn=mp.cos(u);dn=mp.mpf(1)
    else:sn=mp.ellipfun("sn",u,m);cn=mp.ellipfun("cn",u,m);dn=mp.ellipfun("dn",u,m)
    if g["type"]=="one":
        H=mp.mpf(g["H"])
        if abs(1-cn)<mp.mpf('1e-35'):raise ZeroDivisionError
        return g["e"]+H*(1+cn)/(1-cn),-4*H**mp.mpf('1.5')*sn*dn/(1-cn)**2
    G=mp.mpf(g["G"]);gap=mp.mpf(g["D"])
    if regular:return g["e"][2]+gap*sn**2,2*gap*mp.sqrt(G)*sn*cn*dn
    if abs(sn)<mp.mpf('1e-35'):raise ZeroDivisionError
    return g["e"][2]+G/sn**2,-2*G**mp.mpf('1.5')*cn*dn/sn**3


def polynomial(expression,values):
    return [[int(i),int(j),float(value(c,values))] for (i,j),c in s.Poly(s.expand(expression),X,Y).terms()]


def poly_at(terms,x,y):return sum(mp.mpf(c)*x**i*y**j for i,j,c in terms)


def zero_data(g,A,B,order,regular):
    if regular and (g["type"]!="rect" or g.get("D")==0):return None
    if g["type"]=="rational":return None
    L=g["period"]
    if L is None:
        G=mp.mpf(g["G"]);e3=mp.mpf(g["e"][2]);r=(A-e3)/G
        if regular:
            if r>1-mp.mpf('1e-11') or r<-mp.mpf('1e-11'):return None
            if abs(r)<mp.mpf('1e-11'):return dict(position=0)
            if B==0:return None
            h=mp.sign(B)*mp.sqrt(r)
        else:
            if r<1+mp.mpf('1e-11') or B==0:return None
            h=-mp.sign(B)/mp.sqrt(r)
        return dict(position=float(mp.atanh(h)/mp.sqrt(G)))
    best=None
    normx=mp.mpf(g["D"] if regular else g.get("G",g.get("H")))
    normy=normx*mp.sqrt(g.get("G",g.get("H")))
    for j in range(order):
        try:xx,yy=base_at(g,mp.mpf(L)*j/order,regular)
        except ZeroDivisionError:continue
        error=abs((xx-A)/normx)+abs((yy-B)/normy)
        if best is None or error<best[0]:best=(error,j/order)
    if best and best[0]<mp.mpf('2e-5'):return dict(fraction=best[1])
    return None


def stable_fifth(T,point):
    # Undivided principal-part readout, including the three exceptional fibres.
    b=-s.Rational(1,14);c=-(49*t-4)/12936
    r=R;d=(-539*r+91*t-20)/530376
    a1=s.Symbol("a1")
    h=(a1*b+6*r*c-168*b*g2+15120*c*d+6*d*t+42*g2)/7560
    return (-X*Y/2+b*X**2+c*Y+d*X+h).subs(a1,s.Float(point["operator"]["1"],17))


def attach_profiles(payload):
    models={m["slug"]:m for m in build_models()};tate=read("tate")
    fixtures=[];total=0
    with mp.workdps(55):
        for record in payload:
            model=models[record["slug"]]
            record["display_scale"]=-1680 if model.get("kawahara") else 1
            record["variable"]="v" if model.get("kawahara") or model["slug"]=="n2_p3_ks" else "u"
            fixture_indices=set([0,len(record["frames"])//4,len(record["frames"])//2,3*len(record["frames"])//4,len(record["frames"])-1])
            fixture_indices.update(min(range(len(record["frames"])),key=lambda j:abs(record["frames"][j]["t"]-c)) for c in record["specials"])
            for index,frame in enumerate(record["frames"]):
                T=mp.mpf(str(frame["t"]))
                for point in frame["points"]:
                    br=model["branches"][point["branch"] if len(model["branches"])>1 else 0]
                    vals=(T,mp.mpf(str(point.get("root",0))),mp.mpf(str(point["g2"])),mp.mpf(str(point["g3"])))
                    near_event=any(abs(frame["t"]-c)<1e-10 for c in record["specials"])
                    ratio=abs(point["delta"])/max(abs(point["g2"])**3+27*point["g3"]**2,1e-40)
                    node=near_event and ratio<2e-5
                    if model.get("kawahara") and node:
                        vals=(mp.sign(T)*mp.mpf(13)/6,vals[1],mp.mpf(1)/432,mp.sign(T)/46656)
                    elif "root_variable" not in br and not node:
                        G2=value(br["g2"],vals);G3=value(br["g3"],vals)
                        vals=(T,vals[1],G2,G3)
                    g=geometry(vals[2],vals[3],node)
                    spec=dict(scale=record["display_scale"],variable=record["variable"],geometry=g)
                    substitute={br["root_variable"]:R} if "root_variable" in br and br["root_variable"]!=g2 else {}
                    expression=br["profile"]
                    if model["slug"]=="n2_p5_quadratic_fifth":
                        expression=stable_fifth(T,point)
                    if not isinstance(expression,str) and not br.get("profile_kind"):
                        expression=s.sympify(expression).subs(substitute)
                        spec.update(kind="polynomial",terms=polynomial(expression,vals))
                        power=s.expand(expression**(model["n"]-1))
                    elif br.get("physical_profile"):
                        expression=s.sympify(expression)
                        spec.update(kind="polynomial",terms=polynomial(expression,vals))
                        power=expression
                    elif br.get("profile_kind"):
                        power=power_expression(model,br).subs(substitute)
                        terms=polynomial(power,vals)
                        if br["profile_kind"]=="tate_root":
                            origin=tate[f"tate_{model['n']-1}"]
                            lam=value(br["rate"],vals)
                            A=value(s.sympify(origin["point"][0]),vals)*lam**2
                            B=value(s.sympify(origin["point"][1]),vals)*lam**3
                        else:A=(15-40*vals[0])/1200;B=(40*vals[0]-11)/2000
                        order=model["n"]-1
                        zeros={};signs={}
                        for mode in ("regular","axis"):
                            regular=mode=="regular" and g["type"]=="rect"
                            zeros[mode]=zero_data(g,A,B,order,regular)
                            samples=[]
                            length=g["period"] or 8/mp.sqrt(g.get("G",1))
                            for f in (.13,.27,.41):
                                try:xx,yy=base_at(g,length*f,regular);samples.append(poly_at(terms,xx,yy))
                                except ZeroDivisionError:pass
                            test=max(samples,key=abs) if samples else 1
                            signs[mode]=-1 if test<0 else 1
                        spec.update(kind="root",order=order,terms=terms,point=[float(A),float(B)],zeros=zeros,power_sign=signs)
                    else:
                        if expression.startswith("sqrt("):
                            E=-s.expand(s.sympify(expression[5:-1])).subs(X,0)
                            form="generator";B=s.Integer(0);DD=s.Integer(0)
                        else:
                            E=br["half_root"];B=br.get("profile_b",s.Integer(0));DD=br.get("profile_d",s.Integer(0))
                            form="derivative" if model["p"]==4 else "multiply"
                        E=value(E,vals);B=value(B,vals);DD=value(DD,vals)
                        root_index=min(range(3),key=lambda j:abs(E-g["e"][j])) if g["type"]=="rect" else 0
                        spec.update(kind="halfroot",root_index=root_index,root=float(E),form=form,b=float(B),d=float(DD))
                        power=power_expression(model,br).subs(substitute)
                    point["profile"]=spec;total+=1
                    if index in fixture_indices:
                        fixtures.append(dict(slug=model["slug"],n=model["n"],p=model["p"],K=record["K"],t=frame["t"],
                            g2=point["g2"],g3=point["g3"],operator=point["operator"],constant=point["constant"],physical=br.get("physical_profile",False),
                            profile=spec,power=polynomial(power,vals)))
            print("PROFILES",record["slug"],flush=True)
    (HERE/"data"/"profile_fixtures.json").write_text(json.dumps(fixtures,separators=(",",":"),allow_nan=False)+"\n")
    return total


def rebuild():
    # Reuse the existing numerical atlas when only the web application changes.
    datafile=HERE/"data"/"atlas.json"
    if datafile.exists():payload=json.loads(datafile.read_text())
    else:
        html=(HERE/"index.html").read_text();start=re.search(r"\bconst\s+ATLAS\s*=\s*",html).end()
        payload,_=json.JSONDecoder().raw_decode(html[start:])
    total=attach_profiles(payload)
    datafile.write_text(json.dumps(payload,separators=(",",":"),ensure_ascii=False,allow_nan=False)+"\n")
    from journeys import compile_site
    compile_site(payload)
    print(f"Compiled {total} selected-wave profiles",flush=True)


if __name__=="__main__":rebuild()
