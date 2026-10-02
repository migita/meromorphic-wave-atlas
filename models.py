"""Declared coefficient slices, exact formulae, and reproducible plot domains.

The atlas covers every integral pole-balance pair in 2<=p<=6, 2<=n<=7.
It is a collection of stated slices, not a claim that a one-dimensional
picture exhausts every coefficient of a high-dimensional operator space.
"""
from pathlib import Path
import json
import sympy as s

HERE=Path(__file__).resolve().parent
t,X,Y,g2,g3,e,d,a2=s.symbols("t X Y g2 g3 e d a2")
R=s.Rational


def read(name):
    return json.loads((HERE/"data"/f"derived_{name}.json").read_text())


def expressions(record):
    return {s.Symbol(k):s.sympify(v) for k,v in record["substitutions"].items()}


def branch(G2,G3,operator,profile=None,character=1,constant=0,**kw):
    return dict(g2=s.sympify(G2),g3=s.sympify(G3),operator={int(j):s.sympify(c) for j,c in operator.items()},
                profile=profile,character=character,constant=s.sympify(constant),**kw)


def model(slug,n,p,title,equation,parameter,interval,branches,notes,**kw):
    return dict(slug=slug,n=n,p=p,q=p//(n-1),title=title,equation=equation,parameter=parameter,
                interval=interval,branches=branches,notes=notes,count_kind="waves",primary=True,**kw)


def build_models():
    ordinary,twisted,tate,mixed=read("ordinary"),read("twisted"),read("tate"),read("mixed")
    all_models=[]
    add=all_models.append
    add(model("n2_p2_kdv",2,2,"KdV · the free-energy line",r"u''-u-6u^2=0",
        "energy coordinate  h",(-2.2,2.2),[branch(R(1,12),t/216,{0:-1},X-R(1,12))],
        ["One fixed equation carries a continuum of lattices; h = 216 g₃ is its first-integral parameter.",
         "The endpoints h = ±1 are degenerations. The open segment between them has Δ > 0."],fixed_equation=True))
    E=R(1,3)
    add(model("n3_p2_mkdv",3,2,"mKdV · a tilted energy line",r"u''-u-2u^3=0",
        "energy coordinate  h",(-1,2.5),[branch(t,4*E**3-t*E,{0:-1},"sqrt(X-1/3)",2)],
        ["Again the equation is fixed: the energy is free. The plotted lattice is the lattice of u².",
         "u = √(℘ − 1/3) has a character of order two; its ordinary period cell is twice as large."],fixed_equation=True))
    add(model("n2_p3_ks",2,3,"Kuramoto–Sivashinsky",r"v'''+4v''+v'+v^2/2+C=0",
        "integration constant  C",(-26,0),[branch(R(1,12),(t+13)/1080,{2:4,1:1},"-60*Y-60*X-1",constant=t,
            physical_profile=True)],
        ["One meromorphic wave per fixed C. Elliptic waves require b² = 16a in the general KS operator.",
         "C = −18 and −8 are the two degenerations; only the interval between them has Δ > 0."]))

    def tate_model(m,slug,title,interval,notes,**kw):
        rec=tate[f"tate_{m}"]
        # The coefficient of D^(m-1) is set to one. For order five this
        # chart excludes t=3; that excluded scaling is recorded explicitly.
        op={int(j):s.sympify(v) for j,v in rec["operator"].items()}
        lam=s.cancel(1/op[m-1])
        G2=s.factor(lam**4*s.sympify(rec["g2"]))
        G3=s.factor(lam**6*s.sympify(rec["g3"]))
        scaled={j:s.factor(lam**(m-j)*value) for j,value in op.items()}
        br=branch(G2,G3,scaled,rec["F"],m,profile_kind="tate_root",tate_order=m,rate=lam,
                  unscaled_F=rec["F"],unscaled_g2=rec["g2"],unscaled_g3=rec["g3"],unscaled_operator=rec["operator"])
        nonlinear=("+" if m%2 else "-")+str(s.factorial(m))
        return model(slug,m+1,m,title,rf"P(D)u{nonlinear}u^{{{m+1}}}=0,\quad a_{{{m-1}}}=1",
                     "torsion parameter  t",interval,[br],notes,**kw)

    add(tate_model(3,"n4_p3_quartic","Quartic flux · a threefold character",(-.035,.095),
        ["The normalized operator is D³ + D² + D/12 + (54t − 5)/108; one wave per equation.",
         "A parabola of lattices. Its two real degeneration parameters are t = 0 and t = 1/27."]))

    radical=s.sqrt(653016**2*t**4-4*662158224*(1457*t**4-28561))
    kaw=[]
    for sign in (-1,1):
        G2=(653016*t**2+sign*radical)/(2*662158224)
        G3=t*(31*t**2-42588*G2)/4745520
        profile=X**2+t*X/78-g2/10-(31*t**2+507)/851760
        kaw.append(branch(G2,G3,{2:t,0:-1},profile,label=("lower g₂" if sign<0 else "upper g₂")))
    add(model("n2_p4_kawahara",2,4,"Kawahara · the returning loop",r"u''''+a u''-u-840u^2=0",
        "dispersion coefficient  a",(-2.8,2.8),kaw,
        ["Two real lattices for |a| < 13/6. At each endpoint they merge into one degenerate wave.",
         "Beyond |a| = 13/6 the two lattices are complex conjugates and leave this real plane."],
        specials=[-R(13,6),R(13,6)],count_kind_override="waves",kawahara=True))

    def cubic_fourth(sign,primary=True):
        b=t/60; gamma=R(sign,360)-b**2
        Gp2=20*gamma; Gp3=20*b*(b**2-gamma)
        Gm2=20*b**2-R(40,3)*gamma; Gm3=-8*b**3+R(80,3)*b*gamma
        br=[branch(Gp2,Gp3,{2:t,0:sign},X+b,1,label="ordinary"),
            branch(Gm2,Gm3,{2:t,0:sign},"-D(sqrt(X-2*b))",2,label="twisted",half_root=2*b)]
        rec=model(f"n3_p4_cubic_{'positive' if sign>0 else 'negative'}",3,4,
            "Swift–Hohenberg · two lattice curves" if sign>0 else "Cubic fourth order · two tangencies",
            r"u''''+a u''"+("+u" if sign>0 else "-u")+r"-120u^3=0",
            "dispersion coefficient  a",(-4.6,4.6),br,
            ["Both character branches are shown, up to translation and sign.",
             "They coincide at a = ±5/2 on the cusp." if sign>0 else "The twisted branch touches the cusp at a = ±10/√11; the ordinary branch stays at Δ < 0."],
            specials=[-R(5,2),R(5,2)] if sign>0 else [-10/s.sqrt(11),10/s.sqrt(11)])
        rec["primary"]=primary
        return rec
    add(cubic_fourth(1))

    # Use the paper's especially simple a_2 parameter for the order-four character.
    G2=(1600*t**2-240*t-39)/120000
    G3=(40*t-9)*(1600*t**2-2160*t+441)/216000000
    F=X**2+Y/10+(40*t-3)*X/600+(1600*t**2-3120*t+657)/1440000
    add(model("n5_p4_quintic",5,4,"Quintic fourth order · a fourfold contact",
        r"u''''+u'''+a u''+a_1(a)u'+a_0(a)u-24u^5=0",
        "second-derivative coefficient  a",(-.2,.65),
        [branch(G2,G3,{3:1,2:t,1:(20*t-3)/200,0:3*(1200*t**2-1080*t+187)/40000},F,4,profile_kind="root_power")],
        ["One wave per compatible equation, with four ordinary poles and a character of order four.",
         "The contact at a = 11/40 has discriminant multiplicity four. The other cusp parameter is a = 1/4."],
        specials=[R(1,4),R(11,40)],zoom_parameter=(.22,.30)))

    for name,p,interval,title in [("quadratic_fifth",5,(-.2,.68),"Nikolaevskiy type · four algebraic branches"),
                                 ("quadratic_sixth",6,(-.05,.48),"Quadratic sixth order · a three-branch fold")]:
        rec=ordinary[name]; sub=expressions(rec)
        unknown=a2 if p==5 else g2
        poly=s.fraction(s.cancel(s.sympify(rec["remaining"][0])))[0]
        profile=s.sympify(rec["profile"]).subs(sub)
        op={4:1,3:t,2:a2,1:sub[s.Symbol("a1")]} if p==5 else {4:1,2:t}
        G2=sub.get(g2,g2); G3=sub[g3]
        K=s.rf(p,p)*(-1)**p
        eq=(r"u^{(5)}+u^{(4)}+a u'''+a_2u''+a_1u'+15120u^2+C=0" if p==5 else
            r"u^{(6)}+u^{(4)}+a u''-332640u^2+C=0")
        br=branch(G2,G3,op,profile,1,constant=sub[s.Symbol("C")],root_variable=unknown,root_polynomial=poly)
        m=model(f"n2_p{p}_{name}",2,p,title,eq,
                "coefficient  a",interval,[br],
                ["Every colored point gives a compatible equation; its remaining coefficients are read from the algebra.",
                 "Several points at the same a generally have different remaining coefficients, so they are different equations."],
                zoom_parameter=(.18,.36),zoom_limits=([[-2e-5,7e-5],[-9e-8,5e-8]] if p==5 else [[-5e-5,8e-5],[-3e-8,3e-8]]))
        m["count_kind"]="equations"
        add(m)
    # Put the p=5 sextic entry before the p=6 plates when sorting below.
    add(tate_model(5,"n6_p5_sextic","Sextic fifth order · a fivefold character",(-.5,20),
        ["A character of order five gives five poles in an ordinary period cell; each displayed equation has one wave.",
         "The chart fixes a₄ = 1. The value t = 3 has a₄ = 0 and is outside this normalization."],
        breaks=[3],specials=[0,(11-5*s.sqrt(5))/2,(11+5*s.sqrt(5))/2],
        zoom_parameter=(-.13,.045),view_limits=[[-.00015,.00125],[-.000018,.000003]]))

    rec=twisted["twisted_3_6"]; sub=expressions(rec)
    chart={d:(t+R(1,83))/2,e:t}
    poly=s.fraction(s.cancel(s.sympify(rec["remaining"][0]).subs(chart)))[0]
    br=branch(g2,4*t**3-t*g2,{j:s.factor(sub[s.Symbol(f"a{j}")].subs(chart)) for j in (0,2,4)},
              "(X+d)*sqrt(X-e)",2,root_variable=g2,root_polynomial=poly,
              profile_d=chart[d],half_root=t)
    m=model("n3_p6_cubic_twisted",3,6,"Cubic sixth order · a twisted bifurcation",
        r"u^{(6)}+u^{(4)}+a_2u''+a_0u-20160u^3=0",
        "half-period value  e",(-.025,.035),[br],
        ["The two roots for g₂ give coefficient families a₂(e), a₀(e); the displayed branches can be different equations.",
         "u = (℘ + (e + 1/83)/2)√(℘ − e). The ordinary branch requires a₄ = 0; see its separate plate."])
    m["count_kind"]="equations";add(m)

    rec=ordinary["quartic_sixth"];sub=expressions(rec);b=s.Symbol("b")
    chart={b:R(1,168),g2:t}
    br=branch(t,s.factor(sub[g3].subs(chart)),{j:s.factor(sub[s.Symbol(f"a{j}")].subs(chart)) for j in (0,2,4)},X+R(1,168))
    add(model("n4_p6_quartic",4,6,"Quartic sixth order · the parabolic branch",
        r"u^{(6)}+u^{(4)}+a_2u''+a_0u-5040u^4=0",
        "lattice parameter  g₂",(-.0006,.0017),[br],
        ["One wave per compatible equation; u = ℘ + 1/168. The curve is an exact parabola in the lattice plane.",
         "Here a₂ = 168g₂ + 5/28. The remaining coefficient a₀ is determined by the same substitution."]))
    add(tate_model(6,"n7_p6_septic","Septic sixth order · a sixfold character",(-1.5,1.5),
        ["The normalized coefficient a₅ = 1 is fixed; the other coefficients follow an exact one-parameter family.",
         "Characters of order six have six ordinary poles. The real cusp parameters are t = −1, −1/9, and 0."],
        specials=[-1,-R(1,9),0],zoom_parameter=(-1.05,.15)))

    # Companion slices prevent even / non-even and ordinary / twisted sectors
    # from being conflated in the primary, one-plate-per-pair overview.
    add(cubic_fourth(-1,False))
    rec=ordinary["quadratic_fourth_mixed"];sub=expressions(rec)
    m=model("n2_p4_mixed",2,4,"Kawahara–KS · the dissipative curve",
        r"u''''+u'''+a u''+a_1(a)u'-840u^2+C(a)=0",
        "second-derivative coefficient  a",(-.1,.65),
        [branch(sub[g2],sub[g3],{3:1,2:t,1:sub[s.Symbol("a1")]},s.sympify(rec["profile"]).subs(sub),constant=sub[s.Symbol("C")])],
        ["The odd derivative removes the two-wave freedom: one wave per compatible fixed equation.",
         "The locus resolves several different passages through the discriminant."],zoom_parameter=(.20,.42))
    m["primary"]=False;add(m)
    rec=mixed["cubic_fourth_mixed"];sub=expressions(rec);E=sub[e]
    m=model("n3_p4_mixed",3,4,"Swift–Hohenberg · the dispersive curve",
        r"u''''+u'''+a u''+a_1(a)u'+a_0(a)u-120u^3=0",
        "second-derivative coefficient  a",(-.12,.65),
        [branch(sub[g2],s.factor(4*E**3-sub[g2]*E),{3:1,2:t,1:sub[s.Symbol("a1")],0:sub[s.Symbol("a0")]},
            "-D(sqrt(X-e))-sqrt(X-e)/14",2,half_root=E,profile_b=-R(1,14))],
        ["The residue condition fixes a₁ = (196a − 39)/1372, then a₀ and the lattice are fixed.",
         "One wave per equation, with a character of order two."],zoom_parameter=(.2,.36))
    m["primary"]=False;add(m)
    m=model("n3_p6_ordinary",3,6,"Cubic sixth order · the ordinary line",
        r"u^{(6)}+a u''-u-20160u^3=0",
        "second-derivative coefficient  a",(-7,4),
        [branch(-t/252,R(1,4320),{2:t,0:-1},-Y/2)],
        ["With a₄ = 0 the elliptic wave is ordinary: u = −℘′/2, one wave per fixed equation.",
         "Its lattice follows a horizontal line; the pure D⁶ − 1 equation is the point g₂ = 0."])
    m["primary"]=False;add(m)
    for n,p,K,Ec in [(5,4,24,-R(1,30)),(7,6,720,-R(1,105))]:
        rec=twisted[f"twisted_{n}_{p}"];sub=expressions(rec)
        op={j:s.factor(sub[s.Symbol(f"a{j}")].subs({e:Ec,g2:t})) for j in range(0,p,2)}
        m=model(f"n{n}_p{p}_even",n,p,"Quintic · the even-operator line" if n==5 else "Septic · the even-operator line",
            rf"u^{{({p})}}+u^{{({p-2})}}+\cdots-{K}u^{{{n}}}=0",
            "lattice parameter  g₂",(-.009,.019) if n==5 else (-.001,.002),
            [branch(t,4*Ec**3-t*Ec,op,f"sqrt(X-({Ec}))",2)],
            ["Even operators also carry characters of order two. This branch is separate from the full-order character family.",
             "Every point determines the compatible lower coefficients; each such equation has one meromorphic wave."])
        m["primary"]=False;add(m)
    all_models.sort(key=lambda m:(not m["primary"],m["p"],m["n"],m["slug"]))
    expected={(n,p) for p in range(2,7) for n in range(2,8) if p%(n-1)==0}
    assert {(m["n"],m["p"]) for m in all_models if m["primary"]}==expected
    return all_models


def serializable(model):
    if isinstance(model,dict): return {str(k):serializable(v) for k,v in model.items()}
    if isinstance(model,(list,tuple)): return [serializable(v) for v in model]
    if isinstance(model,s.Basic): return str(model)
    return model


if __name__=="__main__":
    models=build_models()
    (HERE/"data"/"models.json").write_text(json.dumps(serializable(models),indent=2,ensure_ascii=False)+"\n")
    print(f"{len(models)} coefficient slices; {sum(m['primary'] for m in models)} distinct (n,p) pairs")
