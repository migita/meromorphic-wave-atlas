"""Elimination polynomials for implicit-branch cusp contacts and folds."""
import json
import sympy as s
from models import HERE,build_models,t,read,g2,g3,a2
from derive import solve_triangular


def factors(expr):
    return [str(s.Poly(f,t).monic().as_expr()) for f,_ in s.factor_list(expr,t)[1] if f.has(t)]


out={}
for model in build_models():
    for br in model["branches"]:
        if "root_variable" not in br:continue
        r=br["root_variable"];P=br["root_polynomial"]
        delta=s.factor(br["g2"]**3-27*br["g3"]**2)
        num=s.fraction(delta)[0]
        disc=s.discriminant(P,r)
        result=s.resultant(P,num,r)
        out[model["slug"]]={"fold_factors":factors(disc),"cusp_factors":factors(result),
            "leading_factors":factors(s.Poly(P,r).LC()),
            "note":"Resultant roots are candidates; check the original equations and exclude poles of the readout."}
        print(model["slug"],json.dumps(out[model["slug"]]),flush=True)

# Recover the finite fibres lost if a vanishing linear g2 coefficient is
# divided out in the fifth-order quadratic readout.
rec=read("ordinary")["quadratic_fifth"]
b,c,d,h,a1,C=s.symbols("b c d h a1 C")
sub,left=solve_triangular([s.sympify(row) for row in rec["rows"]],[b,c,d,a1,h,g3,C])
linear=next(row for row in left if s.degree(row,g2)==1)
L=s.Poly(s.fraction(s.cancel(linear))[0],g2)
r0=s.solve(L.nth(1),a2)[0]
crit=s.factor(L.nth(0).subs(a2,r0))
fibres=[]
for T in s.solve(crit,t):
    if T.is_real is not True:continue
    R0=s.factor(r0.subs(t,T));where={t:T,a2:R0}
    remaining=[s.factor(row.subs(where)) for row in left if s.factor(row.subs(where))!=0]
    poly=s.Poly(s.fraction(s.cancel(remaining[0]))[0],g2).monic().as_expr()
    fibres.append({"t":str(T),"a2":str(R0),"g2_polynomial":str(poly),
       "g3":str(s.factor(sub[g3].subs(where))),"C":str(s.factor(sub[C].subs(where))),
       "a1":str(s.factor(sub[a1].subs(where))),
       "profile_substitutions":{str(k):str(s.factor(sub[k].subs(where))) for k in (b,c,d,h)}})
out["n2_p5_quadratic_fifth"]["exceptional_fibres"]=fibres
print("recovered exceptional fibres:",json.dumps(fibres),flush=True)
(HERE/"data"/"boundaries.json").write_text(json.dumps(out,indent=2)+"\n")
