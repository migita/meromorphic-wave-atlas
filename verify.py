"""Independent checks of the plotted formulae, using the power u^(n-1).

The torsion curves were constructed from divisors and Laurent matching.
Here a separate rational-function calculation checks the full ODE in the
Weierstrass function field. Implicit branches are checked modulo their
polynomial, including complex roots, without choosing numerical roots.
"""
import json
import sys
import time
from pathlib import Path
import sympy as s

from models import HERE,build_models,serializable,t,X,Y,g2,g3


class Curve:
    def __init__(self,G2,G3):
        self.cubic=4*X**3-G2*X-G3
        self.half_derivative=6*X**2-G2/2

    def pair(self,f):
        f=s.Poly(s.expand(f),Y).rem(s.Poly(Y**2-self.cubic,Y)).as_expr().expand()
        return f.coeff(Y,0),f.coeff(Y,1)

    def add(self,a,b): return tuple(s.cancel(x+y) for x,y in zip(a,b))
    def scale(self,a,b): return tuple(s.cancel(x*b) for x in a)
    def multiply(self,a,b):
        return s.cancel(a[0]*b[0]+self.cubic*a[1]*b[1]),s.cancel(a[0]*b[1]+a[1]*b[0])
    def inverse(self,a):
        norm=s.factor(a[0]**2-self.cubic*a[1]**2)
        return s.cancel(a[0]/norm),s.cancel(-a[1]/norm)
    def derivative(self,a):
        return s.cancel(s.diff(a[1],X)*self.cubic+a[1]*self.half_derivative),s.diff(a[0],X)


def power_expression(model,br):
    """Return u^(n-1), or u itself for the quadratic equation."""
    n,p=model["n"],model["p"]
    kind=br.get("profile_kind")
    if kind=="root_power": return br["profile"]
    if kind=="tate_root":
        lam=br["rate"]; m=n-1
        return s.expand(lam**m*s.sympify(br["unscaled_F"]).subs({X:X/lam**2,Y:Y/lam**3}))
    profile=br["profile"]
    if not isinstance(profile,str): return profile**(n-1)
    if br.get("physical_profile"): return (s.sympify(profile)/120)**(n-1)
    if n==3 and p==4 and br["character"]==2:
        E=br["half_root"]; B=br.get("profile_b",s.Integer(0))
        # (-g' + B g)^2, g^2=X-E and E a root of the cubic.
        return X**2+E*X+E**2-br["g2"]/4 - B*Y+B**2*(X-E)
    if n==3 and p==6:
        return (X+br["profile_d"])**2*(X-br["half_root"])
    if profile.startswith("sqrt("):
        return s.sympify(profile[5:-1])**s.Rational(n-1,2)
    raise ValueError((model["slug"],profile))


def zero_on_root_polynomial(value,root_poly,variable):
    num=s.Poly(s.fraction(s.cancel(value))[0],X,Y)
    divisor=s.Poly(root_poly,variable)
    return all(s.Poly(coef,variable).rem(divisor).is_zero for coef in num.coeffs())


def verify_at(model,br,T):
    n,p=model["n"],model["p"]; q=p//(n-1); K=(-1)**p*s.rf(q,p)
    substitutions={t:T}
    G2=s.simplify(br["g2"].subs(substitutions))
    G3=s.simplify(br["g3"].subs(substitutions))
    op={j:s.simplify(v.subs(substitutions)) for j,v in br["operator"].items()}
    F=power_expression(model,br)
    F=s.expand(F.subs(substitutions).subs({g2:G2,g3:G3}))
    C=s.simplify(br["constant"].subs(substitutions))
    if br.get("physical_profile"): C/=120
    curve=Curve(G2,G3)
    Fpair=curve.pair(F)
    if n==2:
        derivatives=[Fpair]
        for _ in range(p): derivatives.append(curve.derivative(derivatives[-1]))
        residual=derivatives[p]
        for j,a in op.items(): residual=curve.add(residual,curve.scale(derivatives[j],a))
        residual=curve.add(residual,curve.scale(curve.multiply(Fpair,Fpair),-K))
        residual=curve.add(residual,(C,s.Integer(0)))
    else:
        psi=curve.scale(curve.multiply(curve.derivative(Fpair),curve.inverse(Fpair)),s.Rational(1,n-1))
        ratios=[(s.Integer(1),s.Integer(0))]
        for _ in range(p): ratios.append(curve.add(curve.derivative(ratios[-1]),curve.multiply(psi,ratios[-1])))
        residual=ratios[p]
        for j,a in op.items(): residual=curve.add(residual,curve.scale(ratios[j],a))
        residual=curve.add(residual,curve.scale(Fpair,-K))
        assert C==0
    if "root_variable" in br:
        poly=br["root_polynomial"].subs(t,T)
        valid=all(zero_on_root_polynomial(r,poly,br["root_variable"]) for r in residual)
    else:
        valid=all(s.simplify(r)==0 for r in residual)
    assert valid,(model["slug"],T,residual)
    return dict(parameter=str(T),identity="P(D)u - K*u^n + C = 0",result="exact zero",
                includes_all_algebraic_roots="root_variable" in br)


def check_degenerations():
    # Check the ODE directly, and prove that the root is meromorphic in z by
    # substituting w=exp(kz) into the nodal Weierstrass parametrization.
    # Every finite nonzero zero/pole of F(w) must have order divisible by m;
    # a power w^(r/m) is an entire exponential in z and is allowed.
    tate=json.loads((HERE/"data"/"derived_tate.json").read_text())
    out=[]
    w=s.Symbol("w")
    for m,values in [(3,[0,s.Rational(1,27)]),(4,[0,-s.Rational(1,16)]),(5,[0]),(6,[0,-1,-s.Rational(1,9)])]:
        for T in values:
            rec=tate[f"tate_{m}"]
            G2=s.sympify(rec["g2"]).subs(t,T);G3=s.sympify(rec["g3"]).subs(t,T)
            F=s.sympify(rec["F"]).subs(t,T)
            assert s.simplify(G2**3-27*G3**2)==0
            e=s.cancel(-3*G3/(2*G2));k=s.sqrt(12*e)
            W=e+12*e*w/(w-1)**2
            Yw=s.cancel(k*w*s.diff(W,w))
            rational=s.cancel(F.subs({X:W,Y:Yw}))
            assert rational!=0 and s.diff(rational,w)!=0
            numerator,denominator=s.fraction(rational)
            factors=[]
            for poly in (numerator,denominator):
                _,rows=s.factor_list(poly,w,extension=True)
                for factor,multiplicity in rows:
                    monomial=s.Poly(factor,w)
                    is_w=monomial.degree()==1 and monomial.nth(0)==0
                    assert is_w or multiplicity%m==0,(m,T,factor,multiplicity)
                    factors.append([str(factor),int(multiplicity)])
            br=dict(g2=G2,g3=G3,profile=F,profile_kind="root_power",constant=s.Integer(0),
                    operator={int(j):s.sympify(c).subs(t,T) for j,c in rec["operator"].items()})
            verify_at(dict(n=m+1,p=m,slug=f"degenerate_tate_{m}"),br,s.sympify(T))
            out.append(dict(n=m+1,p=m,t=str(T),power_in_exponential=str(rational),
                            factors=factors,global_meromorphic_root=True,ode_identity="exact zero"))
    return out


def main():
    results=[]
    start=time.monotonic()
    for model in build_models():
        samples=[s.Rational(0),s.Rational(1,10),s.Rational(-1,5)]
        if model.get("kawahara"): samples=[s.Rational(0),s.Rational(1),s.Rational(13,6)]
        # Use small values for a free lattice invariant; they need not lie in
        # the plotting window for the algebraic identity to be verified.
        checks=[]
        for br in model["branches"]:
            for T in samples:
                checks.append(verify_at(model,br,T))
        results.append(dict(slug=model["slug"],checks=checks))
        print("PASS",model["slug"],len(checks),"exact identities",flush=True)
    exceptional_checks=[]
    data=json.loads((HERE/"data"/"boundaries.json").read_text())
    for fibre in data["n2_p5_quadratic_fifth"]["exceptional_fibres"]:
        T=s.sympify(fibre["t"])
        sub={s.Symbol(k):s.sympify(v) for k,v in fibre["profile_substitutions"].items()}
        b,c,d,h=s.symbols("b c d h")
        profile=(-X*Y/2+b*X**2+c*Y+d*X+h).subs(sub)
        br={"g2":g2,"g3":s.sympify(fibre["g3"]),"operator":{4:s.Integer(1),3:T,2:s.sympify(fibre["a2"]),1:s.sympify(fibre["a1"])},
            "constant":s.sympify(fibre["C"]),"profile":profile,"character":1,"root_variable":g2,"root_polynomial":s.sympify(fibre["g2_polynomial"])}
        model={"n":2,"p":5,"slug":"quadratic_fifth_exceptional"}
        exceptional_checks.append(verify_at(model,br,T))
    print("PASS",len(exceptional_checks),"exceptional fibres",flush=True)
    # Negative control: moving the Kawahara point off its line fails.
    model=next(m for m in build_models() if m.get("kawahara"))
    bad=dict(model["branches"][0]);bad["g3"]=bad["g3"]+1
    try:
        verify_at(model,bad,s.Integer(0))
    except AssertionError:
        negative="off-locus point correctly rejected"
    else:
        raise AssertionError("The negative control was accepted")
    degeneration_checks=check_degenerations()
    output=dict(status="PASS",method=__doc__,models=results,negative_control=negative,
                degeneration_checks=degeneration_checks,exceptional_fibres=exceptional_checks,
                elapsed_seconds=round(time.monotonic()-start,2))
    (HERE/"data"/"verification.json").write_text(json.dumps(output,indent=2)+"\n")
    print("PASS",sum(len(r['checks']) for r in results),"identities;",len(degeneration_checks),"degenerate operators; negative control",flush=True)


if __name__=="__main__":main()
