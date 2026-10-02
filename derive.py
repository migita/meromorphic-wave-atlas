"""Exact algebra behind the lattice atlas; no numerical fitting of curves."""

from pathlib import Path
import json
import sys
import sympy as s

HERE = Path(__file__).resolve().parent
from elliptic_divisor import torsion_function

X, Y, g2, g3, t, e, d = s.symbols("X Y g2 g3 t e d")


def reduce(expr, G2=g2, G3=g3):
    """Coordinate ring of Y**2 = 4 X**3 - g2 X - g3."""
    return s.Poly(s.expand(expr), Y).rem(s.Poly(Y**2 - 4*X**3 + G2*X + G3, Y)).as_expr().expand()


def derivative(expr, G2=g2, G3=g3):
    return reduce(s.diff(expr, X)*Y + s.diff(expr, Y)*(6*X**2-G2/2), G2, G3)


def equations(expr):
    return [s.factor(c) for _, c in sorted(s.Poly(expr, X, Y).terms(), key=lambda item: 2*item[0][0]+3*item[0][1], reverse=True) if c != 0]


def ordinary(n, p, profile, coefficients, constant=0):
    q = p//(n-1)
    K = (-1)**p*s.rf(q,p)
    derivatives = [profile]
    for _ in range(p):
        derivatives.append(derivative(derivatives[-1]))
    residual = reduce(derivatives[p] + sum(coefficients.get(j,0)*derivatives[j] for j in range(p)) - K*reduce(profile**n) + constant)
    return residual


def solve_triangular(rows, variables):
    sub = {}
    for variable in variables:
        rows = [s.factor(row.subs(sub)) for row in rows]
        rows = [row for row in rows if row != 0]
        row = next(row for row in rows if row.has(variable))
        roots = s.solve(row, variable)
        if len(roots) != 1:
            raise ValueError((variable, row, roots))
        sub[variable] = roots[0]
        print("  eliminated", variable, flush=True)
    return {key:s.factor(value.subs(sub)) for key,value in sub.items()}, [s.factor(row.subs(sub)) for row in rows if s.factor(row.subs(sub)) != 0]


def ordinary_families():
    a0,a1,a2,a3,a4,b,c,h,C=s.symbols("a0 a1 a2 a3 a4 b c h C")
    cases = [
        ("quadratic_fifth", 2,5,-X*Y/2+b*X**2+c*Y+d*X+h,{4:1,3:t,2:a2,1:a1},C,[b,c,d,a1,h,g3,g2,C]),
        ("quadratic_sixth",2,6,X**3+b*X**2+c*X+h,{4:1,2:t},C,[b,c,h,g3,C]),
        ("quartic_sixth",4,6,X+b,{4:a4,2:a2,0:a0},0,[a4,a2,a0,g3]),
        ("cubic_sixth_ordinary",3,6,-Y/2,{4:a4,2:a2,0:a0},0,[a4,a2,a0]),
        ("quadratic_fourth_mixed",2,4,X**2+b*Y+c*X+h,{3:1,2:t,1:a1},C,[b,c,a1,h,g2,g3,C]),
    ]
    out={}
    for name,n,p,u,coeffs,constant,variables in cases:
        rows=equations(ordinary(n,p,u,coeffs,constant))
        print(name, "rows:", rows, flush=True)
        try:
            sub,left=solve_triangular(rows,variables)
            print("solved:",sub,"remaining:",left,flush=True)
            out[name]={"profile":str(u),"substitutions":{str(k):str(v) for k,v in sub.items()},"remaining":[str(row) for row in left],"rows":[str(row) for row in rows]}
        except (ValueError,StopIteration) as exc:
            print("triangular stop:",exc,flush=True)
            out[name]={"rows":[str(row) for row in rows]}
    return out


def twisted_even():
    # g=sqrt(X-e), g3=4e^3-e*g2, g''=(2X+e)g.
    # D(g*A+g'*B)=g*(D A+(2X+e)B)+g'*(A+D B).
    a0,a2,a4=s.symbols("a0 a2 a4")
    cubic=4*X**3-g2*X-(4*e**3-e*g2)
    out={}
    for n,p,A in [(3,6,X+d),(5,4,s.Integer(1)),(7,6,s.Integer(1))]:
        G3=4*e**3-e*g2
        # Even derivatives divided by g are polynomials in X.
        def D2(poly):
            return s.expand(s.diff(poly,X,2)*cubic+s.diff(poly,X)*(6*X**2-g2/2+cubic/(X-e))+(2*X+e)*poly).cancel().expand()
        ds=[A]
        for _ in range(p//2): ds.append(D2(ds[-1]))
        coeff={0:a0,2:a2,4:a4}
        K=s.rf(p//(n-1),p)
        residual=s.expand(ds[-1]+sum(coeff[j]*ds[j//2] for j in range(0,p,2))-K*A**n*(X-e)**((n-1)//2))
        rows=list(reversed(s.Poly(residual,X).all_coeffs()))[::-1]
        rows=[s.factor(row) for row in rows if row != 0]
        sub,left=solve_triangular(rows,list(reversed([coeff[j] for j in range(0,p,2)])))
        print("twisted",n,p,sub,left,flush=True)
        out[f"twisted_{n}_{p}"]={"profile_over_sqrt":str(A),"substitutions":{str(k):str(v) for k,v in sub.items()},"remaining":[str(row) for row in left],"rows":[str(row) for row in rows]}
    return out


def tate_families():
    out={}
    for m in (3,4,5,6):
        if m==3: a1,a2,a3=s.Integer(1),s.Integer(0),t
        else:
            c={4:s.Integer(0),5:t,6:t}[m]
            b={4:t,5:t,6:t+t**2}[m]
            a1,a2,a3=1-c,-b,-b
        B2=a1**2+4*a2
        B4=a1*a3
        B6=a3**2
        G2=s.factor((B2**2-24*B4)/12)
        G3=s.factor((-B2**3+36*B2*B4-216*B6)/216)
        A,B=B2/12,a3
        F,chain=torsion_function(G2,G3,(A,B),m,domain=s.QQ.frac_field(t))
        # Imported module has symbols named X,Y as well.
        F=s.factor(F)
        print("Tate",m,"g2",G2,"g3",G3,"F",F,flush=True)
        z=s.Symbol("z")
        wp=z**-2+G2*z**2/20+G3*z**4/28+G2**2*z**6/1200
        fjet=s.series(F.subs({X:wp,Y:s.diff(wp,z)})*z**m,z,0,m+1).removeO().expand()
        # Truncated polynomial arithmetic avoids nested symbolic series growth.
        if s.expand(fjet).coeff(z,0)!=1: raise ValueError("normalization")
        def power_jet(base, exponent, degree):
            result=[s.Integer(1)]+[s.Integer(0)]*degree
            for _ in range(exponent):
                result=[s.expand(sum(result[k]*base[j-k] for k in range(j+1) if j-k<len(base))) for j in range(degree+1)]
            return result
        rootjet=[s.Integer(1)]
        for j in range(1,m+1):
            known=power_jet(rootjet,m,j)[j]
            rootjet.append(s.expand((fjet.coeff(z,j)-known)/m))
        root=sum(s.factor(c)*z**(j-1) for j,c in enumerate(rootjet))
        K=(-1)**m*s.factorial(m)
        power=sum(c*z**(j-m-1) for j,c in enumerate(power_jet(rootjet,m+1,m)))
        residual=s.diff(root,z,m)-K*power
        operator={}
        for j in range(m-1,-1,-1):
            value=s.factor(-s.expand(residual).coeff(z,-j-1)/((-1)**j*s.factorial(j)))
            operator[j]=value
            residual+=value*s.diff(root,z,j)
        print("operator",operator,flush=True)
        disc=s.factor(G2**3-27*G3**2)
        out[f"tate_{m}"]={"g2":str(G2),"g3":str(G3),"F":str(F),"point":[str(A),str(B)],"operator":{str(j):str(v) for j,v in operator.items()},"discriminant":str(disc)}
    return out


def mixed_cubic():
    a0,a1,b=s.symbols("a0 a1 b")
    cubic=4*X**3-g2*X-4*e**3+e*g2
    def add(a,b): return tuple(s.cancel(x+y) for x,y in zip(a,b))
    def scale(a,f): return tuple(s.cancel(f*x) for x in a)
    def mul(a,b): return (s.cancel(a[0]*b[0]+cubic*a[1]*b[1]),s.cancel(a[0]*b[1]+a[1]*b[0]))
    def deriv(a): return (s.cancel(s.diff(a[1],X)*cubic+a[1]*(6*X**2-g2/2)),s.diff(a[0],X))
    log=(s.Integer(0),1/(2*(X-e)))
    u=(b,-1/(2*(X-e)))
    ds=[u]
    for _ in range(4): ds.append(add(deriv(ds[-1]),mul(log,ds[-1])))
    residual=add(add(add(add(ds[4],ds[3]),scale(ds[2],t)),scale(ds[1],a1)),scale(ds[0],a0))
    residual=add(residual,scale(mul(mul(u,u),u),-120*(X-e)))
    rows=[]
    for component in residual:
        rows.extend([s.factor(c) for c in s.Poly(s.fraction(component)[0],X).all_coeffs() if c!=0])
    rows=[rows[j] for j in (0,3,1,4,2)]
    sub,left=solve_triangular(rows,[b,e,a1,a0,g2])
    print(sub,left,flush=True)
    return {"cubic_fourth_mixed":{"profile":"-D(sqrt(X-e)) + b*sqrt(X-e)","substitutions":{str(k):str(v) for k,v in sub.items()},"remaining":[str(x) for x in left],"rows":[str(x) for x in rows]}}


if __name__=="__main__":
    mode=sys.argv[1] if len(sys.argv)>1 else "ordinary"
    answer={"ordinary":ordinary_families,"twisted":twisted_even,"tate":tate_families,"mixed":mixed_cubic}[mode]()
    (HERE/"data"/f"derived_{mode}.json").write_text(json.dumps(answer,indent=2)+"\n")
