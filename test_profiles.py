"""Independent ODE residual checks using local power-series arithmetic.

The browser draws Jacobi functions. This checker instead builds jets directly
from X'=Y, Y'=6X^2-g2/2, recovers the wave's Taylor coefficients, and substitutes
into the displayed high-order ODE. Points too close to zeros for a fractional
power's floating-point jet are counted separately; continuation is tested in JS.
"""
from pathlib import Path
import json
import math
import subprocess

HERE=Path(__file__).resolve().parent


def multiply(a,b,N):
    return [sum(a[k]*b[j-k] for k in range(j+1) if k<len(a) and j-k<len(b)) for j in range(N+1)]


def power(a,n,N):
    result=[1]+[0]*N
    for _ in range(n):result=multiply(result,a,N)
    return result


def polynomial(terms,X,Y,N):
    result=[0]*(N+1)
    for i,j,c in terms:
        monomial=multiply(power(X,i,N),power(Y,j,N),N)
        result=[a+c*b for a,b in zip(result,monomial)]
    return result


def main():
    process=subprocess.run(['node','test_profiles.js'],cwd=HERE,check=False)
    fixtures=json.loads((HERE/'data/profile_fixtures.json').read_text())
    evaluations=json.loads((HERE/'data/profile_evaluations.json').read_text())
    failures=list(evaluations['failures']);ode_errors=[];power_errors=[];curve_errors=[];ill_conditioned=0
    for sample in evaluations['samples']:
        f=fixtures[sample['fixture']];p=f['p'];m=f['n']-1;xx=sample['X'];yy=sample['Y']
        u=complex(*sample['value'])/f['profile']['scale']
        curve_terms=[yy**2,-4*xx**3,f['g2']*xx,f['g3']]
        curve_error=abs(sum(curve_terms))/max(sum(abs(c) for c in curve_terms),1e-150)
        curve_errors.append(curve_error)
        X=[xx,yy]
        for j in range(p):X.append((6*sum(X[k]*X[j-k] for k in range(j+1))-(f['g2']/2 if j==0 else 0))/((j+2)*(j+1)))
        Y=[(k+1)*X[k+1] for k in range(p+1)]
        F=polynomial(f['power'],X,Y,p)
        norm=sum(abs(c*xx**i*yy**j) for i,j,c in f['power'])
        error=abs(u**m-F[0])/max(norm,abs(u)**m,1e-150)
        power_errors.append(error)
        if error>2e-6:failures.append([f['slug'],f['t'],sample['mode'],'power identity',error])
        if f['profile']['kind']=='polynomial':
            jet=polynomial(f['profile']['terms'],X,Y,p)
        elif f['profile']['kind']=='halfroot':
            spec=f['profile'];G0=complex(xx-spec['root'])**.5
            if abs(G0)<1e-14*max(abs(xx)**.5,1e-30):
                ill_conditioned+=1;continue
            gj=[G0]
            for k in range(1,p+2):gj.append((X[k]-sum(gj[j]*gj[k-j] for j in range(1,k)))/(2*G0))
            if spec['form']=='generator':jet=gj[:p+1]
            elif spec['form']=='derivative':jet=[-(k+1)*gj[k+1]+spec['b']*gj[k] for k in range(p+1)]
            else:
                multiplier=list(X);multiplier[0]+=spec['d'];jet=multiply(multiplier,gj,p)
            if abs(-jet[0]-u)<abs(jet[0]-u):jet=[-v for v in jet]
        else:
            if abs(F[0])<1e-8*norm or abs(u)<1e-80:
                ill_conditioned+=1;continue
            jet=[u]
            for k in range(1,p+1):
                known=power(jet,m,k)[k]
                jet.append((F[k]-known)/(m*u**(m-1)))
        nonlinear=.5 if f.get('physical') else -f['K']
        terms=[math.factorial(p)*jet[p],nonlinear*jet[0]**f['n'],f['constant']]
        terms.extend(float(a)*math.factorial(int(j))*jet[int(j)] for j,a in f['operator'].items())
        rate=abs(f['g2'])**.25+abs(f['g3'])**(1/6)
        natural=max(abs(f['K']),1)*rate**(p+p/(f['n']-1))
        residual=abs(sum(terms))/max(sum(abs(c) for c in terms),natural,1e-150)
        ode_errors.append(residual)
        if residual>3e-4:failures.append([f['slug'],f['t'],sample['mode'],'ODE',residual])
    result=dict(status='PASS' if not failures else 'FAIL',fixtures=len(fixtures),samples=len(evaluations['samples']),
        ode_checks=len(ode_errors),ill_conditioned_fractional_jets=ill_conditioned,
        max_ode_residual=max(ode_errors),max_power_error=max(power_errors),max_curve_error=max(curve_errors),
        period_checks=evaluations['periodChecks'],zero_continuation_checks=evaluations['zeroChecks'],failures=failures)
    (HERE/'data/profile_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({**result,'failures':failures[:25]},indent=2))
    if failures:raise SystemExit(1)


if __name__=='__main__':main()
