"""Check JavaScript C4 solutions by an independent full ODE substitution."""
from pathlib import Path
import json
import subprocess
import mpmath as mp
from profiles import geometry,base_at

HERE=Path(__file__).resolve().parent


def main():
    subprocess.run(['node','test_c4.js'],cwd=HERE,check=True)
    data=json.loads((HERE/'data/c4_numeric_fixtures.json').read_text());errors=[];failures=[];lattice_errors=[]
    mp.mp.dps=50
    geometries=[geometry(mp.mpf(f['point']['g2']),mp.mpf(f['point']['g3']),f['point']['profile']['geometry'].get('node',False))
                for f in data['fixtures']]
    for f,g in zip(data['fixtures'],geometries):
        L=f['point']['profile']['geometry']['period']
        assert (L is None)==(g['period'] is None)
        if L is not None:assert abs(L-g['period'])<2e-8*g['period']
    for sample in data['evaluations']:
        f=data['fixtures'][sample['fixture']];c=f['coeff'];p=f['point']
        al=complex(*p['alpha']);be=complex(*p['beta']);X=sample['X'];Y=sample['Y']
        refX,refY=base_at(geometries[sample['fixture']],mp.mpf(sample['x']),sample['mode']=='regular')
        lattice_scale=max(abs(p['g2'])**.5,abs(p['g3'])**(1/3),1e-100)
        lattice_error=max(abs(X-float(refX))/max(abs(float(refX)),lattice_scale),
                          abs(Y-float(refY))/max(abs(float(refY)),lattice_scale**1.5))
        lattice_errors.append(lattice_error)
        if lattice_error>2e-8:failures.append(dict(name=f['name'],lattice_error=lattice_error))
        v=complex(*sample['value']);v1=al*Y;v2=al*(6*X*X-p['g2']/2)
        v4=al*(120*X**3-18*p['g2']*X-12*p['g3'])
        expected=al*X+be
        assert abs(v-expected)<2e-10*max(abs(v),abs(expected),1e-20)
        terms=[v4,c['A']*v*v2,c['B']*v1*v1,c['C']*v**3,c['a2']*v2,c['a0']*v,c['a00']]
        rate=max(abs(p['g2'])**.25,abs(p['g3'])**(1/6))
        size=abs(al)*rate**2+abs(be)
        natural=abs(al)*rate**6+abs(c['A'])*size*abs(al)*rate**4+abs(c['B'])*abs(al)**2*rate**6+abs(c['C'])*size**3
        natural+=abs(c['a2'])*abs(al)*rate**4+abs(c['a0'])*size+abs(c['a00'])
        scale=max(sum(abs(t) for t in terms),natural,1e-100)
        error=abs(sum(terms))/scale;errors.append(error)
        if error>5e-7:failures.append(dict(name=f['name'],coefficients=c,mode=sample['mode'],error=error))
    result=dict(status='PASS' if not failures else 'FAIL',fixtures=len(data['fixtures']),profile_evaluations=len(errors),
        sweep_points=data['sweepPoints'],max_scaled_ode_residual=max(errors),
        max_scaled_lattice_error=max(lattice_errors),independent_period_checks=len(geometries),failures=failures)
    (HERE/'data/c4_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    if failures:raise SystemExit(1)


if __name__=='__main__':main()
