"""Independent high-precision profile, ODE and period checks for the limit view."""
import json
import math
import subprocess
from pathlib import Path
import mpmath as mp
import sympy as s
from test_profiles import polynomial, multiply

HERE = Path(__file__).resolve().parent


def main():
    source = """
const W=require('./wave_math.js'),j=require('./data/journey.json');
const out=j.records.map(r=>({values:[-23.1,-12,-2.5,0,2.5,8,23.1].map(x=>({x,value:W.at(r.profile,x),...W.base(r.profile.geometry,x,'regular')})),
stats:(()=>{const p=W.generate(r.profile,'regular',2,801,.5);return {amplitude:p.amplitude,mean:p.mean,period:p.period,constant:p.constant};})()}));
process.stdout.write(JSON.stringify(out));
"""
    result = subprocess.run(['node', '-e', source], cwd=HERE, check=True, capture_output=True, text=True)
    values = json.loads(result.stdout)
    records = json.loads((HERE/'data/journey.json').read_text())['records']
    profile_errors, residuals, period_errors, mean_errors, amplitude_errors = [], [], [], [], []
    mp.mp.dps = 110
    # A direct hyperbolic substitution independently verifies the solitary ODE.
    T = s.Symbol('T')
    pulse = 2-s.Rational(35,12)*(1-T*T)**2
    derivative = lambda f: s.expand((1-T*T)*s.diff(f,T)/(2*s.sqrt(6)))
    assert s.expand(derivative(derivative(derivative(derivative(pulse))))-s.Rational(13,6)*derivative(derivative(pulse))-pulse+pulse*pulse/2)==0
    for index, (r, browser) in enumerate(zip(records, values)):
        if r['distance']:
            g2, g3, a = map(mp.mpf, [r['g2_decimal'], r['g3_decimal'], r['a_decimal']])
            radius = mp.sqrt(g2/12)
            theta = mp.acos(3*mp.sqrt(3)*g3/g2**mp.mpf('1.5'))/3
            e = [2*radius*mp.cos(theta-j*2*mp.pi/3) for j in (0,1,-1)]
            G, D = e[0]-e[2], e[1]-e[2]
            m = D/G
            period = 2*mp.ellipk(m)/mp.sqrt(G)
            def reference(z):
                X = e[2]+D*mp.ellipfun('sn', mp.sqrt(G)*z, m)**2
                return -1680*(X*X+a*X/78-g2/10-(31*a*a+507)/851760)
            period_errors.append(abs(float(period)-browser['stats']['period'])/float(period))
            if index in [0, 40, 80, 120, 159]:
                # Quadrature does not use the compiled analytic mean formula.
                mean = mp.quad(reference, [0, period/8, period/4, period/2])*2/period
                mean_errors.append(abs(float(mean)-r['mean']))
                mean_errors.append(abs(float(mean)-browser['stats']['mean']))
        else:
            reference = lambda z: 2-mp.mpf(35)/12/mp.cosh(z/(2*mp.sqrt(6)))**4
            assert browser['stats']['period'] is None and browser['stats']['mean'] is None
        amplitude_errors.append(abs(browser['stats']['amplitude']-r['amplitude']))
        assert not browser['stats']['constant'], 'A short display window was mistaken for a constant.'
        for sample in browser['values']:
            profile_errors.append(abs(float(reference(mp.mpf(sample['x'])))-sample['value'][0]))
            x, y = sample['X'], sample['Y']
            X = [x,y]
            for k in range(4):
                X.append((6*sum(X[j]*X[k-j] for j in range(k+1))-(r['g2']/2 if k==0 else 0))/((k+2)*(k+1)))
            Y = [(k+1)*X[k+1] for k in range(5)]
            jet = [r['profile']['scale']*v for v in polynomial(r['profile']['terms'],X,Y,4)]
            terms = [24*jet[4],r['t']*2*jet[2],-jet[0],jet[0]**2/2]
            residuals.append(abs(sum(terms))/max(sum(map(abs,terms)),1))
    assert max(profile_errors)<2e-10
    assert max(residuals)<2e-10
    assert max(period_errors)<2e-12
    assert max(mean_errors)<2e-8
    assert max(amplitude_errors)<2e-10
    report=dict(status='PASS',states=len(records),profile_values=len(profile_errors),
                exact_pulse_identity=True,short_window_statistics=True,
                max_profile_error=max(profile_errors),max_scaled_ode_residual=max(residuals),
                max_relative_period_error=max(period_errors),max_mean_error=max(mean_errors),
                max_amplitude_error=max(amplitude_errors))
    (HERE/'data/journey_checks.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
