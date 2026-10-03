"""Compile the interactive site and a well resolved Kawahara pulse approach.

Near the double root, keep the complementary Jacobi parameter mc explicitly:
computing 1-m in JavaScript would erase the long-period states. No pictures
are generated. All coefficients remain in the displayed physical normalization.
"""
import copy
import json
from pathlib import Path
import mpmath as mp

HERE = Path(__file__).resolve().parent


def build_journey():
    records = []
    with mp.workdps(110):
        critical = mp.mpf(13)/6
        start_distance = critical-mp.mpf('1.45')
        for i in range(160):
            distance = start_distance*mp.power(10, -mp.mpf(i)*30/159)
            a = -critical+distance
            radical = mp.sqrt(653016**2*a**4-4*662158224*(1457*a**4-28561))
            g2 = (653016*a**2+radical)/(2*662158224)
            g3 = a*(31*a**2-42588*g2)/4745520
            r = mp.sqrt(g2/12)
            theta = mp.acos(3*mp.sqrt(3)*g3/g2**mp.mpf('1.5'))/3
            e = [2*r*mp.cos(theta-j*2*mp.pi/3) for j in (0, 1, -1)]
            G = e[0]-e[2]
            D = e[1]-e[2]
            mc = (e[0]-e[1])/G
            m = D/G
            K = mp.ellipk(m)
            period = 2*K/mp.sqrt(G)
            b = -g2/10-(31*a*a+507)/851760
            value = lambda x: -1680*(x*x+a*x/78+b)
            extremes = [value(e[2]), value(e[1])]
            vertex = -a/156
            if e[2] < vertex < e[1]:
                extremes.append(value(vertex))
            mean = -1680*(g2/12+a*(e[0]-G*mp.ellipe(m)/K)/78+b)
            geometry = dict(type='rect', delta=1, node=False, e=list(map(float, e)),
                            G=float(G), D=float(D), m=float(m), mc=float(mc), K=float(K), period=float(period))
            spec = dict(kind='polynomial', scale=-1680, variable='v', geometry=geometry,
                        terms=[[2, 0, 1], [1, 0, float(a/78)], [0, 0, float(b)]])
            records.append(dict(t=float(a), g2=float(g2), g3=float(g3), delta=float(g2**3-27*g3**2),
                                distance=float(distance), a_decimal=mp.nstr(a, 100),
                                g2_decimal=mp.nstr(g2, 100), g3_decimal=mp.nstr(g3, 100),
                                period=float(period), amplitude=float(max(extremes)-min(extremes)),
                                mean=float(mean), profile=spec))
    atlas = json.loads((HERE/'data/atlas.json').read_text())
    model = next(m for m in atlas if m['slug'] == 'n2_p4_kawahara')
    pulse = min(model['frames'], key=lambda f: abs(f['t']+13/6))['points'][0]
    records.append(dict(t=-13/6, g2=pulse['g2'], g3=pulse['g3'], delta=0,
                        distance=0, a_decimal='-13/6', period=None, amplitude=35/12,
                        mean=None, profile=copy.deepcopy(pulse['profile'])))
    result = dict(slug='n2_p4_kawahara', width=60, records=records)
    (HERE/'data/journey.json').write_text(json.dumps(result, separators=(',', ':'), allow_nan=False)+'\n')
    return result


def compile_site(payload=None, journey=None):
    if payload is None:
        payload = json.loads((HERE/'data/atlas.json').read_text())
    if journey is None:
        path = HERE/'data/journey.json'
        journey = json.loads(path.read_text()) if path.exists() else build_journey()
    # Keep data in HTML so file:// viewing works without a web server or fetch.
    template = (HERE/'gallery_template.html').read_text()
    for placeholder, value in [('__ATLAS_DATA__', payload), ('__JOURNEY_DATA__', journey)]:
        template = template.replace(placeholder, json.dumps(value, separators=(',', ':'), ensure_ascii=False, allow_nan=False))
    (HERE/'index.html').write_text(template)
    print('Compiled interactive atlas and', len(journey['records']), 'Kawahara limit states', flush=True)


if __name__ == '__main__':
    compile_site(journey=build_journey())
