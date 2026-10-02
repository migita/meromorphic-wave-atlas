const fs=require('node:fs'),C=require('./c4_math.js'),W=require('./wave_math.js');
function assert(value,message){if(!value)throw Error(message);}
const fixtures=[],states=[];
for(const key of Object.keys(C.presets)){
  const s=C.preset(key);
  for(const a2 of [-1,0,1])for(const a0 of [-2,-.1,1])for(const a00 of [0,1/7]){
    const c={...s.coeff,a2,a0,a00},free={beta:.13,g2:.07,g3:.0001},r=C.solve(c,free);states.push({c,free});
    for(const p of r.points)fixtures.push({name:key,coeff:c,point:p});
  }
}
const cases=[
 {A:1,B:1,C:0,a2:1,a0:0,a00:0},
 {A:3,B:2,C:1,a2:0,a0:-1,a00:0},
 {A:3,B:2,C:1,a2:0,a0:1,a00:0},
 {A:3,B:2,C:1,a2:0,a0:-1,a00:.1},
 {A:3,B:2,C:1,a2:0,a0:-1,a00:1},
 {A:4,B:1,C:4/3,a2:0,a0:-1,a00:0},
 {A:0,B:1,C:-4/27,a2:1,a0:.25,a00:0},
 {A:0,B:1,C:-.5,a2:1,a0:0,a00:4/27},
 {A:0,B:0,C:30,a2:1,a0:1,a00:0},
 {A:0,B:1,C:0,a2:0,a0:-1,a00:.2},
 {A:0,B:1,C:1/30,a2:1,a0:-1,a00:.15},
];
for(const coeff of cases)for(const g2 of [.03,.1])for(const g3 of [-.001,.001]){
  const free={beta:.12,g2,g3};states.push({c:coeff,free});
  for(const point of C.solve(coeff,free).points)fixtures.push({name:'resonant/custom',coeff,point});
}
const lax=C.oneGap(),base=C.solve(lax.coeff,lax.free);
assert(base.points.length===2,'Lax one-gap should show a family member and an isolated branch');
const branch=base.points.find(p=>p.alpha[0]===-2);
assert(branch.free.includes('beta')&&Math.abs(branch.g2-.5)<1e-12&&branch.g3===0,'Wrong one-gap coordinates');
const pulse=C.solve(lax.coeff,{...lax.free,beta:1/6}).points.find(p=>Math.abs(p.alpha[0]+2)<1e-12);
for(const x of [0,1,3,8])assert(Math.abs(W.at(pulse.profile,x,'regular')[0]-.5/Math.cosh(x/2)**2)<1e-11,'Incorrect C4 pulse');
assert(C.solve(C.preset('lax').coeff).points.every(p=>Math.abs(p.alpha[0]+2)>1e-8),'Resonant branch present for a2 nonzero');
assert(C.solve({A:0,B:0,C:0,a2:0,a0:1,a00:0}).points.length===0,'Linear equation has no affine-wp pole');
assert(C.solve({A:1,B:1,C:0,a2:1,a0:0,a00:.01}).points.length===0,'Missed compatibility condition');
const all=C.solve(cases[0]);assert(all.points.length===1&&all.points[0].free.length===2,'Missing two-dimensional family');
assert(C.parseNumber('20/3')===20/3&&C.parseNumber('0/0')===null&&C.parseNumber('1;alert(1)')===null,'Coefficient parser');
assert(C.solve(cases[cases.length-1]).points.length===1,'Repeated leading root must be counted once');
assert(C.fromParams(new URLSearchParams('preset=__proto__&sweep=constructor')).sweep==='beta','Invalid scene selectors');
const evaluations=[];
for(let i=0;i<fixtures.length;i++){
  const f=fixtures[i],g=f.point.profile.geometry,L=g.period||8/Math.sqrt(g.G||g.H||1);
  for(const requested of ['regular','axis']){
    const mode=W.actualMode(f.point.profile,requested);
    for(const fraction of [.137,.313]){
      const x=L*fraction,b=W.base(g,x,mode),v=W.at(f.point.profile,x,mode);
      assert(v&&v.every(Number.isFinite),'Nonfinite C4 profile');evaluations.push({fixture:i,x,mode,X:b.X,Y:b.Y,value:v});
    }
  }
}
let sweepPoints=0;
for(const s of [C.oneGap(),...Object.keys(C.presets).map(k=>C.preset(k))])for(const key of ['A','B','C','a2','a0','a00']){
  const state=JSON.parse(JSON.stringify(s));state.sweep=key;state.range=C.defaultRange(state,key);const m=C.build(state,121);
  assert(m.limits.flat().every(Number.isFinite),'Nonfinite C4 plot limits');
  for(const f of m.frames)for(const p of f.points){assert(p.residual<2e-8,'Unverified C4 sweep point');const v=W.at(p.profile,(p.profile.geometry.period||1)*.19,'regular');assert(v&&v.every(Number.isFinite),'Nonfinite sweep profile');sweepPoints++;}
}
fs.writeFileSync('data/c4_numeric_fixtures.json',JSON.stringify({fixtures,evaluations,sweepPoints}));
console.log(JSON.stringify({fixtures:fixtures.length,profileEvaluations:evaluations.length,sweepPoints,status:'PASS'}));
