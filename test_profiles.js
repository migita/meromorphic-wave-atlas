/* Evaluate independent fixtures for the Python jet/ODE checker. */
const fs=require('node:fs'),M=require('./wave_math.js');
const fixtures=JSON.parse(fs.readFileSync('data/profile_fixtures.json','utf8'));
const samples=[],failures=[];let periodChecks=0,zeroChecks=0;
for(let index=0;index<fixtures.length;index++){
  const f=fixtures[index],g=f.profile.geometry,L=g.period||8/Math.sqrt(g.G||g.H||1);
  for(const requested of ['regular','axis']){
    const mode=M.actualMode(f.profile,requested),drawn=M.generate(f.profile,requested);
    if(!drawn.range.every(Number.isFinite)||!drawn.re.some(Number.isFinite))failures.push([f.slug,f.t,mode,'empty']);
    for(const fraction of [.137,.319,.571]){
      const x=L*fraction,base=M.base(g,x,mode),value=M.at(f.profile,x,mode);
      samples.push({fixture:index,mode,x,X:base.X,Y:base.Y,value});
      if(g.period){
        const other=M.at(f.profile,x+g.period*M.periodMultiplier(f.profile,mode),mode);
        const error=Math.hypot(value[0]-other[0],value[1]-other[1])/Math.max(Math.hypot(...value),Math.hypot(...other),1e-25);
        if(error>3e-7)failures.push([f.slug,f.t,mode,'period',error]);periodChecks++;
      }
    }
    if(f.profile.kind==='root'&&f.profile.zeros[mode]){
      const zero=f.profile.zeros[mode],x=zero.position??zero.fraction*g.period,h=L*1e-5;
      const a=M.at(f.profile,x-h,mode),b=M.at(f.profile,x+h,mode);
      if(a&&b){
        const error=Math.hypot(a[0]+b[0],a[1]+b[1])/Math.max(Math.hypot(...a),Math.hypot(...b),1e-20);
        if(error>.02)failures.push([f.slug,f.t,mode,'zero continuation',error]);zeroChecks++;
      }
    }
  }
}
fs.writeFileSync('data/profile_evaluations.json',JSON.stringify({samples,failures,periodChecks,zeroChecks}));
console.log(JSON.stringify({samples:samples.length,periodChecks,zeroChecks,failures:failures.slice(0,20),totalFailures:failures.length}));
if(failures.length)process.exitCode=1;
