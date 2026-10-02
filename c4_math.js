/* C4 affine-Weierstrass sector: v = alpha*wp + beta.
 * The undivided coefficient equations are solved before a lattice is drawn.
 * Resonant denominators produce compatibility conditions and free coordinates,
 * not divisions by zero. Coefficients and free slider coordinates are real.
 */
(function(root){
  'use strict';
  const W=typeof module!=='undefined'?require('./wave_math.js'):root.WaveMath;
  const z=x=>Array.isArray(x)?x:[x,0],add=(a,b)=>{a=z(a);b=z(b);return [a[0]+b[0],a[1]+b[1]];},neg=a=>{a=z(a);return [-a[0],-a[1]];},sub=(a,b)=>add(a,neg(b));
  const mul=(a,b)=>{a=z(a);b=z(b);return [a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]];},abs=a=>{a=z(a);return Math.hypot(...a);};
  function div(a,b){a=z(a);b=z(b);if(Math.abs(b[0])>=Math.abs(b[1])){const r=b[1]/b[0],d=b[0]+b[1]*r;return [(a[0]+a[1]*r)/d,(a[1]-a[0]*r)/d];}const r=b[0]/b[1],d=b[1]+b[0]*r;return [(a[0]*r+a[1])/d,(a[1]*r-a[0])/d];}
  const pow=(a,n)=>{let r=[1,0];for(let j=0;j<n;j++)r=mul(r,a);return r;},eps=2e-12,zero=(a,scale)=>abs(a)<=eps*Math.abs(scale),real=a=>Math.abs(z(a)[1])<=2e-9*Math.max(abs(a),1e-30);
  const coefKeys=['A','B','C','a2','a0','a00'],freeKeys=['beta','g2','g3'];
  const labels={A:'A · v v″',B:'B · (v′)²',C:'C · v³',a2:'a₂ · v″',a0:'a₀ · v',a00:'a₀₀ · constant',beta:'β · free wave background',g2:'g₂ · free lattice coordinate',g3:'g₃ · free lattice coordinate'};
  const presets={
    lax:{name:'Lax core',core:[10,5,10],a0:-1,range:[-2,.25]},
    sk:{name:'Sawada–Kotera core',core:[15,0,15],a0:-2,range:[-3,.25]},
    kk:{name:'Kaup–Kupershmidt core',core:[10,7.5,20/3],a0:-1,range:[-8,.5]},
    cubic:{name:'Pure cubic',core:[0,0,-30],a0:.3,range:[-.6,.6]},
    generic:{name:'Generic core (5,30,30)',core:[5,30,30],a0:-1,range:[-3.2,.4]}
  };
  const hasPreset=key=>Object.prototype.hasOwnProperty.call(presets,key);
  function preset(key){const p=hasPreset(key)?presets[key]:presets.lax;return {coeff:{A:p.core[0],B:p.core[1],C:p.core[2],a2:1,a0:p.a0,a00:0},free:{beta:0,g2:.1,g3:.001},sweep:'a0',range:p.range.slice(),focus:key==='kk'||key==='generic'};}
  function oneGap(){const r=preset('lax');r.coeff.a2=0;r.sweep='beta';r.range=[-.21,.21];r.focus=false;return r;}
  function presetKey(c){return Object.keys(presets).find(key=>presets[key].core.every((v,j)=>Math.abs(v-c[coefKeys[j]])<1e-11*Math.max(1,Math.abs(v))))||'custom';}
  function parseNumber(text){const parts=String(text).trim().replace(/−/g,'-').split('/');if(parts.length>2||parts.some(p=>!/^[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:e[-+]?\d+)?$/i.test(p.trim())))return null;const v=Number(parts[0])/(parts.length===2?Number(parts[1]):1);return Number.isFinite(v)&&Math.abs(v)<=1e6?v:null;}
  const format=x=>Number(x).toPrecision(6).replace(/\.?0+(e|$)/,'$1');
  function complexText(a){a=z(a);const r=Math.abs(a[0])<1e-12*Math.max(abs(a),1)?0:a[0],i=Math.abs(a[1])<1e-12*Math.max(abs(a),1)?0:a[1];if(!i)return format(r);return (r?format(r)+(i<0?' − ':' + '):i<0?'−':'')+format(Math.abs(i))+'i';}
  function quadratic(a,b,c){
    if(a===0)return b===0?[]:[{value:[-c/b,0],id:b>0?0:1}];
    let d=b*b-4*a*c;if(zero(d,Math.abs(b*b)+Math.abs(4*a*c)))d=0;
    let values;if(d>=0){const r=Math.sqrt(d),p=-b+r,m=-b-r;values=[Math.abs(p)<Math.abs(m)?2*c/m:p/(2*a),Math.abs(m)<Math.abs(p)?2*c/p:m/(2*a)].map(v=>[v,0]);}
    else values=[[-b/(2*a),Math.sqrt(-d)/(2*a)],[-b/(2*a),-Math.sqrt(-d)/(2*a)]];
    return values.map((value,id)=>({value,id})).filter((v,i,a)=>v.value.every(Number.isFinite)&&!a.slice(0,i).some(q=>abs(sub(q.value,v.value))<eps*Math.max(abs(q.value),abs(v.value),1e-30)));
  }
  function cubicRealLattice(a,b,c){
    if(a===0)return b===0?[]:[[-c/b,0]];
    if(c===0)return [[0,0],...quadratic(a,0,b).map(x=>x.value)].filter((v,i,all)=>!all.slice(0,i).some(q=>abs(sub(v,q))<eps*Math.max(abs(v),1)));
    const p=b/a,q=c/a,D=q*q/4+p*p*p/27;
    if(D>0){const u=Math.cbrt(-q/2+Math.sqrt(D)),v=Math.cbrt(-q/2-Math.sqrt(D)),r=u+v,i=Math.sqrt(3)*(u-v)/2;return [[r,0],[-r/2,i],[-r/2,-i]];}
    const r=2*Math.sqrt(Math.max(0,-p/3));if(r===0)return [[Math.cbrt(-q),0]];
    const theta=Math.acos(Math.max(-1,Math.min(1,-q/(2*Math.sqrt(-p*p*p/27)))))/3;
    return [0,1,2].map(j=>[r*Math.cos(theta-j*2*Math.PI/3),0]).filter((v,i,all)=>!all.slice(0,i).some(q=>abs(sub(v,q))<eps*Math.max(abs(v),1)));
  }
  function constraints(c,alpha,beta,g2,g3){
    const {A,B,C,a2,a0,a00}=c,b2=pow(beta,2),b3=pow(beta,3),aa=abs(alpha),bb=abs(beta);
    const rows=[
      [add(add(mul(C,pow(alpha,2)),mul(6*A+4*B,alpha)),120),Math.abs(C)*aa*aa+Math.abs(6*A+4*B)*aa+120],
      [add(mul(add(2*A,mul(C,alpha)),beta),2*a2),(2*Math.abs(A)+Math.abs(C)*aa)*bb+2*Math.abs(a2)],
      [sub(mul(add(36,mul(A+2*B,alpha)),g2),add(mul(6*C,b2),2*a0)),(36+Math.abs(A+2*B)*aa)*abs(g2)+6*Math.abs(C)*bb*bb+2*Math.abs(a0)],
      [sub(mul(mul(2,alpha),mul(add(12,mul(B,alpha)),g3)),add(neg(mul(mul(alpha,add(mul(A,beta),a2)),g2)),add(add(mul(2*C,b3),mul(2*a0,beta)),2*a00))),
       2*aa*(12+Math.abs(B)*aa)*abs(g3)+aa*(Math.abs(A)*bb+Math.abs(a2))*abs(g2)+2*Math.abs(C)*bb**3+2*Math.abs(a0)*bb+2*Math.abs(a00)]
    ];
    return rows.map(([v,norm])=>abs(v)/Math.max(norm,1e-250));
  }
  function solve(c,free={beta:0,g2:.1,g3:.001}){
    const {A,B,C,a2,a0,a00}=c,points=[],messages=[],families=[],outside=[];
    const roots=quadratic(C,6*A+4*B,120);
    if(!roots.length)messages.push('No nonzero α solves the double-pole balance for this core.');
    for(const {value:alpha,id} of roots){
      const dB=add(2*A,mul(C,alpha)),d2=add(36,mul(A+2*B,alpha)),h3=add(12,mul(B,alpha)),d3=mul(mul(2,alpha),h3);
      const zb=zero(dB,2*Math.abs(A)+Math.abs(C)*abs(alpha)),z2=zero(d2,36+Math.abs(A+2*B)*abs(alpha)),z3=zero(h3,12+Math.abs(B)*abs(alpha));
      let betas,betaFree=false;
      if(!zb)betas=[div(-2*a2,dB)];
      else if(a2!==0){messages.push(`α = ${complexText(alpha)} requires a₂ = 0 in this sector.`);continue;}
      else if(z2&&C!==0)betas=quadratic(6*C,0,2*a0).map(x=>x.value);
      else if(z3&&!z2){
        const a=sub(2*C,div(mul(6*A*C,alpha),d2)),b=sub(2*a0,div(mul(2*A*a0,alpha),d2));
        if(!real(a)||!real(b)){messages.push('A complex resonant compatibility condition is unresolved in this real-coordinate chart.');continue;}
        betas=cubicRealLattice(a[0],b[0],2*a00);
      }else{betas=[[free.beta,0]];betaFree=true;}
      for(let j=0;j<betas.length;j++){
        const beta=betas[j],freeCoords=betaFree?['beta']:[],rhs2=add(mul(6*C,pow(beta,2)),2*a0);
        let G2;
        if(!z2)G2=div(rhs2,d2);
        else if(!zero(rhs2,6*Math.abs(C)*abs(beta)**2+2*Math.abs(a0)))continue;
        else{G2=[free.g2,0];freeCoords.push('g2');}
        const linear=mul(alpha,add(mul(A,beta),a2)),constant=add(add(mul(2*C,pow(beta,3)),mul(2*a0,beta)),2*a00);
        let G3;
        if(!z3)G3=div(sub(constant,mul(linear,G2)),d3);
        else{
          const remainder=sub(constant,mul(linear,G2)),scale=2*Math.abs(C)*abs(beta)**3+2*Math.abs(a0)*abs(beta)+2*Math.abs(a00)+abs(linear)*abs(G2);
          if(!zero(remainder,scale)){
            if(freeCoords.includes('g2')&&!zero(linear,abs(alpha)*(Math.abs(A)*abs(beta)+Math.abs(a2)))){G2=div(constant,linear);freeCoords.splice(freeCoords.indexOf('g2'),1);}
            else{messages.push(`α = ${complexText(alpha)} fails the constant-term compatibility condition.`);continue;}
          }
          G3=[free.g3,0];freeCoords.push('g3');
        }
        const branch=id*10+j;
        if(freeCoords.length)families.push({branch,alpha,beta,free:freeCoords.slice()});
        if(!real(G2)||!real(G3)){outside.push(branch);continue;}
        const g2=G2[0],g3=G3[0];
        if(![...alpha,...beta,g2,g3].every(Number.isFinite)||Math.abs(g2)>1e60||Math.abs(g3)>1e90){messages.push('A branch exceeds the numerical plotting range.');continue;}
        const errors=constraints(c,alpha,beta,[g2,0],[g3,0]);
        if(errors.some(e=>!Number.isFinite(e)||e>2e-8)){messages.push(`α = ${complexText(alpha)} has an unresolved near-singular coefficient condition.`);continue;}
        const geometry=W.geometry(g2,g3),terms=[[1,0,alpha[0]],[0,0,beta[0]]],imag_terms=[[1,0,alpha[1]],[0,0,beta[1]]];
        points.push({g2,g3,delta:g2*g2*g2-27*g3*g3,alpha,beta,branch,lead:id,character:1,free:freeCoords,operator:{2:a2,0:a0},constant:a00,
          profile:{kind:'polynomial',terms,imag_terms,scale:1,variable:'v',geometry},residual:Math.max(...errors)});
      }
    }
    return {points,families,outside,messages:[...new Set(messages)],coefficients:{...c}};
  }
  function getValue(state,key=state.sweep){return coefKeys.includes(key)?state.coeff[key]:state.free[key];}
  function setValue(state,key,value){if(coefKeys.includes(key))state.coeff[key]=value;else state.free[key]=value;}
  function defaultRange(state,key){
    const v=getValue(state,key),core=presets[presetKey(state.coeff)];
    if(key==='beta')return [-.21,.21];if(key==='g2')return [-.2,.2];if(key==='g3')return [-.03,.03];
    if(key==='a0')return core?core.range.slice():[v-2*Math.max(1,Math.abs(v)),v+2*Math.max(1,Math.abs(v))];
    if(key==='a2')return [-2,2];if(key==='a00')return [-.5,.5];
    const span=Math.max(Math.abs(v),10);return [v-span,v+span];
  }
  function build(state,count=481){
    const [lo,hi]=state.range,selected=getValue(state),values=Array.from({length:count},(_,j)=>lo+(hi-lo)*j/(count-1));
    values.push(selected);if(lo<=0&&hi>=0)values.push(0);
    if(state.sweep==='beta'&&presetKey(state.coeff)==='lax'&&state.coeff.a2===0&&state.coeff.a0===-1&&state.coeff.a00===0)values.push(-1/6,1/6);
    const parameters=[...new Set(values.filter(t=>Number.isFinite(t)&&t>=lo&&t<=hi))].sort((a,b)=>a-b),frames=[],ids=new Set();let previous=[],nextId=100;
    for(const t of parameters){
      const c={...state.coeff},f={...state.free};(coefKeys.includes(state.sweep)?c:f)[state.sweep]=t;const result=solve(c,f),pairs=[],usedOld=new Set(),usedNew=new Set();
      for(let j=0;j<result.points.length;j++)for(let k=0;k<previous.length;k++)if(result.points[j].lead===previous[k].lead)pairs.push({j,k,d:abs(sub(result.points[j].beta,previous[k].beta))/(1+abs(previous[k].beta)+abs(result.points[j].beta))});
      pairs.sort((a,b)=>a.d-b.d);
      for(const {j,k} of pairs)if(!usedOld.has(k)&&!usedNew.has(j)){result.points[j].branch=previous[k].branch;usedOld.add(k);usedNew.add(j);}
      for(let j=0;j<result.points.length;j++)if(!usedNew.has(j)){if(ids.has(result.points[j].branch))result.points[j].branch=nextId++;}
      result.points.forEach(p=>ids.add(p.branch));frames.push({t,...result});previous=result.points;
    }
    const tracks=[...ids].sort((a,b)=>a-b).map(branch=>({branch,character:1,breaks:[],points:frames.map(f=>{const p=f.points.find(p=>p.branch===branch);return p?[p.g2,p.g3,f.t]:null;})}));
    const selectedIndex=parameters.reduce((best,t,j)=>Math.abs(t-selected)<Math.abs(parameters[best]-selected)?j:best,0);
    const model={slug:'c4',kind:'c4',n:null,p:4,q:2,K:0,primary:false,variable:'v',display_scale:1,title:'C4 · coefficient explorer',parameter:labels[state.sweep],interval:[lo,hi],frames,tracks,selectedIndex,specials:[],count_kind:'representatives',notes:[],figure:'c4',focus:state.focus};
    model.limits=limits(model,null);return model;
  }
  function limits(model,branch,point){
    let tracks=model.focus&&branch!==null?model.tracks.filter(t=>t.branch===branch):model.tracks;
    let points=tracks.flatMap(t=>t.points.filter(Boolean));if(!points.length)points=model.frames.flatMap(f=>f.points.map(p=>[p.g2,p.g3]));
    if(!points.length)return [[-.1,.1],[-.015,.015]];
    return [0,1].map(axis=>{const a=points.map(p=>p[axis]).filter(Number.isFinite).sort((a,b)=>a-b);let lo=a[Math.floor(a.length*.01)],hi=a[Math.min(a.length-1,Math.floor(a.length*.99))];if(point){lo=Math.min(lo,point[axis===0?'g2':'g3']);hi=Math.max(hi,point[axis===0?'g2':'g3']);}let span=hi-lo;if(span<1e-10*Math.max(Math.abs(hi),Math.abs(lo),1e-30)){span=Math.max(Math.abs(hi)*.5,axis===0?.01:.001);lo-=span/2;hi+=span/2;}return [lo-span*.08,hi+span*.08];});
  }
  function equation(c){let out='v⁗';for(const [key,term] of [['A','v v″'],['B','(v′)²'],['C','v³'],['a2','v″'],['a0','v'],['a00','']]){const x=c[key];if(x===0)continue;out+=' '+(x<0?'−':'+')+' '+(Math.abs(x)===1&&term?'':format(Math.abs(x))+' ')+term;}return out.trim()+' = 0';}
  function fromParams(params){const state=oneGap();const key=params.get('preset');if(hasPreset(key))Object.assign(state,preset(key));for(const k of [...coefKeys,...freeKeys])if(params.has(k)){const v=parseNumber(params.get(k));if(v!==null)setValue(state,k,v);}if([...coefKeys,...freeKeys].includes(params.get('sweep')))state.sweep=params.get('sweep');state.range=defaultRange(state,state.sweep);const lo=parseNumber(params.get('min')),hi=parseNumber(params.get('max'));if(lo!==null&&hi!==null&&lo<hi)state.range=[lo,hi];if(params.has('parameter')){const v=parseNumber(params.get('parameter'));if(v!==null)setValue(state,state.sweep,v);}const v=getValue(state);state.range=[Math.min(state.range[0],v),Math.max(state.range[1],v)];if(params.has('focus'))state.focus=params.get('focus')==='1';return state;}
  function toParams(params,state){for(const k of [...coefKeys,...freeKeys])params.set(k,String(getValue(state,k)));params.set('sweep',state.sweep);params.set('min',String(state.range[0]));params.set('max',String(state.range[1]));params.set('focus',state.focus?'1':'0');}
  function csv(model){const rows=[['parameter','alpha_re','alpha_im','beta_re','beta_im','g2','g3','A','B','C','a2','a0','a00']];for(const f of model.frames)for(const p of f.points)rows.push([f.t,...p.alpha,...p.beta,p.g2,p.g3,...coefKeys.map(k=>f.coefficients[k])]);return rows.map(r=>r.join(',')).join('\n')+'\n';}
  root.C4={presets,preset,presetKey,oneGap,labels,coefKeys,freeKeys,parseNumber,format,complexText,solve,build,limits,equation,fromParams,toParams,getValue,setValue,defaultRange,csv,constraints};
  if(typeof module!=='undefined')module.exports=root.C4;
})(typeof globalThis!=='undefined'?globalThis:window);
