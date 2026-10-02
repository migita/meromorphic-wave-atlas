/* Real Weierstrass slices and analytic character profiles.
 * A finite real lattice need not carry a smooth real-valued solution.
 * Zero/pole continuation is explicit; principal square roots are not glued
 * across a zero. No numerical time integration is used to draw these curves.
 */
(function(root){
  'use strict';
  const PI=Math.PI,mod=(x,p)=>((x%p)+p)%p;
  function jacobi(u,m,mc=1-m){
    if(m===0)return {sn:Math.sin(u),cn:Math.cos(u),dn:1};
    if(mc===0){const cn=1/Math.cosh(u);return {sn:Math.tanh(u),cn,dn:cn};}
    const a=[1],c=[Math.sqrt(Math.max(0,m))];let b=Math.sqrt(mc),i=0,two=1;
    while(i<18&&Math.abs(c[i]/a[i])>2e-15){const old=a[i];c.push((old-b)/2);a.push((old+b)/2);b=Math.sqrt(old*b);i++;two*=2;}
    const K=PI/(2*a[i]);u=mod(u+2*K,4*K)-2*K;
    let phi=two*a[i]*u;
    while(i>0){const t=c[i]*Math.sin(phi)/a[i];phi=(phi+Math.asin(Math.max(-1,Math.min(1,t))))/2;i--;}
    const sn=Math.sin(phi),cn=Math.cos(phi),dn=Math.sqrt(cn*cn+mc*sn*sn);
    return {sn,cn,dn};
  }
  function regularAvailable(g){return g.type==='rect';}
  function actualMode(spec,preferred){return preferred==='regular'&&regularAvailable(spec.geometry)?'regular':'axis';}
  function cellScale(g){return g.period||8/Math.sqrt(g.G||g.H||1);}
  function nearPole(g,x,mode,tolerance=1e-9){
    if(mode==='regular')return false;
    if(g.period){const t=mod(x+g.period/2,g.period)-g.period/2;return Math.abs(t)<tolerance*g.period;}
    return Math.abs(x)<tolerance*cellScale(g);
  }
  function base(g,x,mode){
    if(g.type==='rational')return {X:1/(x*x),Y:-2/(x*x*x),sn:0,cn:0,dn:0};
    let u=Math.sqrt(g.G||g.H)*x;
    if(g.type==='one')u*=2;
    const j=jacobi(u,g.m,g.mc),{sn,cn,dn}=j;
    if(g.type==='one'){
      const den=1-cn,H=g.H;
      return {...j,X:g.e+H*(1+cn)/den,Y:-4*Math.pow(H,1.5)*sn*dn/(den*den)};
    }
    if(mode==='regular')return {...j,X:g.e[2]+g.D*sn*sn,Y:2*g.D*Math.sqrt(g.G)*sn*cn*dn};
    return {...j,X:g.e[2]+g.G/(sn*sn),Y:-2*Math.pow(g.G,1.5)*cn*dn/(sn*sn*sn)};
  }
  function polynomial(terms,X,Y){let result=0;for(const [i,j,c] of terms)result+=c*Math.pow(X,i)*Math.pow(Y,j);return result;}
  function geometry(g2,g3){
    if(!Number.isFinite(g2)||!Number.isFinite(g3))throw Error('Nonfinite lattice invariants');
    if(g2===0&&g3===0)return {type:'rational',delta:0,period:null,scale:1};
    const scale=Math.max(Math.sqrt(Math.abs(g2)),Math.cbrt(Math.abs(g3))),a=g2/scale/scale,b=g3/scale/scale/scale;
    const disc=a*a*a-27*b*b,den=Math.abs(a*a*a)+27*b*b;
    const ellipticK=mc=>{let x=1,y=Math.sqrt(Math.max(0,mc));for(let j=0;j<20&&Math.abs(x-y)>2e-15*x;j++){const next=(x+y)/2;y=Math.sqrt(x*y);x=next;}return Math.PI/(2*x);};
    if(Math.abs(disc)<=2e-12*den&&g2>=0){
      const e=Math.sign(-g3)*Math.sqrt(g2/12);
      if(e>0)return {type:'rect',delta:0,node:true,e:[e,e,-2*e],G:3*e,D:3*e,m:1,mc:0,K:null,period:null};
      const G=-3*e;return {type:'rect',delta:0,node:true,e:[-2*e,e,e],G,D:0,m:0,mc:1,K:PI/2,period:PI/Math.sqrt(G)};
    }
    if(disc>0){
      const r=Math.sqrt(a/12),theta=Math.acos(Math.max(-1,Math.min(1,3*Math.sqrt(3)*b/Math.pow(a,1.5))))/3;
      const e=[0,1,-1].map(j=>2*r*Math.cos(theta-j*2*PI/3)*scale),G=e[0]-e[2],D=e[1]-e[2],mc=(e[0]-e[1])/G,m=D/G,K=ellipticK(mc);
      return {type:'rect',delta:1,e,G,D,m,mc,K,period:2*K/Math.sqrt(G)};
    }
    const rad=Math.sqrt(Math.max(0,b*b/64-a*a*a/1728)),en=Math.cbrt(b/8+rad)+Math.cbrt(b/8-rad),hn=Math.sqrt(3*en*en-a/4);
    const e=en*scale,H=hn*scale,m=.5-3*en/(4*hn),mc=.5+3*en/(4*hn),K=ellipticK(mc);
    return {type:'one',delta:-1,e,H,m,mc,K,period:2*K/Math.sqrt(H)};
  }
  function halfroot(spec,x,mode,basis){
    const g=spec.geometry;if(g.type==='rational')return {v:1/x,d:-1/(x*x),imag:false};
    const {sn,cn,dn}=basis;
    if(g.type==='one')return {v:Math.sqrt(g.H)*sn/(1-cn),d:-2*g.H*dn/(1-cn),imag:false};
    const G=g.G,s=Math.sqrt(G),D=Math.sqrt(Math.max(0,g.D)),i=spec.root_index;
    if(mode==='regular'){
      if(i===2)return {v:D*sn,d:D*s*cn*dn,imag:false};
      if(i===1)return {v:D*cn,d:-D*s*sn*dn,imag:true};
      return {v:s*dn,d:-g.D*sn*cn,imag:true};
    }
    if(i===2)return {v:s/sn,d:-G*cn*dn/(sn*sn),imag:false};
    if(i===1)return {v:s*dn/sn,d:-G*cn/(sn*sn),imag:false};
    return {v:s*cn/sn,d:-G*dn/(sn*sn),imag:false};
  }
  function crossings(spec,x,mode){
    const g=spec.geometry,zero=spec.zeros[mode];let count=0;
    if(g.period){
      if(mode==='axis')count+=Math.floor(x/g.period);
      if(zero)count+=Math.floor(x/g.period-zero.fraction);
    }else{
      if(mode==='axis')count+=(x>0?1:0);
      if(zero)count+=(x>zero.position?1:0);
    }
    return count;
  }
  function rootValue(spec,x,mode,basis){
    const {X,Y}=basis,n=spec.order;
    const F=polynomial(spec.terms,X,Y),other=polynomial(spec.terms,X,-Y);
    let magnitude=Math.pow(Math.abs(F),1/n);
    if(mode==='regular'&&spec.geometry.node&&spec.geometry.D===0){
      const gap=X-spec.point[0],scale=Math.max(Math.abs(X),Math.abs(spec.point[0]),1e-30);
      if(Math.abs(gap)<scale*1e-9)return [0,0];
      magnitude=Math.sqrt(Math.abs(gap));
      const phase=spec.power_sign[mode]<0?(n%2?PI:PI/n):0;
      return [magnitude*Math.cos(phase),magnitude*Math.sin(phase)];
    }
    // div(F)=nQ-nO implies F(X,Y)F(X,-Y)=(-1)^n(X-A)^n.
    // Use its noncancelling side at a high-order zero of F.
    if(Math.abs(other)>Math.abs(F)&&Math.abs(other)>0){
      if(n%2)return [-(X-spec.point[0])/(Math.sign(other)*Math.pow(Math.abs(other),1/n)),0];
      magnitude=Math.abs(X-spec.point[0])/Math.pow(Math.abs(other),1/n);
    }
    if(n%2)return [Math.sign(F)*magnitude,0];
    const reference=cellScale(spec.geometry)*.071;
    const parity=crossings(spec,x,mode)-crossings(spec,reference,mode);
    magnitude*=Math.abs(parity%2)===1?-1:1;
    const phase=spec.power_sign[mode]<0?PI/n:0;
    return [magnitude*Math.cos(phase),magnitude*Math.sin(phase)];
  }
  function at(spec,x,preferred='regular'){
    const mode=actualMode(spec,preferred);
    if(nearPole(spec.geometry,x,mode))return null;
    const b=base(spec.geometry,x,mode);let re=0,im=0;
    if(spec.kind==='polynomial'){re=polynomial(spec.terms,b.X,b.Y);if(spec.imag_terms)im=polynomial(spec.imag_terms,b.X,b.Y);}
    else if(spec.kind==='halfroot'){
      const v=halfroot(spec,x,mode,b);
      let value=spec.form==='generator'?v.v:spec.form==='derivative'?-v.d+spec.b*v.v:(b.X+spec.d)*v.v;
      if(v.imag)im=value;else re=value;
    }else [re,im]=rootValue(spec,x,mode,b);
    re*=spec.scale;im*=spec.scale;
    return Number.isFinite(re)&&Number.isFinite(im)?[re,im]:null;
  }
  function periodMultiplier(spec,mode){
    if(spec.kind==='polynomial')return 1;
    if(spec.kind==='halfroot')return spec.geometry.type==='rect'&&spec.root_index!==0?2:1;
    if(spec.order%2)return 1;
    return ((spec.zeros[mode]?1:0)+(mode==='axis'?1:0))%2?2:1;
  }
  function generate(spec,preferred='regular',cells=2,count=801){
    const g=spec.geometry,mode=actualMode(spec,preferred),basePeriod=g.period;
    const width=basePeriod?basePeriod*cells:2*cellScale(g);
    const x=[],re=[],im=[],poles=[];let anyImag=false;
    if(mode==='axis'){
      if(basePeriod){for(let k=Math.ceil(-width/(2*basePeriod));k<=Math.floor(width/(2*basePeriod));k++)poles.push(k*basePeriod);}
      else poles.push(0);
    }
    for(let j=0;j<count;j++){
      const xx=-width/2+width*j/(count-1);x.push(xx);
      const v=nearPole(g,xx,mode,.0015)?null:at(spec,xx,mode);
      re.push(v?v[0]:null);im.push(v?v[1]:null);
    }
    const realValues=re.filter(Number.isFinite),imagValues=im.filter(Number.isFinite);
    const maximum=Math.max(...realValues.map(Math.abs),...imagValues.map(Math.abs),1e-300);
    anyImag=Math.max(...imagValues.map(Math.abs),0)>maximum*1e-9;
    const all=anyImag?realValues.concat(imagValues):realValues;
    let low=Math.min(...all),high=Math.max(...all),clipped=false;
    if(mode==='axis'&&poles.length){
      const sorted=all.slice().sort((a,b)=>a-b),q=a=>sorted[Math.floor(a*(sorted.length-1))];
      low=Math.min(0,q(.08));high=Math.max(0,q(.92));
      const span=high-low||Math.max(Math.abs(high),1e-6);low-=span*.12;high+=span*.12;clipped=true;
    }
    let span=high-low;
    const variation=Math.max(Math.max(...realValues)-Math.min(...realValues),Math.max(...imagValues)-Math.min(...imagValues));
    const constant=mode==='regular'&&variation<1e-9*Math.max(Math.abs(high),1e-20);
    if(span<1e-12*Math.max(Math.abs(high),1e-20)){span=Math.max(Math.abs(high)*.2,1e-4);low-=span/2;high+=span/2;}
    else{low-=span*.1;high+=span*.1;}
    let period=basePeriod?basePeriod*periodMultiplier(spec,mode):null,mean=null,amplitude=null;
    if(mode==='regular'&&!anyImag){
      amplitude=Math.max(...realValues)-Math.min(...realValues);
      if(basePeriod){
        const n=1024;let sum=0;
        for(let j=0;j<n;j++)sum+=at(spec,(j+.37)*period/n,mode)[0];
        mean=sum/n;
      }
      if(constant){mean=realValues[0];amplitude=0;period=null;}
    }
    if(constant)period=null;
    return {x,re,im,mode,anyImag,poles,clipped,constant,period,mean,amplitude,range:[low,high],width};
  }
  root.WaveMath={jacobi,base,polynomial,geometry,halfroot,at,generate,actualMode,periodMultiplier};
  if(typeof module!=='undefined')module.exports=root.WaveMath;
})(typeof globalThis!=='undefined'?globalThis:window);
