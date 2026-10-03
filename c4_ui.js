/* Coefficient controls for the explicitly constructed C4 affine-wp sector. */
(function(root){
  'use strict';
  let state=C4.oneGap(),index=-1;
  const el=id=>document.getElementById(id),clone=x=>JSON.parse(JSON.stringify(x));
  function inputText(value){
    const plain=String(value);if(plain.length<12||value===0)return plain;
    for(let q=2;q<=120;q++){const n=Math.round(value*q);if(Math.abs(value-n/q)<1e-14*Math.abs(value))return n+'/'+q;}
    return plain;
  }
  function available(){return new Set(C4.solve(state.coeff,state.free).families.flatMap(f=>f.free));}
  function controls(){
    el('c4-preset').value=C4.presetKey(state.coeff);el('c4-sweep').value=state.sweep;
    for(const key of [...C4.coefKeys,...C4.freeKeys])el('c4-'+key).value=inputText(C4.getValue(state,key));
    el('c4-min').value=inputText(state.range[0]);el('c4-max').value=inputText(state.range[1]);el('c4-focus').checked=state.focus;
    const free=available();
    for(const key of C4.freeKeys){el('c4-field-'+key).hidden=!free.has(key);el('c4-sweep').querySelector(`option[value="${key}"]`).disabled=!free.has(key);}
    el('c4-free-title').hidden=!free.size;
  }
  function rebuild(resetBranch=false){
    const free=available();
    if(C4.freeKeys.includes(state.sweep)&&!free.has(state.sweep)){state.sweep=free.values().next().value||'a0';state.range=C4.defaultRange(state,state.sweep);}
    const value=C4.getValue(state);state.range=[Math.min(state.range[0],value),Math.max(state.range[1],value)];
    if(state.range[0]===state.range[1])state.range=[value-1,value+1];
    ATLAS[index]=C4.build(state);const m=ATLAS[index];
    if(resetBranch)selectedBranch=null;
    controls();el('paramLabel').textContent=m.parameter;slider.max=m.frames.length-1;slider.value=m.selectedIndex;
    el('c4-error').textContent='';update();
  }
  function setPreset(key){state=C4.preset(key);preferredProfileMode='regular';rebuild(true);}
  function oneGap(pulse=false){state=C4.oneGap();if(pulse)state.free.beta=1/6;preferredProfileMode='regular';selectedBranch=0;rebuild(false);}
  function register(atlas){
    index=atlas.length;atlas.push(C4.build(state));
    el('c4-preset').innerHTML=Object.entries(C4.presets).map(([key,p])=>`<option value="${key}">${p.name}</option>`).join('')+'<option value="custom">Custom coefficients</option>';
    el('c4-sweep').innerHTML='<optgroup label="Equation coefficients">'+C4.coefKeys.map(k=>`<option value="${k}">${C4.labels[k]}</option>`).join('')+'</optgroup><optgroup label="Free wave coordinates">'+C4.freeKeys.map(k=>`<option value="${k}">${C4.labels[k]}</option>`).join('')+'</optgroup>';
    el('c4-preset').onchange=()=>{if(C4.presets[el('c4-preset').value])setPreset(el('c4-preset').value);};
    el('c4-sweep').onchange=()=>{state.sweep=el('c4-sweep').value;state.range=C4.defaultRange(state,state.sweep);rebuild();};
    for(const key of [...C4.coefKeys,...C4.freeKeys])el('c4-'+key).onchange=()=>{
      const value=C4.parseNumber(el('c4-'+key).value);
      if(value===null){el('c4-error').textContent='Enter a finite real number or fraction, with magnitude at most 10⁶.';return;}
      C4.setValue(state,key,value);rebuild();
    };
    for(const id of ['c4-min','c4-max'])el(id).onchange=()=>{
      const lo=C4.parseNumber(el('c4-min').value),hi=C4.parseNumber(el('c4-max').value);
      if(lo===null||hi===null||lo>=hi){el('c4-error').textContent='The sweep minimum must be smaller than its maximum.';return;}
      state.range=[lo,hi];C4.setValue(state,state.sweep,Math.max(lo,Math.min(hi,C4.getValue(state))));rebuild();
    };
    el('c4-focus').onchange=()=>{state.focus=el('c4-focus').checked;ATLAS[index].focus=state.focus;update();};
    controls();
  }
  function show(){el('c4-controls').hidden=false;controls();}
  function hide(){el('c4-controls').hidden=true;}
  function onFrame(model,frame,point){
    C4.setValue(state,state.sweep,frame.t);el('c4-'+state.sweep).value=inputText(frame.t);
    el('c4-preset').value=C4.presetKey(frame.coefficients);model.limits=C4.limits(model,point?point.branch:null,point);
    const free=new Set(frame.families.flatMap(f=>f.free));
    for(const key of C4.freeKeys){el('c4-field-'+key).hidden=!free.has(key);el('c4-sweep').querySelector(`option[value="${key}"]`).disabled=!free.has(key);}
    el('c4-free-title').hidden=!free.size;
  }
  function count(frame){
    const n=frame.points.length,free=new Set(frame.families.flatMap(f=>f.free));
    if(!n)return free.size?'Continuous affine-℘ family; no real lattice at the selected free coordinates':'No real-lattice affine-℘ representative at these coefficients';
    return `${n} displayed affine-℘ representative${n===1?'':'s'}`+(free.size?' · continuous family present':'');
  }
  function pointLabel(point){return `α = ${C4.complexText(point.alpha)}, β = ${C4.complexText(point.beta)}`;}
  function renderNotes(frame){
    const free=[...new Set(frame.families.flatMap(f=>f.free))];
    const lines=['C4: v⁗ + A v v″ + B (v′)² + C v³ + a₂ v″ + a₀ v + a₀₀ = 0.',
      'This constructor shows v = α℘ + β. Multiple-pole elliptic and genus-two solutions are outside this selected sector.'];
    if(free.length)lines.push('Free wave coordinates: '+free.map(k=>C4.labels[k].split(' · ')[0]).join(', ')+'. Markers show the chosen members; a free coordinate sweep changes the wave while holding the equation fixed.');
    if(frame.outside.length)lines.push(`${frame.outside.length} constructed branch${frame.outside.length===1?' has':'es have'} complex lattice invariants at these free coordinates and lies outside this real plane.`);
    lines.push(...frame.messages);
    if(C4.presetKey(frame.coefficients)==='cubic'&&frame.coefficients.a00===0)lines.push('The two opposite-sign affine profiles use the same lattice and are listed separately. The twisted cubic sector is available in the pure-power atlas.');
    el('notes').replaceChildren(...lines.map(text=>{const p=document.createElement('p');p.textContent=text;return p;}));
    const p=document.createElement('p'),a=document.createElement('a');a.href='data/c4_formulas.json';a.textContent='C4 coefficient identities ↗';p.append(a,document.createTextNode(' · '));
    const button=document.createElement('button');button.className='text-button';button.textContent='Download this sweep as CSV';
    button.onclick=()=>{const url=URL.createObjectURL(new Blob([C4.csv(ATLAS[index])],{type:'text/csv;charset=utf-8'})),link=document.createElement('a');link.href=url;link.download='c4-'+state.sweep+'.csv';link.click();setTimeout(()=>URL.revokeObjectURL(url),1000);};p.append(button);el('notes').append(p);
  }
  function shortcuts(){
    el('shortcuts').innerHTML='<button id="c4-one-gap">Lax one-gap curve</button><button id="c4-pulse">Lax sech² pulse · β = 1/6</button><button id="c4-two-branches">Kaup–Kupershmidt: two branches</button>';
    el('c4-one-gap').onclick=()=>oneGap(false);el('c4-pulse').onclick=()=>oneGap(true);el('c4-two-branches').onclick=()=>setPreset('kk');
  }
  function restore(params){state=C4.fromParams(params);const free=available();if(C4.freeKeys.includes(state.sweep)&&!free.has(state.sweep)){state.sweep=free.values().next().value||'a0';state.range=C4.defaultRange(state,state.sweep);}ATLAS[index]=C4.build(state);}
  root.C4Panel={register,show,hide,onFrame,count,pointLabel,renderNotes,shortcuts,restore,setPreset,oneGap,rebuild,
    toParams:p=>C4.toParams(p,state),state:()=>clone(state),index:()=>index};
})(typeof globalThis!=='undefined'?globalThis:window);
