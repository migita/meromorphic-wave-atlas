// Optional offline gallery smoke check. Node >= 22 and a Chromium binary.
// Set ATLAS_CHROME when Chromium is not in the usual Playwright cache.
import fs from 'node:fs';
import os from 'node:os';
import path from 'node:path';
import {spawn} from 'node:child_process';
import {fileURLToPath,pathToFileURL} from 'node:url';

const here=path.dirname(fileURLToPath(import.meta.url));
const cache=path.join(os.homedir(),'.cache','ms-playwright');
const candidates=fs.existsSync(cache)?fs.readdirSync(cache).filter(n=>n.startsWith('chromium-')).map(n=>path.join(cache,n,'chrome-linux64','chrome')):[];
const chrome=process.env.ATLAS_CHROME||candidates.find(p=>fs.existsSync(p));
if(!chrome)throw new Error('Set ATLAS_CHROME to a Chromium executable.');
const profile=fs.mkdtempSync(path.join(os.tmpdir(),'wave-atlas-browser-'));
const proc=spawn(chrome,['--headless','--no-sandbox','--disable-gpu','--no-first-run','--remote-debugging-port=0',`--user-data-dir=${profile}`],{stdio:['ignore','ignore','pipe']});
let socket;
try{
  const endpoint=await new Promise((resolve,reject)=>{
    let text='';const timer=setTimeout(()=>reject(new Error('Browser startup timeout')),20000);
    proc.stderr.on('data',buf=>{text+=buf;const m=text.match(/DevTools listening on (ws:\/\/[^\s]+)/);if(m){clearTimeout(timer);resolve(m[1]);}});
    proc.on('error',reject);
  });
  socket=new WebSocket(endpoint);await new Promise((resolve,reject)=>{socket.onopen=resolve;socket.onerror=reject;});
  let next=0;const pending=new Map(),errors=[];
  socket.onmessage=event=>{const msg=JSON.parse(event.data);if(msg.id){const request=pending.get(msg.id);pending.delete(msg.id);msg.error?request.reject(msg.error):request.resolve(msg.result);}else if(msg.method==='Runtime.exceptionThrown')errors.push(msg.params);};
  const send=(method,params={},sessionId)=>new Promise((resolve,reject)=>{const id=++next;pending.set(id,{resolve,reject});socket.send(JSON.stringify({id,method,params,...(sessionId?{sessionId}:{})}));});
  const {targetId}=await send('Target.createTarget',{url:'about:blank'});
  const {sessionId}=await send('Target.attachToTarget',{targetId,flatten:true});
  const call=(method,params={})=>send(method,params,sessionId);
  const evaluate=async expression=>{const r=await call('Runtime.evaluate',{expression,returnByValue:true,awaitPromise:true});if(r.exceptionDetails)throw new Error(JSON.stringify(r.exceptionDetails));return r.result.value;};
  await call('Runtime.enable');await call('Page.enable');
  await call('Emulation.setDeviceMetricsOverride',{width:1440,height:1120,deviceScaleFactor:1,mobile:false});
  const pageURL=process.env.ATLAS_URL||pathToFileURL(path.join(here,'index.html')).href;
  await call('Page.navigate',{url:pageURL});
  for(let tries=0;tries<300;tries++){
    if(await evaluate('document.readyState === "complete" && typeof ATLAS !== "undefined"'))break;
    await new Promise(resolve=>setTimeout(resolve,100));
  }
  const counts=await evaluate('({slices:ATLAS.length,pairs:ATLAS.filter(m=>m.primary).length,cards:document.querySelectorAll(".card").length})');
  if(counts.slices!==19||counts.pairs!==13||counts.cards!==19)throw new Error(JSON.stringify(counts));
  const checked=[];
  for(let i=0;i<counts.slices;i++){
    const result=await evaluate(`(()=>{choose(${i});const m=ATLAS[${i}],out=[];for(const f of [0,.5,1]){slider.value=Math.floor((m.frames.length-1)*f);slider.dispatchEvent(new Event('input'));const points=m.frames[+slider.value].points;if(points.length&&(!lastProfile||!lastProfile.range.every(Number.isFinite)))throw Error('Missing wave profile');out.push({count:document.getElementById('count').textContent,profile:lastProfile?{mode:lastProfile.mode,constant:lastProfile.constant,period:lastProfile.period}:null});}slider.value=Math.floor(m.frames.length/2);update();const points=m.frames[+slider.value].points;for(const pt of points){selectedBranch=pt.branch;update();if(points[selected].branch!==pt.branch||!lastProfile)throw Error('Wrong selected branch');}return {slug:m.slug,readouts:out,branches:points.length};})()`);
    checked.push(result);
  }
  const kawahara=await evaluate(`(()=>{choose(ATLAS.findIndex(m=>m.slug==='n2_p4_kawahara'));const m=ATLAS[current];return [-13/6,13/6].map(T=>{let index=0;for(let i=1;i<m.frames.length;i++)if(Math.abs(m.frames[i].t-T)<Math.abs(m.frames[index].t-T))index=i;slider.value=index;update();return {parameter:m.frames[index].t,count:m.frames[index].points.length,delta:m.frames[index].points[0]?.delta};});})()`);
  if(kawahara.some(r=>r.count!==1||Math.abs(r.delta)>1e-19))throw new Error('Kawahara coalescences: '+JSON.stringify(kawahara));
  const limits=await evaluate(`(()=>{preferredProfileMode='regular';toParameter(-13/6);const pulse=lastProfile;const middle=pulse.re[Math.floor(pulse.re.length/2)];if(Math.abs(middle+11/12)>1e-10||pulse.constant||pulse.period!==null)throw Error('Incorrect pulse');toParameter(13/6);if(!lastProfile.constant||Math.abs(lastProfile.mean-2)>1e-10)throw Error('Incorrect upper regular limit');preferredProfileMode='axis';drawWave();if(lastProfile.constant||lastProfile.poles.length<1)throw Error('Missing trigonometric poles');preferredProfileMode='regular';return {pulse_minimum:middle,upper_regular_mean:2,trigonometric_poles:lastProfile.poles.length};})()`);
  const links=await evaluate('Array.from(document.querySelectorAll("a[href]")).map(a=>a.getAttribute("href"))');
  const missing=Array.from(new Set(links)).filter(href=>!href.startsWith('#')&&!/^https?:/.test(href)).filter(href=>!fs.existsSync(path.resolve(here,href.split('#')[0])));
  if(missing.length)throw new Error('Missing local links: '+JSON.stringify(missing));
  await evaluate(`choose(ATLAS.findIndex(m=>m.slug==='n2_p4_kawahara'))`);
  const screenshots=path.join(here,'data','browser');fs.mkdirSync(screenshots,{recursive:true});
  const screen=await call('Page.captureScreenshot',{format:'png'});fs.writeFileSync(path.join(screenshots,'desktop.png'),Buffer.from(screen.data,'base64'));
  await call('Emulation.setDeviceMetricsOverride',{width:390,height:844,deviceScaleFactor:1,mobile:true});
  await evaluate('new Promise(resolve=>requestAnimationFrame(()=>requestAnimationFrame(resolve)))');
  const mobile=await evaluate('({width:innerWidth,document:document.documentElement.scrollWidth,cards:document.querySelectorAll(".card").length})');
  if(mobile.document>391||mobile.width>391)throw new Error('Mobile horizontal overflow: '+JSON.stringify(mobile));
  const mobileScreen=await call('Page.captureScreenshot',{format:'png'});fs.writeFileSync(path.join(screenshots,'mobile.png'),Buffer.from(mobileScreen.data,'base64'));
  if(errors.length)throw new Error('JavaScript errors: '+JSON.stringify(errors));
  const report={status:'PASS',public_url:process.env.ATLAS_URL||null,...counts,checked,kawahara,limits,mobile,missing_links:missing,javascript_errors:errors};
  const reportPath=process.env.ATLAS_URL?'live_site_check.json':'browser_check.json';
  fs.writeFileSync(path.join(here,'data',reportPath),JSON.stringify(report,null,2)+'\n');
  console.log('PASS: 19 families with solution plots, 57 slider positions, branch changes, pulse/constant/trigonometric limits, links, desktop and mobile.');
}finally{
  if(socket)socket.close();proc.kill();
  await new Promise(resolve=>proc.once('exit',resolve));
  fs.rmSync(profile,{recursive:true,force:true});
}
