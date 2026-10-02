// Optional offline interaction check; Node >=22 and local Chromium.
import fs from 'node:fs';
import os from 'node:os';
import path from 'node:path';
import {spawn} from 'node:child_process';
import {fileURLToPath,pathToFileURL} from 'node:url';
const here=path.dirname(fileURLToPath(import.meta.url));
const cache=path.join(os.homedir(),'.cache','ms-playwright');
const chrome=process.env.ATLAS_CHROME||fs.readdirSync(cache).filter(n=>n.startsWith('chromium-')).map(n=>path.join(cache,n,'chrome-linux64','chrome')).find(p=>fs.existsSync(p));
if(!chrome)throw new Error('Set ATLAS_CHROME to a Chromium executable.');
const temporary=fs.mkdtempSync(path.join(os.tmpdir(),'kawahara-physical-'));
const proc=spawn(chrome,['--headless','--no-sandbox','--use-gl=angle','--use-angle=swiftshader','--enable-unsafe-swiftshader',
  '--no-first-run','--remote-debugging-port=0',`--user-data-dir=${temporary}`],{stdio:['ignore','ignore','pipe']});
let socket;
try{
  const endpoint=await new Promise((resolve,reject)=>{let text='';const timer=setTimeout(()=>reject(new Error('Browser startup timeout')),20000);
    proc.stderr.on('data',buf=>{text+=buf;const match=text.match(/DevTools listening on (ws:\/\/[^\s]+)/);if(match){clearTimeout(timer);resolve(match[1]);}});proc.on('error',reject);});
  socket=new WebSocket(endpoint);await new Promise((resolve,reject)=>{socket.onopen=resolve;socket.onerror=reject;});
  let next=0;const pending=new Map(),errors=[];
  socket.onmessage=event=>{const message=JSON.parse(event.data);if(message.id){const p=pending.get(message.id);pending.delete(message.id);message.error?p.reject(message.error):p.resolve(message.result);}else if(message.method==='Runtime.exceptionThrown')errors.push(message.params);};
  const send=(method,params={},sessionId)=>new Promise((resolve,reject)=>{const id=++next;pending.set(id,{resolve,reject});socket.send(JSON.stringify({id,method,params,...(sessionId?{sessionId}:{})}));});
  const {targetId}=await send('Target.createTarget',{url:'about:blank'});
  const {sessionId}=await send('Target.attachToTarget',{targetId,flatten:true});
  const call=(method,params={})=>send(method,params,sessionId);
  const evaluate=async expression=>{const result=await call('Runtime.evaluate',{expression,returnByValue:true,awaitPromise:true});if(result.exceptionDetails)throw new Error(JSON.stringify(result.exceptionDetails));return result.result.value;};
  await call('Runtime.enable');await call('Page.enable');
  await call('Emulation.setDeviceMetricsOverride',{width:1480,height:1030,deviceScaleFactor:1,mobile:false});
  await call('Page.navigate',{url:pathToFileURL(path.join(here,'kawahara_physical.html')).href});
  for(let tries=0;tries<150;tries++){
    if(await evaluate('document.readyState==="complete" && !!document.getElementById("kawahara-physical")?._fullLayout?.scene?._scene'))break;
    await new Promise(resolve=>setTimeout(resolve,100));
  }
  const count=await evaluate('document.getElementById("kawahara-physical").layout.updatemenus[0].buttons.length');
  const cases=[];
  for(let i=0;i<count;i++){
    const result=await evaluate(`(async()=>{const graph=document.getElementById('kawahara-physical'),button=graph.layout.updatemenus[0].buttons[${i}];await Plotly.update(graph,...button.args);const y=Array.from(graph.data[2].y);return {label:button.label,wavelength:graph.data[1].x[0],amplitude:graph.data[1].y[0],mean:graph.data[1].z[0],meanLine:graph.data[3].y[0],sampledAmplitude:Math.max(...y)-Math.min(...y),profilePoints:y.length};})()`);
    if(result.profilePoints!==1601||Math.abs(result.mean-result.meanLine)>1e-12||Math.abs(result.amplitude-result.sampledAmplitude)>1e-3)throw new Error(JSON.stringify(result));
    cases.push(result);
  }
  await evaluate(`(async()=>{const graph=document.getElementById('kawahara-physical');await Plotly.update(graph,...graph.layout.updatemenus[0].buttons[5].args);})()`);
  await evaluate('new Promise(resolve=>requestAnimationFrame(()=>requestAnimationFrame(resolve)))');
  const shot=await call('Page.captureScreenshot',{format:'png'});fs.writeFileSync(path.join(here,'browser_preview.png'),Buffer.from(shot.data,'base64'));
  const links=await evaluate('Array.from(document.querySelectorAll("a[href]")).map(a=>a.getAttribute("href"))');
  const missing=links.filter(link=>!fs.existsSync(path.resolve(here,link)));
  if(missing.length||errors.length)throw new Error(JSON.stringify({missing,errors}));
  fs.writeFileSync(path.join(here,'browser_checks.json'),JSON.stringify({status:'PASS',cases,missing_links:missing,javascript_errors:errors},null,2)+'\n');
  console.log(`PASS: all ${count} profile selectors, coordinated 3-D marker, mean line, sampled amplitude, and local links.`);
}finally{
  if(socket)socket.close();proc.kill();await new Promise(resolve=>proc.once('exit',resolve));fs.rmSync(temporary,{recursive:true,force:true});
}
