PAGE = r'''<!doctype html>
<html lang="fr"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>gueye — interface interactive</title>
<style>
:root{--bg:#f6f7f9;--panel:#fff;--ink:#1c2330;--mute:#667085;--line:#d9dde4;--acc:#2456d6;--accink:#fff;--warn:#b54708;--err:#b42318;--ok:#067647}
@media (prefers-color-scheme:dark){:root:not([data-theme=light]){--bg:#10141b;--panel:#171c26;--ink:#e6e9ef;--mute:#98a2b3;--line:#2b3342;--acc:#6d93ff;--accink:#0b1020;--warn:#f79009;--err:#f97066;--ok:#47cd89}}
:root[data-theme=dark]{--bg:#10141b;--panel:#171c26;--ink:#e6e9ef;--mute:#98a2b3;--line:#2b3342;--acc:#6d93ff;--accink:#0b1020;--warn:#f79009;--err:#f97066;--ok:#47cd89}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:14px/1.4 system-ui,Segoe UI,sans-serif}
header{display:flex;align-items:center;gap:12px;padding:8px 16px;border-bottom:1px solid var(--line);background:var(--panel)}
header h1{font-size:16px;margin:0}header .sp{flex:1}
main{display:grid;grid-template-columns:330px 1fr;gap:12px;padding:12px;max-width:1500px;margin:0 auto}
@media (max-width:900px){main{grid-template-columns:1fr}}
.card{background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:10px 12px;margin-bottom:10px}
.card h2{font-size:12px;text-transform:uppercase;letter-spacing:.06em;color:var(--mute);margin:0 0 8px}
label{display:block;color:var(--mute);font-size:12px;margin:6px 0 2px}
input,select,button{font:inherit;color:var(--ink);background:var(--bg);border:1px solid var(--line);border-radius:6px;padding:4px 6px}
input[type=number]{width:72px}input[type=range]{width:100%;padding:0}
button{cursor:pointer}button.pri{background:var(--acc);color:var(--accink);border-color:var(--acc)}
.row{display:flex;gap:6px;align-items:center;flex-wrap:wrap}.row>*{flex:0 0 auto}
.bc{display:grid;grid-template-columns:44px 1fr 1fr;gap:6px;align-items:start;margin-bottom:6px}
.bc .ax{color:var(--mute);padding-top:5px}.bcp{display:flex;flex-direction:column;gap:3px}
.tabs{display:flex;gap:4px;margin-bottom:8px}.tabs button{border-radius:6px 6px 0 0}.tabs button.on{background:var(--acc);color:var(--accink);border-color:var(--acc)}
#view{position:relative}canvas{width:100%;display:block;background:var(--panel);border:1px solid var(--line);border-radius:8px}
#stat{min-height:20px;font-size:12px;color:var(--mute);margin:6px 0}.warn{color:var(--warn)}.err{color:var(--err)}.ok{color:var(--ok)}
#msg{margin:6px 0;padding:8px 10px;border-radius:8px;display:none;border:1px solid var(--err);color:var(--err)}
.hint{font-size:12px;color:var(--mute)}code{background:var(--bg);padding:0 4px;border-radius:4px}
.res span{display:inline-block;margin:2px 4px 2px 0;padding:1px 6px;border:1px solid var(--line);border-radius:10px;cursor:pointer;font-size:12px}
</style></head><body>
<header><h1>gueye</h1><span class="hint">inversion semi-analytique de Laplace / Helmholtz · 1D–3D</span><span class="sp"></span>
<select id="preset" title="Exemples"></select><button id="theme" title="Clair / sombre">◐</button></header>
<main>
<aside>
<div class="card"><h2>Domaine</h2>
 <div class="row"><label>Dimension</label><select id="dim"><option>1</option><option selected>2</option><option>3</option></select></div>
 <div id="geo"></div>
</div>
<div class="card"><h2>Équation</h2>
 <label>k (nombre d'onde) : <b id="kv"></b></label><input type="range" id="k" min="0" max="30" step="0.01" value="9">
 <div class="row"><label style="margin:0">k =</label><input type="number" id="kn" step="0.01" value="9"><label style="margin:0">μ =</label><input type="number" id="mu" step="0.1" value="0" min="0"></div>
 <div class="hint">−Δu + μ²u − k²u = source ponctuelle. Garde-fou actif : une résonance exacte est signalée.</div>
</div>
<div class="card"><h2>Conditions aux limites</h2><div id="bcs"></div>
 <div class="hint">D Dirichlet · N Neumann · R Robin (α complexe) · A transparente exacte (dernier axe) · P PML (m mailles, réflexion R₀)</div></div>
<div class="card"><h2>Source et obstacle</h2>
 <div class="row" id="srcrow"></div>
 <div class="hint">Cliquez sur le champ pour placer la source ; Maj+clic pour le récepteur (balayage).</div>
 <div id="obsbox"><label><input type="checkbox" id="obs"> Obstacle de Dirichlet circulaire (2D)</label>
 <div class="row"><label style="margin:0">x</label><input type="number" id="ox" step="0.05"><label style="margin:0">y</label><input type="number" id="oy" step="0.05"><label style="margin:0">r</label><input type="number" id="or" step="0.05" min="0.02"></div></div>
</div>
<div class="card"><h2>Affichage</h2>
 <div class="row"><select id="qty"><option value="re">Partie réelle</option><option value="im">Partie imaginaire</option><option value="abs">Module</option><option value="log">Module (log)</option></select>
 <label style="margin:0"><input type="checkbox" id="auto" checked> calcul auto</label><button id="go" class="pri">Calculer</button></div>
 <div id="slicebox" style="display:none"><label>Coupe 3D : plan <select id="sax"><option value="z">z</option><option value="y">y</option><option value="x">x</option></select> <b id="sv"></b></label><input type="range" id="slice" min="0" max="10" value="5"></div>
</div>
</aside>
<section>
<div class="tabs"><button data-t="field" class="on">Champ</button><button data-t="sweep">Balayage en fréquence</button><button data-t="wave">Onde (temps)</button></div>
<div id="msg"></div>
<div id="view"><canvas id="cv" width="900" height="620"></canvas></div>
<div id="stat"></div>
<div id="sweepctl" style="display:none" class="card"><div class="row"><label style="margin:0">k min</label><input type="number" id="kmin" value="1" step="0.5"><label style="margin:0">k max</label><input type="number" id="kmax" value="15" step="0.5"><label style="margin:0">points</label><input type="number" id="kn2" value="200" step="10"><button id="dosweep" class="pri">Balayer</button></div><div class="hint">Récepteur : Maj+clic sur l'onglet Champ. Les traits verticaux sont les résonances exactes du schéma (si l'opérateur est auto-adjoint, sans PML/transparente).</div><div class="res" id="resl"></div></div>
<div id="wavectl" style="display:none" class="card"><div class="row"><label style="margin:0">c</label><input type="number" id="wc" value="1" step="0.1"><label style="margin:0">f₀</label><input type="number" id="wf" value="2" step="0.5"><label style="margin:0">T</label><input type="number" id="wT" value="3" step="0.5"><label style="margin:0">amort. η</label><input type="number" id="we" value="0" step="0.05"><button id="dowave" class="pri">Simuler</button><button id="play">▶ / ❚❚</button></div><input type="range" id="wt" min="0" max="0" value="0"><div class="hint">Impulsion de Ricker en la source (Newmark). Les PML ne sont pas disponibles pas à pas : utilisez D, N, R ou A.</div></div>
</section>
</main>
<script>
const $=id=>document.getElementById(id);
const FACES={1:[["gauche","droite"]],2:[["gauche","droite"],["bas","haut"]],3:[["gauche","droite"],["avant","arrière"],["bas","haut"]]};
const AX=["x","y","z"];
let S={dim:2,size:[3,2,1],cells:[150,100,30],bc:null,src:[0.9,0.7,0.5],rec:[2.0,1.2,0.5],data:null,tab:"field",wave:null,sweep:null,req:0,timer:null,playing:null};
const PRESETS={
"Cavité Dirichlet":{dim:2,size:[3,2],cells:[150,100],k:9,bc:[[{type:"D"},{type:"D"}],[{type:"D"},{type:"D"}]],src:[0.9,0.7]},
"Neumann + Robin":{dim:2,size:[3,2],cells:[150,100],k:6,bc:[[{type:"N"},{type:"R",re:3,im:0}],[{type:"N"},{type:"D"}]],src:[0.9,0.7]},
"Guide ouvert en haut (A)":{dim:2,size:[4,2],cells:[200,100],k:12,bc:[[{type:"D"},{type:"D"}],[{type:"D"},{type:"A"}]],src:[2,0.4]},
"Domaine ouvert (PML partout)":{dim:2,size:[4,3],cells:[160,120],k:16,bc:[[{type:"P",m:16,R0:1e-6},{type:"P",m:16,R0:1e-6}],[{type:"P",m:16,R0:1e-6},{type:"P",m:16,R0:1e-6}]],src:[1.2,1.0]},
"Obstacle en domaine ouvert":{dim:2,size:[4,3],cells:[160,120],k:16,bc:[[{type:"P",m:16,R0:1e-6},{type:"P",m:16,R0:1e-6}],[{type:"P",m:16,R0:1e-6},{type:"P",m:16,R0:1e-6}]],src:[1.2,1.0],obs:{cx:2.4,cy:1.6,r:0.35}},
"1D Robin":{dim:1,size:[1],cells:[400],k:20,bc:[[{type:"D"},{type:"R",re:2,im:0}]],src:[0.3]},
"3D cavité":{dim:3,size:[1,1,1],cells:[32,32,32],k:6,bc:[[{type:"D"},{type:"D"}],[{type:"D"},{type:"D"}],[{type:"D"},{type:"D"}]],src:[0.3,0.4,0.5]}
};
function el(t,a,h){const e=document.createElement(t);for(const k in a||{})e[k]=a[k];if(h)e.innerHTML=h;return e}
function buildGeo(){const g=$("geo");g.innerHTML="";for(let i=0;i<S.dim;i++){const r=el("div",{className:"row"});r.append(el("label",{style:"margin:0;width:14px"},AX[i]),
 Object.assign(el("input",{type:"number",value:S.size[i],step:0.1,min:0.1,id:"L"+i,title:"longueur"}),{}),el("span",{className:"hint"},"cells"),el("input",{type:"number",value:S.cells[i],step:10,min:4,id:"M"+i}));g.append(r)}
 for(let i=0;i<S.dim;i++){$("L"+i).onchange=()=>{S.size[i]=+$("L"+i).value;changed()};$("M"+i).onchange=()=>{S.cells[i]=+$("M"+i).value;changed()}}}
function bcHTML(i,side,b){const last=(i==S.dim-1);const id=`b${i}${side}`;
 const opts=["D","N","R"].concat(last?["A"]:[]).concat(["P"]).map(t=>`<option ${b.type==t?"selected":""}>${t}</option>`).join("");
 return `<div class="bcp"><select id="${id}t">${opts}</select><span id="${id}x"></span></div>`}
function buildBC(){const c=$("bcs");c.innerHTML="";for(let i=0;i<S.dim;i++){const r=el("div",{className:"bc"});r.append(el("div",{className:"ax"},AX[i]));
 for(let s=0;s<2;s++){const d=el("div");d.innerHTML=`<div class="hint">${FACES[S.dim][i][s]}</div>`+bcHTML(i,s,S.bc[i][s]);r.append(d)}c.append(r);
 for(let s=0;s<2;s++){const id=`b${i}${s}`;$(id+"t").onchange=()=>{S.bc[i][s]={type:$(id+"t").value,m:16,R0:1e-6,re:1,im:0};paramUI(i,s);changed()};paramUI(i,s)}}}
function paramUI(i,s){const b=S.bc[i][s],id=`b${i}${s}`,x=$(id+"x");x.innerHTML="";
 const num=(k,l,st)=>{const w=el("label",{style:"margin:0;display:inline-block"},l+" ");const n=el("input",{type:"number",value:b[k],step:"any",style:"width:"+(k=="R0"?"96px":"62px")});n.onchange=()=>{b[k]=+n.value;changed()};w.append(n);return w};
 if(b.type=="R"){b.re=b.re??1;b.im=b.im??0;x.append(num("re","α′",0.5),num("im","α″",0.5))}
 if(b.type=="P"){b.m=b.m??16;b.R0=b.R0??1e-6;x.append(num("m","m",2),num("R0","R₀",1e-6))}}
function buildSrc(){const r=$("srcrow");r.innerHTML="";for(let i=0;i<S.dim;i++){const n=el("input",{type:"number",step:0.05,value:S.src[i].toFixed(3)});n.onchange=()=>{S.src[i]=+n.value;changed()};r.append(el("label",{style:"margin:0"},AX[i]),n)}
 $("obsbox").style.display=S.dim==2?"block":"none";$("slicebox").style.display=S.dim==3?"block":"none"}
function load(p){S.dim=p.dim;S.size=p.size.concat([1,1,1]).slice(0,3);S.cells=p.cells.concat([30,30,30]).slice(0,3);S.bc=JSON.parse(JSON.stringify(p.bc));S.src=p.src.concat([0.5,0.5,0.5]).slice(0,3);
 S.rec=S.src.map(v=>v);$("dim").value=p.dim;$("k").value=$("kn").value=p.k;$("mu").value=0;$("obs").checked=!!p.obs;if(p.obs){$("ox").value=p.obs.cx;$("oy").value=p.obs.cy;$("or").value=p.obs.r}else{$("ox").value=p.size[0]*0.6;$("oy").value=(p.size[1]||1)*0.55;$("or").value=0.3}
 buildGeo();buildBC();buildSrc();updK();solve()}
function defaultBC(d){const D=()=>({type:"D"});const b=[];for(let i=0;i<d;i++)b.push([D(),D()]);return b}
function req(){const d=S.dim;return {size:S.size.slice(0,d),cells:S.cells.slice(0,d),bc:S.bc.slice(0,d),k:+$("kn").value,mu:+$("mu").value,src:S.src.slice(0,d),rec:S.rec.slice(0,d),
 obstacle:($("obs").checked&&d==2)?{cx:+$("ox").value,cy:+$("oy").value,r:+$("or").value}:null}}
async function post(u,b){const r=await fetch(u,{method:"POST",body:JSON.stringify(b)});return r.json()}
function msg(t,cls){const m=$("msg");if(!t){m.style.display="none";return}m.style.display="block";m.className="";m.style.color=cls=="warn"?"var(--warn)":"var(--err)";m.style.borderColor=m.style.color;m.textContent=t}
function b64(s){const bin=atob(s),n=bin.length,u=new Uint8Array(n);for(let i=0;i<n;i++)u[i]=bin.charCodeAt(i);return new Float32Array(u.buffer)}
let busy=false,again=false;
async function solve(){if(S.tab=="sweep"||S.tab=="wave"){return}
 if(busy){again=true;return}busy=true;const id=++S.req;
 try{$("stat").textContent="calcul…";const r=await post("/api/solve",req());if(id!=S.req&&again){}
  if(r.error){msg(r.error);$("stat").textContent="";S.data=null;draw();return}
  msg(r.warnings&&r.warnings.length?"Avertissement : "+r.warnings[0]:"", "warn");
  S.data={dim:r.dim,shape:r.shape,re:b64(r.re),im:b64(r.im),axes:r.axes,ms:r.ms};
  const ns=S.dim==3?S.data.shape[$("sax").selectedIndex==0?0:($("sax").value=="y"?1:2)]:1;setSliceMax();
  draw();$("stat").textContent=`${r.shape.join(" × ")} nœuds · ${r.ms} ms`}
 catch(e){msg("Erreur de communication avec le serveur : "+e)}finally{busy=false;if(again){again=false;solve()}}}
function changed(){buildSrc();if($("auto").checked)debounce(solve,150)}
function debounce(f,t){clearTimeout(S.timer);S.timer=setTimeout(f,t)}
// --- couleurs
const RDBU=[[5,48,97],[33,102,172],[103,169,207],[209,229,240],[247,247,247],[253,219,199],[239,138,98],[178,24,43],[103,0,31]];
const VIR=[[68,1,84],[59,82,139],[33,145,140],[94,201,98],[253,231,37]];
function cmap(L,t){t=Math.max(0,Math.min(1,t));const x=t*(L.length-1),i=Math.min(L.length-2,Math.floor(x)),f=x-i;return L[i].map((c,j)=>c+(L[i+1][j]-c)*f)}
function values(){const D=S.data;if(!D)return null;const q=$("qty").value;const n=D.re.length;const v=new Float32Array(n);
 for(let i=0;i<n;i++){const a=D.re[i],b=D.im[i];v[i]=q=="re"?a:q=="im"?b:q=="abs"?Math.hypot(a,b):Math.log10(Math.hypot(a,b)+1e-300)}return v}
function pct(v,p){const a=Array.from(v).filter(Number.isFinite).sort((x,y)=>x-y);return a[Math.min(a.length-1,Math.floor(p*(a.length-1)))]}
function sliceOf(v){const D=S.data;if(D.dim<3)return{w:D.shape[D.dim-1],h:D.dim==2?D.shape[0]:1,get:(i,j)=>v[j*D.shape[D.dim-1]+i]};
 const [nz,ny,nx]=D.shape,ax=$("sax").value,k=+$("slice").value;
 if(ax=="z")return{w:nx,h:ny,get:(i,j)=>v[(k*ny+j)*nx+i]};if(ax=="y")return{w:nx,h:nz,get:(i,j)=>v[(j*ny+k)*nx+i]};return{w:ny,h:nz,get:(i,j)=>v[(j*ny+i)*nx+k]}}
function setSliceMax(){if(S.dim!=3||!S.data)return;const [nz,ny,nx]=S.data.shape,ax=$("sax").value,n=ax=="z"?nz:ax=="y"?ny:nx;const s=$("slice");s.max=n-1;if(+s.value>n-1)s.value=n>>1;
 const L={x:S.size[0],y:S.size[1],z:S.size[2]}[ax];$("sv").textContent=`${ax} = ${(+s.value/(n-1)*L).toFixed(2)}`}
const M={l:56,r:70,t:14,b:36};
function draw(){const cv=$("cv"),g=cv.getContext("2d");const css=getComputedStyle(document.documentElement);const ink=css.getPropertyValue("--ink"),mute=css.getPropertyValue("--mute"),panel=css.getPropertyValue("--panel");
 const W=cv.width=cv.clientWidth*devicePixelRatio,H=cv.height=Math.round(cv.clientWidth*(S.tab=="field"?0.62:0.5))*devicePixelRatio;cv.style.height=H/devicePixelRatio+"px";g.setTransform(devicePixelRatio,0,0,devicePixelRatio,0,0);
 const w=W/devicePixelRatio,h=H/devicePixelRatio;g.fillStyle=panel;g.fillRect(0,0,w,h);g.font="12px system-ui";
 if(S.tab=="sweep")return drawSweep(g,w,h,ink,mute);if(S.tab=="wave")return drawWave(g,w,h,ink,mute);
 if(!S.data){g.fillStyle=mute;g.fillText("Aucun champ (voir le message).",M.l,40);return}
 const v=values();const q=$("qty").value;
 if(S.dim==1)return draw1D(g,w,h,v,ink,mute);
 const sl=sliceOf(v);const pw=w-M.l-M.r,ph=h-M.t-M.b;
 const vals=[];for(let j=0;j<sl.h;j++)for(let i=0;i<sl.w;i++)vals.push(sl.get(i,j));
 let lo,hi,L=RDBU;if(q=="re"||q=="im"){const m=pct(vals.map(Math.abs),0.99)||1;lo=-m;hi=m}else{lo=pct(vals,0.01);hi=pct(vals,0.995);L=VIR}if(hi==lo)hi=lo+1;
 const off=document.createElement("canvas");off.width=sl.w;off.height=sl.h;const og=off.getContext("2d"),im=og.createImageData(sl.w,sl.h);
 for(let j=0;j<sl.h;j++)for(let i=0;i<sl.w;i++){const c=cmap(L,(sl.get(i,j)-lo)/(hi-lo)),p=((sl.h-1-j)*sl.w+i)*4;im.data[p]=c[0];im.data[p+1]=c[1];im.data[p+2]=c[2];im.data[p+3]=255}
 og.putImageData(im,0,0);g.imageSmoothingEnabled=true;g.drawImage(off,M.l,M.t,pw,ph);
 S.geom={x0:M.l,y0:M.t,pw,ph,sl,lo,hi};ticks(g,pw,ph,ink,mute);colorbar(g,w,ph,lo,hi,L,ink);marks(g,pw,ph)
 const ax=axesOfSlice();}
function axesOfSlice(){if(S.dim==2)return[0,1];const a=$("sax").value;return a=="z"?[0,1]:a=="y"?[0,2]:[1,2]}
function ticks(g,pw,ph,ink,mute){const[a,b]=axesOfSlice();g.fillStyle=mute;g.strokeStyle=mute;g.lineWidth=1;g.strokeRect(M.l,M.t,pw,ph);
 for(let i=0;i<=5;i++){const x=M.l+pw*i/5;g.fillText((S.size[a]*i/5).toFixed(2),x-12,M.t+ph+16);const y=M.t+ph-ph*i/5;g.fillText((S.size[b]*i/5).toFixed(2),M.l-36,y+4)}
 g.fillStyle=ink;g.fillText(AX[a],M.l+pw/2,M.t+ph+31);g.fillText(AX[b],6,M.t+ph/2)}
function colorbar(g,w,ph,lo,hi,L,ink){const x=w-M.r+14,wd=14;for(let j=0;j<ph;j++){const c=cmap(L,1-j/ph);g.fillStyle=`rgb(${c[0]|0},${c[1]|0},${c[2]|0})`;g.fillRect(x,M.t+j,wd,1)}
 g.fillStyle=ink;g.fillText(fmt(hi),x+wd+4,M.t+10);g.fillText(fmt(lo),x+wd+4,M.t+ph);g.fillText(fmt((lo+hi)/2),x+wd+4,M.t+ph/2+4)}
function fmt(v){return Math.abs(v)>=100||(Math.abs(v)<0.01&&v!=0)?v.toExponential(1):v.toPrecision(3)}
function toPix(p){const[a,b]=axesOfSlice(),G=S.geom;return[G.x0+G.pw*p[a]/S.size[a],G.y0+G.ph-G.ph*p[b]/S.size[b]]}
function marks(g,pw,ph){const o=toPix(S.src);g.strokeStyle="#111";g.fillStyle="#ffd43b";g.lineWidth=1.5;g.beginPath();g.arc(o[0],o[1],6,0,7);g.fill();g.stroke();
 const r=toPix(S.rec);g.fillStyle="#74c0fc";g.beginPath();g.rect(r[0]-5,r[1]-5,10,10);g.fill();g.stroke();
 if($("obs").checked&&S.dim==2){const c=toPix([+$("ox").value,+$("oy").value]);g.strokeStyle="#111";g.setLineDash([4,3]);g.beginPath();g.ellipse(c[0],c[1],+$("or").value/S.size[0]*pw,+$("or").value/S.size[1]*ph,0,0,7);g.stroke();g.setLineDash([])}}
function draw1D(g,w,h,v,ink,mute){const n=v.length,pw=w-M.l-20,ph=h-M.t-M.b;let lo=Math.min(...v),hi=Math.max(...v);if(lo==hi)hi=lo+1;const p=(hi-lo)*0.08;lo-=p;hi+=p;
 g.strokeStyle=mute;g.strokeRect(M.l,M.t,pw,ph);g.strokeStyle="#2456d6";g.lineWidth=1.6;g.beginPath();for(let i=0;i<n;i++){const x=M.l+pw*i/(n-1),y=M.t+ph-ph*(v[i]-lo)/(hi-lo);i?g.lineTo(x,y):g.moveTo(x,y)}g.stroke();
 g.fillStyle=mute;for(let i=0;i<=5;i++){g.fillText((S.size[0]*i/5).toFixed(2),M.l+pw*i/5-12,M.t+ph+16)}g.fillText(fmt(hi),4,M.t+10);g.fillText(fmt(lo),4,M.t+ph);
 const sx=M.l+pw*S.src[0]/S.size[0];g.strokeStyle="#e8a600";g.beginPath();g.moveTo(sx,M.t);g.lineTo(sx,M.t+ph);g.stroke();S.geom={x0:M.l,y0:M.t,pw,ph,n}}
function coordOf(ev){const cv=$("cv"),r=cv.getBoundingClientRect(),x=ev.clientX-r.left,y=ev.clientY-r.top,G=S.geom;if(!G||S.tab!="field")return null;
 const fx=(x-G.x0)/G.pw,fy=1-(y-G.y0)/G.ph;if(fx<0||fx>1||(S.dim>1&&(fy<0||fy>1)))return null;return{fx,fy}}
function pointFrom(c){const p=S.src.slice();if(S.dim==1)p[0]=c.fx*S.size[0];else{const[a,b]=axesOfSlice();p[a]=c.fx*S.size[a];p[b]=c.fy*S.size[b]}return p}
$("cv").onclick=ev=>{const c=coordOf(ev);if(!c)return;const p=pointFrom(c);if(ev.shiftKey){S.rec=p}else{S.src=p}buildSrc();ev.shiftKey?draw():(msg(""),solve())};
$("cv").onmousemove=ev=>{const c=coordOf(ev);if(!c||!S.data)return;const v=values(),G=S.geom;let val,pos;
 if(S.dim==1){const i=Math.round(c.fx*(v.length-1));val=`${$("qty").value} = ${fmt(v[i])}`;pos=`x = ${(c.fx*S.size[0]).toFixed(3)}`}
 else{const sl=G.sl,i=Math.round(c.fx*(sl.w-1)),j=Math.round(c.fy*(sl.h-1));const p=pointFrom(c);pos=p.slice(0,S.dim).map((u,k)=>AX[k]+" = "+u.toFixed(3)).join(", ");const re=S.data.re,im=S.data.im;val=`${$("qty").value} = ${fmt(sl.get(i,j))}`}
 $("stat").textContent=`${pos}   ${val}   · ${S.data.ms} ms`};
// --- balayage
async function doSweep(){msg("");$("stat").textContent="balayage…";const b=req();b.kmin=+$("kmin").value;b.kmax=+$("kmax").value;b.n=+$("kn2").value;const r=await post("/api/sweep",b);
 if(r.error){msg(r.error);return}S.sweep=r;const rl=$("resl");rl.innerHTML="";if(r.resonances)r.resonances.forEach(k=>{const s=el("span",{},"k = "+k.toFixed(3));s.onclick=()=>{$("kn").value=$("k").value=k.toFixed(4);tab("field");solve()};rl.append(s)});
 $("stat").textContent=`${r.ks.length} fréquences · ${r.ms} ms`;draw()}
function drawSweep(g,w,h,ink,mute){const r=S.sweep;const pw=w-M.l-20,ph=h-M.t-M.b;g.strokeStyle=mute;g.strokeRect(M.l,M.t,pw,ph);
 if(!r){g.fillStyle=mute;g.fillText("Réglez k min / k max puis « Balayer » : |G(récepteur ; source)| en fonction de k.",M.l+10,40);return}
 const ys=r.mag.filter(v=>v!=null).map(v=>Math.log10(v+1e-300));let lo=Math.min(...ys),hi=Math.max(...ys);if(hi==lo)hi=lo+1;
 const X=k=>M.l+pw*(k-r.ks[0])/(r.ks[r.ks.length-1]-r.ks[0]),Y=v=>M.t+ph-ph*(Math.log10(v+1e-300)-lo)/(hi-lo);
 if(r.resonances){g.strokeStyle="#e8a60088";r.resonances.forEach(k=>{g.beginPath();g.moveTo(X(k),M.t);g.lineTo(X(k),M.t+ph);g.stroke()})}
 g.strokeStyle="#2456d6";g.lineWidth=1.6;g.beginPath();let on=false;r.mag.forEach((v,i)=>{if(v==null){on=false;return}const x=X(r.ks[i]),y=Y(v);on?g.lineTo(x,y):g.moveTo(x,y);on=true});g.stroke();
 g.fillStyle=mute;for(let i=0;i<=6;i++){const k=r.ks[0]+(r.ks[r.ks.length-1]-r.ks[0])*i/6;g.fillText(k.toFixed(1),X(k)-8,M.t+ph+16)}g.fillText("k",M.l+pw/2,M.t+ph+31);g.fillText("log₁₀|G|",4,M.t+12);g.fillText(hi.toFixed(1),4,M.t+26);g.fillText(lo.toFixed(1),4,M.t+ph)}
// --- onde
async function doWave(){msg("");$("stat").textContent="simulation…";const b=req();b.c=+$("wc").value;b.f0=+$("wf").value;b.T=+$("wT").value;b.eta=+$("we").value;const r=await post("/api/wave",b);
 if(r.error){msg(r.error);$("stat").textContent="";return}S.wave={shape:r.shape,t:r.t,f:b64(r.frames),dim:r.dim};const s=$("wt");s.max=r.t.length-1;s.value=0;$("stat").textContent=`${r.t.length} images · ${r.ms} ms`;playWave(true)}
function drawWave(g,w,h,ink,mute){const W=S.wave;if(!W){g.fillStyle=mute;g.fillText("Réglez l'impulsion puis « Simuler ».",M.l,40);return}const sh=W.shape,fi=+$("wt").value;let nx,ny,off;
 if(W.dim==2){ny=sh[1];nx=sh[2];off=fi*ny*nx}else{const nz=sh[1];ny=sh[2];nx=sh[3];off=fi*nz*ny*nx+(nz>>1)*ny*nx}
 let m=0;const nfr=sh[0],per=(W.dim==2?ny*nx:sh[1]*ny*nx);for(let f=0;f<nfr;f+=Math.max(1,nfr>>3))for(let i=0;i<per;i+=7)m=Math.max(m,Math.abs(W.f[f*per+i]));m=m||1;
 const pw=w-M.l-M.r,ph=h-M.t-M.b,off2=document.createElement("canvas");off2.width=nx;off2.height=ny;const og=off2.getContext("2d"),im=og.createImageData(nx,ny);
 for(let j=0;j<ny;j++)for(let i=0;i<nx;i++){const c=cmap(RDBU,0.5+0.5*W.f[off+j*nx+i]/(0.6*m)),p=((ny-1-j)*nx+i)*4;im.data[p]=c[0];im.data[p+1]=c[1];im.data[p+2]=c[2];im.data[p+3]=255}
 og.putImageData(im,0,0);g.drawImage(off2,M.l,M.t,pw,ph);g.strokeStyle=mute;g.strokeRect(M.l,M.t,pw,ph);g.fillStyle=ink;g.fillText(`t = ${W.t[fi].toFixed(3)}`+(W.dim==3?"  (plan z médian)":""),M.l+6,M.t+16)}
function playWave(on){clearInterval(S.playing);S.playing=null;if(!on||!S.wave)return;S.playing=setInterval(()=>{const s=$("wt");s.value=(+s.value+1)%(+s.max+1);draw()},60)}
// --- onglets / événements
function tab(t){S.tab=t;document.querySelectorAll(".tabs button").forEach(b=>b.classList.toggle("on",b.dataset.t==t));$("sweepctl").style.display=t=="sweep"?"block":"none";$("wavectl").style.display=t=="wave"?"block":"none";if(t!="wave")playWave(false);msg("");if(t=="field")solve();else draw()}
document.querySelectorAll(".tabs button").forEach(b=>b.onclick=()=>tab(b.dataset.t));
$("dosweep").onclick=doSweep;$("dowave").onclick=doWave;$("play").onclick=()=>playWave(!S.playing);$("wt").oninput=draw;
function updK(){$("kv").textContent=(+$("kn").value).toFixed(3)}
$("k").oninput=()=>{$("kn").value=$("k").value;updK();changed()};$("kn").onchange=()=>{$("k").value=$("kn").value;updK();changed()};
$("mu").onchange=changed;$("go").onclick=solve;$("qty").onchange=draw;$("slice").oninput=()=>{setSliceMax();draw()};$("sax").onchange=()=>{setSliceMax();draw()};
["obs","ox","oy","or"].forEach(i=>$(i).onchange=changed);
$("dim").onchange=()=>{const d=+$("dim").value;S.dim=d;S.bc=defaultBC(d);S.src=S.size.slice(0,d).map(v=>v*0.4).concat([0.5,0.5,0.5]).slice(0,3);if(d==3)S.cells=[32,32,32];if(d==1)S.cells[0]=400;buildGeo();buildBC();buildSrc();changed()};
$("theme").onclick=()=>{const r=document.documentElement;const dark=getComputedStyle(r).getPropertyValue("--bg").trim()=="#10141b";r.dataset.theme=dark?"light":"dark";draw()};
const ps=$("preset");Object.keys(PRESETS).forEach(n=>ps.append(el("option",{},n)));ps.onchange=()=>load(PRESETS[ps.value]);
addEventListener("resize",draw);
load(PRESETS["Cavité Dirichlet"]);
</script></body></html>'''
