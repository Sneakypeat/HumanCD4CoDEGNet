"""Generate a self-advancing teleprompter/slide deck (self-contained HTML, figures base64-inlined)."""
import base64, os, json
FIG=os.path.expanduser("~/CoDEG_Tcell/artifacts/figures")
def uri(p): return "data:image/png;base64,"+base64.b64encode(open(f"{FIG}/{p}","rb").read()).decode()
f1,f2,f3,f5=uri("fig_arch1_topology.png"),uri("fig_arch2_rewiring.png"),uri("fig_arch3_generality.png"),uri("fig_arch5_disease.png")

SLIDES=[
 {"idx":"","label":"Claude Hackathon · Life Sciences","type":"title","eyebrow":"Researcher track",
  "headline":"HumanCD4CoDEGNet","sub":"Syed Sabih ur Rehman · PhD student, United Arab Emirates University",
  "script":"Hi, I'm <b>Syed Sabih ur Rehman, a PhD student at United Arab Emirates University</b>. My project is HumanCD4CoDEGNet."},
 {"idx":"","label":"The hook","type":"text",
  "headline":"The genes that switch on when an immune cell activates are druggable disease genes.",
  "script":"When a human T cell switches on to fight an infection, a specific set of control genes comes online — and they turn out to be disproportionately <b>druggable</b>, and disproportionately genes that <b>cause disease</b> when they break. Genes like ZAP70 and ITK. I found them by mapping the entire causal control system of a human immune cell. Here's how."},
 {"idx":"","label":"The data","type":"text",
  "headline":"A causal map: 22M CD4⁺ T cells, every gene knocked down, at rest and after activation.",
  "script":"You need causality, not correlation — you perturb genes and watch what moves. The Marson lab released exactly that: a genome-scale Perturb-seq atlas — twenty-two million primary human CD4 T cells, every gene knocked down, at rest and after activation. For the first time we can draw the causal wiring diagram of a human immune cell, and redraw it in each state."},
 {"idx":"01","label":"The architecture","type":"fig","img":f1,
  "headline":"Hub-dominated: the top 5% of regulators drive ~78% of all causal effects.",
  "script":"So I drew it. The network has a very particular shape: it's <b>hub-dominated</b>. The top five percent of regulators drive nearly eighty percent of all causal effects; most genes do almost nothing when knocked down. That's exactly the shape a Pritchard-lab theory paper predicted last year from heritability alone — confirmed causally here for the first time."},
 {"idx":"02","label":"The twist","type":"fig","img":f2,
  "headline":"Shape-invariant, identity-labile: the shape holds while up to 59% of hubs turn over.",
  "script":"When the cell activates, does the architecture change? It splits in two. The <b>shape</b> doesn't change — the concentration stays pinned at rest, eight hours, forty-eight hours. But the <b>identity</b> of the hubs turns over almost completely — up to fifty-nine percent replaced. The cell keeps the shape of its control while swapping out who's in control. The new hubs are the T-cell receptor signalosome — ZAP70, LAT, ITK — switching on."},
 {"idx":"03","label":"The payoff","type":"fig","img":f5,
  "headline":"State-specific hubs are ~2× enriched for druggable disease genes: ZAP70, ITK, LCK, PTPRC, IL12RB2.",
  "script":"Here's why that matters. Those activation-specific hubs are about <b>twice as enriched</b> for genes that cause disease when mutated as the stable hubs — and they're druggable: ZAP70, ITK, LCK, PTPRC, IL12RB2. The control layer activation builds is exactly where druggable disease risk concentrates — targets that only matter in the activated state."},
 {"idx":"04","label":"Generality & close","type":"fig","img":f3,
  "headline":"Not a T-cell quirk — and it survives every confound. Reproducible from public data.",
  "script":"It's not a T-cell quirk — the same architecture appears in a completely different cell type. And it survives every confound I tested — reproducible from public data, code on GitHub. It's a map of how a human immune cell controls itself, and a shortlist of where to intervene. Thanks for watching."},
]
DUR=[12,26,22,30,35,30,20]

TEMPLATE=r"""<title>HumanCD4CoDEGNet — presenter / teleprompter</title>
<style>
:root{--bg:#0B0E13;--fg:#EAECEF;--muted:#9AA4B0;--blue:#6FB3D6;--orange:#E8935A;--plate:#FBFBFC;--hair:rgba(255,255,255,.13);}
*{box-sizing:border-box;margin:0;padding:0}
html,body{height:100%}
body{background:var(--bg);color:var(--fg);font-family:system-ui,-apple-system,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;overflow:hidden}
.stage{height:100vh;display:flex;flex-direction:column;padding:clamp(14px,2.4vw,30px);gap:8px}
.top{display:flex;align-items:center;justify-content:space-between;gap:12px;font-family:ui-monospace,Menlo,Consolas,monospace;font-size:.8rem;color:var(--muted);flex:0 0 auto}
.top .label{color:var(--blue);letter-spacing:.08em;text-transform:uppercase}
.top .idx{color:var(--orange);font-weight:700}
.clock{font-variant-numeric:tabular-nums}
.main{flex:1 1 auto;display:flex;flex-direction:column;align-items:center;justify-content:center;gap:clamp(10px,2vw,20px);min-height:0;text-align:center}
.plate{background:var(--plate);border-radius:12px;padding:14px;max-width:min(1080px,92vw);max-height:54vh;display:flex}
.plate img{max-width:100%;max-height:100%;object-fit:contain;display:block}
.headline{font-size:clamp(1.25rem,2.7vw,2.1rem);font-weight:650;line-height:1.16;max-width:26ch;text-wrap:balance}
.title .big{font-size:clamp(2.3rem,6vw,4rem);font-weight:700;letter-spacing:-.02em}
.title .sub{color:var(--muted);font-size:clamp(1rem,2vw,1.35rem);margin-top:12px}
.title .eyebrow{font-family:ui-monospace,monospace;color:var(--blue);text-transform:uppercase;letter-spacing:.16em;font-size:.8rem;margin-bottom:16px}
.prompt{flex:0 0 auto;background:rgba(255,255,255,.05);border:1px solid var(--hair);border-radius:12px;padding:clamp(12px,1.8vw,20px);font-size:clamp(1.05rem,1.65vw,1.5rem);line-height:1.5;max-width:74ch;margin:0 auto;color:#EDF0F3}
.prompt.hidden{display:none}
.prompt b{color:var(--orange)}
.bar{flex:0 0 auto;height:4px;background:var(--hair);border-radius:2px;overflow:hidden}
.bar>i{display:block;height:100%;width:0;background:linear-gradient(90deg,var(--blue),var(--orange))}
.controls{flex:0 0 auto;display:flex;align-items:center;justify-content:center;gap:10px;flex-wrap:wrap}
button{background:rgba(255,255,255,.08);color:var(--fg);border:1px solid var(--hair);border-radius:8px;padding:8px 14px;font-size:.9rem;cursor:pointer}
button:hover{background:rgba(255,255,255,.16)}
.dots{flex:0 0 auto;display:flex;gap:7px;justify-content:center}
.dots i{width:9px;height:9px;border-radius:50%;background:var(--hair);cursor:pointer}
.dots i.on{background:var(--orange)}
.hint{flex:0 0 auto;color:var(--muted);font-size:.72rem;text-align:center;font-family:ui-monospace,monospace}
</style>
<div class="stage">
  <div class="top"><span><span class="idx" id="sidx"></span> <span class="label" id="slabel"></span></span>
    <span class="clock"><span id="clk">0:00</span> / 2:55</span></div>
  <div class="main" id="main"></div>
  <div class="prompt" id="prompt"></div>
  <div class="bar"><i id="barfill"></i></div>
  <div class="controls">
    <button id="prev">◀ Prev</button><button id="play">▶ Auto-play</button><button id="next">Next ▶</button>
    <button id="toggle">Hide script (T)</button><button id="reset">↺ Reset</button><button id="full">⛶ Full-screen</button>
  </div>
  <div class="dots" id="dots"></div>
  <div class="hint">→ / Space next · ← prev · P play-pause · T toggle script · R reset · F full-screen — auto-play paces you to 2:55</div>
</div>
<script>
const S=__SLIDES__, DUR=__DUR__;
let i=0,playing=false,timer=null,elapsed=0,slideStart=0;
const $=id=>document.getElementById(id);
const fmt=s=>{s=Math.max(0,Math.round(s));return Math.floor(s/60)+':'+String(s%60).padStart(2,'0')};
function render(){const s=S[i];$('sidx').textContent=s.idx||'';$('slabel').textContent=s.label;let h='';
 if(s.type==='title')h='<div class="title"><div class="eyebrow">'+(s.eyebrow||'')+'</div><div class="big">'+s.headline+'</div><div class="sub">'+(s.sub||'')+'</div></div>';
 else if(s.type==='fig')h='<div class="plate"><img src="'+s.img+'" alt=""></div><div class="headline">'+s.headline+'</div>';
 else h='<div class="headline" style="max-width:32ch">'+s.headline+'</div>';
 $('main').innerHTML=h;$('prompt').innerHTML=s.script;
 [...document.querySelectorAll('.dots i')].forEach((d,k)=>d.classList.toggle('on',k===i));}
function go(n){i=Math.max(0,Math.min(S.length-1,n));slideStart=elapsed;render();}
function next(){i<S.length-1?go(i+1):stop();}
function prev(){go(i-1);}
function tick(){elapsed+=0.2;$('clk').textContent=fmt(elapsed);const d=DUR[i]||20;
 $('barfill').style.width=Math.min(100,(elapsed-slideStart)/d*100)+'%';
 if(playing&&elapsed-slideStart>=d){i<S.length-1?go(i+1):stop();}}
function play(){playing=true;$('play').textContent='⏸ Pause';if(!timer)timer=setInterval(tick,200);}
function stop(){playing=false;$('play').textContent='▶ Auto-play';if(timer){clearInterval(timer);timer=null;}}
function reset(){stop();elapsed=0;slideStart=0;$('clk').textContent='0:00';go(0);}
$('prev').onclick=prev;$('next').onclick=next;$('reset').onclick=reset;
$('play').onclick=()=>playing?stop():play();
$('toggle').onclick=()=>{const p=$('prompt');p.classList.toggle('hidden');$('toggle').textContent=p.classList.contains('hidden')?'Show script (T)':'Hide script (T)';};
$('full').onclick=()=>{(document.fullscreenElement?document.exitFullscreen():document.documentElement.requestFullscreen()).catch(()=>{});};
S.forEach((s,k)=>{const d=document.createElement('i');d.onclick=()=>go(k);$('dots').appendChild(d);});
addEventListener('keydown',e=>{const k=e.key.toLowerCase();
 if(e.key==='ArrowRight'||e.key===' '){e.preventDefault();next();}else if(e.key==='ArrowLeft')prev();
 else if(k==='p'){playing?stop():play();}else if(k==='t')$('toggle').click();else if(k==='r')reset();else if(k==='f')$('full').click();});
render();
</script>"""
html=TEMPLATE.replace("__SLIDES__",json.dumps(SLIDES)).replace("__DUR__",json.dumps(DUR))
for out in [os.path.expanduser("~/CoDEG_Tcell/video/presenter.html"),
            "/private/tmp/claude-501/-Users-syedsabihurrehman-Desktop-Metagenome-Project/8c660892-f51b-4520-a39b-07443b4e6eb2/scratchpad/presenter.html"]:
    os.makedirs(os.path.dirname(out),exist_ok=True); open(out,"w").write(html)
print("wrote presenter.html | bytes:",len(html))
