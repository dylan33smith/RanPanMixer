import markdown, pathlib, re, html
D = pathlib.Path(__file__).resolve().parent / "src"
OUT = pathlib.Path(__file__).resolve().parent / "primer.html"

md = markdown.Markdown(extensions=["tables", "fenced_code", "toc", "attr_list", "sane_lists"],
                       extension_configs={"toc": {"toc_depth": "2-3", "anchorlink": False}})

order = ["01-landscape.md", "02-problem.md", "03-panmixer.md", "04-ours.md"]
bodies, tocs = [], []
for i, f in enumerate(order, 1):
    src = (D / f).read_text()
    md.reset()
    bodies.append(f'<section class="sec" id="sec{i}">' + md.convert(src) + "</section>")
    tocs.append(md.toc_tokens)

# Build our own nav from the h2/h3 structure
nav = ['<nav class="toc" aria-label="Contents"><p class="toc-h">Contents</p><ol class="toc-l1">']
for i, tt in enumerate(tocs, 1):
    for t in tt:
        nav.append(f'<li><a href="#{t["id"]}">{html.escape(re.sub("^Section [0-9]+ . ", "", t["name"]))}</a>')
        kids = [c for c in t["children"]]
        if kids:
            nav.append('<ol class="toc-l2">')
            for c in kids:
                nav.append(f'<li><a href="#{c["id"]}">{html.escape(c["name"])}</a></li>')
            nav.append("</ol>")
        nav.append("</li>")
nav.append("</ol></nav>")
nav = "".join(nav)

body = "\n".join(bodies)
# Callout detection: paragraphs/blockquotes that open with a role word get a semantic class.
body = re.sub(r"<blockquote>\s*<p><strong>(Definition|Worked example|Example|Caution|Defect|Open question|Note|Key idea)([^<]*)</strong>",
              lambda m: f'<blockquote class="callout c-{m.group(1).lower().split()[0]}"><p><strong>{m.group(1)}{m.group(2)}</strong>', body)
# Wrap tables for horizontal overflow
body = body.replace("<table>", '<div class="tw"><table>').replace("</table>", "</table></div>")

HEAD = """<title>Pangenome Privacy Primer</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans:wght@400;500;600;700&family=Source+Serif+4:opsz,wght@8..60,400;8..60,600&display=swap">
<style>
:root{
  --paper:#F4F6F7; --surface:#FFFFFF; --ink:#12171C; --ink-2:#4A5560; --ink-3:#6E7A86;
  --rule:#DCE2E7; --rule-2:#EDF1F4;
  --accent:#17587F; --accent-soft:#E4EEF5;
  --warn:#9C4E15; --warn-soft:#F7ECE2;
  --ok:#2C6349; --ok-soft:#E4F0EA;
  --code-bg:#EEF2F5;
  --sans:"IBM Plex Sans",ui-sans-serif,system-ui,sans-serif;
  --serif:"Source Serif 4",Georgia,"Times New Roman",serif;
  --mono:"IBM Plex Mono",ui-monospace,SFMono-Regular,Menlo,monospace;
}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){
  --paper:#0F1317; --surface:#161B21; --ink:#E7ECF1; --ink-2:#AFBAC5; --ink-3:#8B97A3;
  --rule:#28313A; --rule-2:#1E252C;
  --accent:#71B4DE; --accent-soft:#162833;
  --warn:#D4884A; --warn-soft:#2A1F16;
  --ok:#6FB694; --ok-soft:#152720;
  --code-bg:#1B222A;
}}
:root[data-theme="dark"]{
  --paper:#0F1317; --surface:#161B21; --ink:#E7ECF1; --ink-2:#AFBAC5; --ink-3:#8B97A3;
  --rule:#28313A; --rule-2:#1E252C;
  --accent:#71B4DE; --accent-soft:#162833;
  --warn:#D4884A; --warn-soft:#2A1F16;
  --ok:#6FB694; --ok-soft:#152720;
  --code-bg:#1B222A;
}
body{background:var(--paper);color:var(--ink);font-family:var(--serif);font-size:17px;line-height:1.62;}
.wrap{max-width:1180px;margin:0 auto;padding-block:0 4rem;padding-left:20px;padding-right:20px;}
.masthead{border-bottom:2px solid var(--ink);padding-block:2.6rem 1.2rem;margin-bottom:2rem;}
.eyebrow{font-family:var(--mono);font-size:.7rem;letter-spacing:.14em;text-transform:uppercase;color:var(--accent);margin:0 0 .9rem;}
h1{font-family:var(--sans);font-weight:700;font-size:clamp(2rem,5vw,3.1rem);line-height:1.06;letter-spacing:-.02em;margin:0 0 .7rem;text-wrap:balance;}
.standfirst{font-size:1.12rem;color:var(--ink-2);max-width:60ch;margin:0 0 1.4rem;}
.meta{display:flex;flex-wrap:wrap;gap:.5rem 1.6rem;font-family:var(--mono);font-size:.72rem;color:var(--ink-3);text-transform:uppercase;letter-spacing:.07em;}
.hero-fig{margin:1.8rem 0 .4rem;}
.hero-fig figcaption{font-family:var(--mono);font-size:.7rem;color:var(--ink-3);margin-top:.5rem;letter-spacing:.03em;}
.layout{display:grid;grid-template-columns:minmax(0,1fr);gap:2.5rem;}
@media(min-width:1000px){.layout{grid-template-columns:250px minmax(0,1fr);gap:3.5rem;}
  .toc{position:sticky;top:calc(env(safe-area-inset-top,0px) + 1.2rem);max-height:86vh;overflow-y:auto;align-self:start;}}
.toc{font-family:var(--sans);font-size:.83rem;line-height:1.4;border-top:1px solid var(--rule);padding-top:1rem;}
.toc-h{font-family:var(--mono);font-size:.66rem;letter-spacing:.14em;text-transform:uppercase;color:var(--ink-3);margin:0 0 .8rem;}
.toc ol{list-style:none;margin:0;padding:0;}
.toc-l1>li{margin-bottom:.75rem;}
.toc-l1>li>a{font-weight:600;color:var(--ink);}
.toc-l2{margin:.35rem 0 0 .1rem;border-left:1px solid var(--rule);padding-left:.7rem;}
.toc-l2 li{margin:.18rem 0;}
.toc-l2 a{color:var(--ink-2);}
.toc a{text-decoration:none;display:block;}
.toc a:hover{color:var(--accent);}
.prose{min-width:0;}
.sec{border-top:1px solid var(--rule);padding-top:2.4rem;margin-top:3rem;}
.sec:first-child{border-top:0;margin-top:0;padding-top:0;}
h2{font-family:var(--sans);font-weight:700;font-size:clamp(1.55rem,3.4vw,2.15rem);line-height:1.14;letter-spacing:-.015em;margin:0 0 1.3rem;text-wrap:balance;}
h3{font-family:var(--sans);font-weight:600;font-size:1.22rem;line-height:1.26;margin:2.4rem 0 .7rem;color:var(--ink);text-wrap:balance;}
h4{font-family:var(--sans);font-weight:600;font-size:1.02rem;margin:1.7rem 0 .5rem;color:var(--ink-2);}
p,li{max-width:68ch;}
p{margin:0 0 1.05rem;}
ul,ol{margin:0 0 1.15rem;padding-left:1.3rem;}
li{margin:.34rem 0;}
strong{font-weight:600;}
a{color:var(--accent);}
code{font-family:var(--mono);font-size:.855em;background:var(--code-bg);padding:.12em .36em;border-radius:3px;word-break:break-word;}
pre{background:var(--code-bg);border:1px solid var(--rule);border-radius:6px;padding:1rem 1.1rem;overflow-x:auto;margin:0 0 1.2rem;}
pre code{background:none;padding:0;font-size:.82rem;line-height:1.55;}
.tw{overflow-x:auto;margin:0 0 1.4rem;border:1px solid var(--rule);border-radius:6px;background:var(--surface);}
table{border-collapse:collapse;width:100%;font-family:var(--sans);font-size:.85rem;min-width:520px;}
th,td{text-align:left;padding:.62rem .8rem;border-bottom:1px solid var(--rule-2);vertical-align:top;}
th{font-weight:600;font-size:.72rem;letter-spacing:.06em;text-transform:uppercase;color:var(--ink-2);background:var(--rule-2);white-space:nowrap;}
tbody tr:last-child td{border-bottom:0;}
td code{font-size:.8em;}
blockquote{margin:0 0 1.3rem;padding:.9rem 1.1rem;border-left:3px solid var(--accent);background:var(--accent-soft);border-radius:0 5px 5px 0;}
blockquote p:last-child{margin-bottom:0;}
blockquote.c-defect,blockquote.c-caution{border-left-color:var(--warn);background:var(--warn-soft);}
blockquote.c-worked,blockquote.c-example{border-left-color:var(--ok);background:var(--ok-soft);}
hr{border:0;border-top:1px solid var(--rule);margin:2.4rem 0;}
h3[id^="in-plain-words"],h3:where(:not(:first-child)){scroll-margin-top:1.5rem;}
.sec h3{scroll-margin-top:1.5rem;}
@media(prefers-reduced-motion:no-preference){html{scroll-behavior:smooth;}}
</style>"""

HERO = """<figure class="hero-fig">
<svg viewBox="0 0 700 118" width="100%" role="img" aria-label="A variation graph: shared sequence nodes with two alternative alleles at a variant site, and one individual's path highlighted through it.">
  <defs><marker id="ah" markerWidth="7" markerHeight="7" refX="6" refY="2.4" orient="auto"><path d="M0,0 L0,4.8 L6,2.4 z" fill="currentColor"/></marker></defs>
  <g font-family="IBM Plex Mono, monospace" font-size="12" color="var(--ink-3)">
    <g stroke="var(--rule)" stroke-width="7" fill="none" stroke-linecap="round">
      <path d="M40,59 H150"/><path d="M196,59 C225,59 232,32 262,32"/><path d="M196,59 C225,59 232,86 262,86"/>
      <path d="M318,32 C348,32 352,59 384,59"/><path d="M318,86 C348,86 352,59 384,59"/><path d="M384,59 H494"/>
      <path d="M540,59 C568,59 574,32 604,32"/><path d="M540,59 C568,59 574,86 604,86"/>
    </g>
    <g stroke="var(--accent)" stroke-width="3.4" fill="none" stroke-linecap="round">
      <path d="M40,59 H150"/><path d="M196,59 C225,59 232,32 262,32"/><path d="M318,32 C348,32 352,59 384,59"/>
      <path d="M384,59 H494"/><path d="M540,59 C568,59 574,86 604,86" marker-end="url(#ah)"/>
    </g>
    <g fill="var(--surface)" stroke="var(--ink-2)" stroke-width="1.2">
      <rect x="150" y="45" width="46" height="28" rx="5"/><rect x="262" y="18" width="56" height="28" rx="5"/>
      <rect x="262" y="72" width="56" height="28" rx="5"/><rect x="494" y="45" width="46" height="28" rx="5"/>
      <rect x="604" y="18" width="52" height="28" rx="5"/><rect x="604" y="72" width="52" height="28" rx="5"/>
    </g>
    <g fill="var(--ink)" text-anchor="middle" font-size="12.5">
      <text x="173" y="64">GAT</text><text x="290" y="37">C</text><text x="290" y="91">T</text>
      <text x="517" y="64">CCA</text><text x="630" y="37">A</text><text x="630" y="91">AG</text>
    </g>
    <g fill="var(--ink-3)" font-size="9.5" text-anchor="middle">
      <text x="290" y="11">allele 0 (REF)</text><text x="290" y="112">allele 1 (ALT)</text>
    </g>
  </g>
</svg>
<figcaption>A variation graph. Shared sequence sits in nodes; a variant site branches into alternative alleles; one person's genome is a <em>path</em> (highlighted). Releasing that path is what this project is about.</figcaption>
</figure>"""

page = f"""{HEAD}
<div class="wrap">
<header class="masthead">
  <p class="eyebrow">RanPanMixer · Orientation</p>
  <h1>Releasing a Genome Without Revealing the Person</h1>
  <p class="standfirst">A ground-up primer on pangenome graphs, genomic privacy attacks, the PanMixer tool, and the randomized path-release mechanism this project is building. Written for a reader with basic biology and basic computer science, and nothing else.</p>
  <div class="meta"><span>Four sections</span><span>~53,600 words</span><span>Compiled 18 Sep 2026</span></div>
  {HERO}
</header>
<div class="layout">
{nav}
<main class="prose">
{body}
</main>
</div>
</div>"""
OUT.write_text(page)
print(f"wrote {OUT}  {len(page):,} bytes")
