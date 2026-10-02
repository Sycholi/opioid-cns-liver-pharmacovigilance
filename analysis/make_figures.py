# make_figures.py —— 双库研究主图集(3图, 面板级PDF + overview预览)
# 图1 数据流+总体信号 / 图2 分药梯度森林 / 图3 时间谱+严重结局
# 配色: Okabe-Ito 蓝#0072B2(全库) 朱#D55E00(肝病层) 绿#009E73(JADER)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd, os

BASE = os.path.expanduser("~/nurse_faers/analysis")
OUT1, OUT2, OUT3 = "Figure1", "Figure2", "Figure3"
for d in (OUT1, OUT2, OUT3): os.makedirs(d, exist_ok=True)
plt.rcParams.update({"font.size": 7.5, "axes.linewidth": 0.6,
                     "xtick.major.width": 0.5, "ytick.major.width": 0.5,
                     "font.family": "DejaVu Sans"})
C_ALL, C_LIV, C_JAD = "#0072B2", "#D55E00", "#009E73"

ror = pd.read_csv(f"{BASE}/faers_ror.csv")
jor = pd.read_csv(f"{BASE}/jader_ror.csv")
ttf = pd.read_csv(f"{BASE}/faers_tto.csv")
ttj = pd.read_csv(f"{BASE}/jader_tto.csv")
ocf = pd.read_csv(f"{BASE}/faers2_outcome_reports.csv").rename(columns={"ocg":"oc"})
ocl = pd.read_csv(f"{BASE}/faers2_outcome_liver_reports.csv").rename(columns={"ocg":"oc"})
_r2 = pd.read_csv(f"{BASE}/faers2_liver_ror.csv")
_r2 = _r2[_r2["def"] == "primary_union"].drop(columns=["def"]).assign(layer="肝病层", outcome="CNS5")
ror = pd.concat([ror[~((ror.layer == "肝病层") & (ror.outcome == "CNS5"))], _r2], join="inner").reset_index(drop=True)
ttf = pd.DataFrame([
    dict(n=4145, shape=0.4233, scale=5.894, median=2.480, p_le7d=0.786, p_le30d=0.881),
    dict(n=36,   shape=0.4431, scale=5.590, median=2.445, p_le7d=0.750, p_le30d=0.889)])
ttj = pd.DataFrame([dict(n=108, shape=0.5403, scale=19.570, median=9.932, p_le7d=0.537, p_le30d=0.759)])

# ---------- 图1A 数据流 ----------
fig, ax = plt.subplots(figsize=(7.2, 2.5))
ax.axis("off")
def box(x, y, w, h, text, fc="#ffffff", ec="#333333", fs=6.8, lw=0.6):
    ax.add_patch(plt.Rectangle((x, y), w, h, fc=fc, ec=ec, lw=lw, zorder=2))
    ax.text(x + w/2, y + h/2, text, ha="center", va="center", fontsize=fs, zorder=3)
def arrow(x1, y1, x2, y2):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="-|>", color="#555555", lw=0.7))
box(0.5, 0.62, 2.3, 0.62, "FAERS 2004Q1–2026Q2\n90 quarterly files", fc="#e8f1f8")
box(3.6, 0.62, 2.3, 0.62, "JADER (PMDA) 202608\n1,049,684 reports", fc="#e9f6f0")
box(0.5, 0.10, 2.3, 0.34, "20,692,684 reports after\ndeduplication", fs=6.2)
box(3.6, 0.10, 2.3, 0.34, "aligned opioid-suspected\nreports (n = 9,654)", fs=6.2)
arrow(1.65, 0.62, 1.65, 0.44); arrow(4.75, 0.62, 4.75, 0.44)
box(6.55, 0.62, 1.15, 0.34, "CNS-depression\nPT set", fs=6.2)
arrow(5.9, 0.78, 6.55, 0.78)
ax.text(1.65, 1.36, "Primary-suspect opioid reports (n = 340,888)\n× CNS-depression events (n = 16,141)",
        ha="center", fontsize=6.6, color=C_ALL)
ax.text(1.65, -0.16, "Liver population (n = 163,023):\nstrict liver indications ∪ liver PTs",
        ha="center", fontsize=6.2, color=C_LIV)
ax.text(4.75, 1.36, "Liver population by\nunderlying disease\n(n = 22,570)", ha="center",
        fontsize=6.2, color=C_LIV)
ax.plot([0.3, 7.85], [1.02, 1.02], lw=0, ls="")
ax.set_xlim(0.3, 7.85); ax.set_ylim(-0.3, 1.62)
fig.savefig(f"{OUT1}/Figure1A_flow.pdf", bbox_inches="tight")
plt.close(fig)

# ---------- 图1B 总体信号森林 ----------
rows = [
    ("FAERS  overall",  ror.query("layer=='全库' and drug=='全阿片类' and outcome=='CNS5'").iloc[0], C_ALL),
    ("FAERS  liver",    ror.query("layer=='肝病层' and drug=='全阿片类' and outcome=='CNS5'").iloc[0], C_LIV),
    ("JADER  overall",  jor.query("layer=='全库' and drug=='全阿片类(对齐集)' and outcome=='CNS5'").iloc[0], C_JAD),
]
fig, ax = plt.subplots(figsize=(3.6, 1.75))
for i, (lab, r, c) in enumerate(rows):
    y = len(rows) - 1 - i
    ax.plot([r.lo, r.hi], [y, y], color=c, lw=1.6, solid_capstyle="round", zorder=2)
    ax.scatter([r.ROR], [y], s=14, color=c, edgecolors="white", linewidths=0.4, zorder=3)
    ax.text(9.6, y, f"{r.ROR:.2f} ({r.lo:.2f}–{r.hi:.2f})", va="center", fontsize=6.8)
    ax.text(-0.03, y, f"{lab}  (n={int(r.a):,})", transform=ax.get_yaxis_transform(),
            ha="right", va="center", fontsize=6.8)
ax.axvline(1, color="#888888", lw=0.6, ls="--")
ax.set_xscale("log"); ax.set_xlim(0.9, 9.5)
ax.set_xticks([1, 2, 3, 5, 8]); ax.set_xticklabels(["1", "2", "3", "5", "8"])
ax.minorticks_off()
ax.grid(axis="x", which="major", color="#dddddd", lw=0.4)
ax.set_yticks([]); ax.set_xlabel("ROR (95% CI), log scale", fontsize=7)
for s in ("left", "top", "right"): ax.spines[s].set_visible(False)
ax.spines["bottom"].set_bounds(1.8, 9.0)
fig.savefig(f"{OUT1}/Figure1B_overall_forest.pdf", bbox_inches="tight")
plt.close(fig)

# ---------- 图2 分药梯度森林 ----------
order = ["tramadol", "morphine", "hydromorphone", "oxycodone", "fentanyl", "hydrocodone"]
cn = {"tramadol": "Tramadol", "morphine": "Morphine", "hydromorphone": "Hydromorphone",
      "oxycodone": "Oxycodone", "fentanyl": "Fentanyl", "hydrocodone": "Hydrocodone",
      "all": "All opioids"}
tag = {"tramadol": "not recommended (Child-Pugh C)",
       "morphine": "reduce (start low)", "hydromorphone": "reduce 1/4–1/2",
       "oxycodone": "reduce dose", "fentanyl": "avoid severe (transdermal)",
       "hydrocodone": "no adjustment (mild/mod)"}
fig, ax = plt.subplots(figsize=(7.2, 4.0))
ys, y = [], 0
for drug in order[::-1]:
    specs = [("FAERS overall", ror.query(f"layer=='全库' and drug=='{drug}'"), C_ALL, 0),
             ("FAERS liver", ror.query(f"layer=='肝病层' and drug=='{drug}'"), C_LIV, 1),
             ("JADER overall", jor.query(f"drug=='{drug}' and layer=='全库'"), C_JAD, 2)]
    for nm, sub, c, k in specs:
        if len(sub) == 0 or pd.isna(sub.iloc[0].ROR): continue
        r = sub.iloc[0]
        yo = y + (1 - k) * 0.28
        ax.plot([r.lo, r.hi], [yo, yo], color=c, lw=1.3, solid_capstyle="round")
        ax.scatter([r.ROR], [yo], s=16, color=c, zorder=3)
        ax.text(45, yo, f"{r.ROR:.2f} ({r.lo:.2f}–{r.hi:.2f})", va="center",
                fontsize=6.2, color=c)
    ax.text(-0.015, y + 0.52, f"{cn[drug]}", transform=ax.get_yaxis_transform(),
            ha="right", va="center", fontsize=7.5, fontweight="bold")
    ax.text(-0.015, y - 0.13, f"PI: {tag[drug]}", transform=ax.get_yaxis_transform(),
            ha="right", va="center", fontsize=5.8, color="#444444")
    y += 1.4
ax.axvline(1, color="#888888", lw=0.6, ls="--")
ax.set_xscale("log"); ax.set_xlim(0.55, 43)
ax.set_xticks([1, 2, 5, 10, 20, 40]); ax.set_xticklabels(["1", "2", "5", "10", "20", "40"])
ax.set_ylim(-0.55, y - 0.55); ax.set_yticks([])
ax.set_xlabel("ROR (95% CI) for CNS-depression events, log scale", fontsize=7)
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([], [], color=C_ALL, marker="o", ms=3.5, lw=1.2, label="FAERS overall"),
                   Line2D([], [], color=C_LIV, marker="o", ms=3.5, lw=1.2, label="FAERS liver population"),
                   Line2D([], [], color=C_JAD, marker="o", ms=3.5, lw=1.2, label="JADER overall")],
          loc="lower right", fontsize=6.5, frameon=False)
for s in ("left", "top", "right"): ax.spines[s].set_visible(False)
ax.spines["bottom"].set_bounds(0.55, 43)
fig.savefig(f"{OUT2}/Figure2A_drug_forest.pdf", bbox_inches="tight")
plt.close(fig)

# ---------- 图3A Weibull 时间谱 ----------
fig, ax = plt.subplots(figsize=(4.7, 3.3))
t = np.linspace(0.01, 90, 500)
def w_cdf(t, k, lam): return 1 - np.exp(-(t / lam) ** k)
groups = [("FAERS overall", ttf.iloc[0], C_ALL), ("FAERS liver", ttf.iloc[1], C_LIV),
          ("JADER overall", ttj.iloc[0], C_JAD)]
for nm, r, c in groups:
    ax.plot(t, 100 * w_cdf(t, r["shape"], r["scale"]), color=c, lw=1.5,
            label=f"{nm} (n={int(r['n'])}, median {r['median']:.1f} d)")
for xd in (7, 30):
    ax.axvline(xd, color="#aaaaaa", lw=0.5, ls=":")
    ax.text(xd, 104, f"{xd} d", ha="center", fontsize=6, color="#777777")
r1 = ttf.iloc[1]
ax.annotate(f"{100*r1['p_le7d']:.0f}% ≤ 7 d\n{100*r1['p_le30d']:.0f}% ≤ 30 d",
            xy=(7, 100 * w_cdf(7, r1["shape"], r1["scale"])),
            xytext=(22, 32), fontsize=6.5, color=C_LIV,
            arrowprops=dict(arrowstyle="-", color=C_LIV, lw=0.6))
ax.set_xlim(0, 90); ax.set_ylim(0, 100)
ax.set_xlabel("Days from opioid initiation to CNS-depression onset", fontsize=7)
ax.set_ylabel("Cumulative proportion of events (%)", fontsize=7)
ax.legend(fontsize=6.2, frameon=False, loc="lower right")
fig.savefig(f"{OUT3}/Figure3A_tto_weibull.pdf", bbox_inches="tight")
plt.close(fig)

# ---------- 图3B 严重结局构成 ----------
order_oc = ["住院", "死亡", "危及生命", "致残", "需干预", "先天异常"]
en = {"住院": "Hospitalisation", "死亡": "Death", "危及生命": "Life-threatening",
      "致残": "Disability", "需干预": "Required intervention", "先天异常": "Congenital anomaly"}
gf = ocf.set_index("oc"); gl = ocl.set_index("oc")
x = np.arange(len(order_oc)); w = 0.38
fig, ax = plt.subplots(figsize=(2.9, 3.3))
ax.bar(x - w/2, [100 * gf.loc[o, "pct"] / 100 for o in order_oc], w,
       color=C_ALL, label="FAERS overall")
ax.bar(x + w/2, [100 * gl.loc[o, "pct"] / 100 for o in order_oc], w,
       color=C_LIV, label="FAERS liver population")
for xi, o in zip(x, order_oc):
    ax.text(xi - w/2, gf.loc[o, "pct"] + 0.8, f"{gf.loc[o, 'pct']:.1f}",
            ha="center", fontsize=5.4, color=C_ALL)
    ax.text(xi + w/2, gl.loc[o, "pct"] + 0.8, f"{gl.loc[o, 'pct']:.1f}",
            ha="center", fontsize=5.4, color=C_LIV)
ax.set_xticks(x); ax.set_xticklabels([en[o] for o in order_oc], rotation=45,
                                     ha="right", fontsize=6)
ax.set_ylabel("Proportion of reports (%)", fontsize=7); ax.set_ylim(0, 88); ax.margins(x=0.06)
ax.legend(fontsize=6.2, frameon=False)
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.savefig(f"{OUT3}/Figure3B_outcome.pdf", bbox_inches="tight")
plt.close(fig)

# ---------- overview 预览(保持比例, 不拉伸) ----------
import subprocess
def preview(folder, files, out):
    try:
        from pypdf import PdfWriter
        wr = PdfWriter()
        for f in files: wr.append(f"{folder}/{f}")
        wr.write(f"{folder}/{out}")
    except ImportError:
        import subprocess
        subprocess.run(["pdfunite"] + [f"{folder}/{f}" for f in files] +
                       [f"{folder}/{out}"], check=True)
preview(OUT1, ["Figure1A_flow.pdf", "Figure1B_overall_forest.pdf"],
        "Figure1_overview_preview.pdf")
preview(OUT2, ["Figure2A_drug_forest.pdf"], "Figure2_overview_preview.pdf")
preview(OUT3, ["Figure3A_tto_weibull.pdf", "Figure3B_outcome.pdf"],
        "Figure3_overview_preview.pdf")
print("figures done:", os.listdir(OUT1), os.listdir(OUT2), os.listdir(OUT3))
