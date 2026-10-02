# make_figures4.py —— Figure 4: INSPIRE 验证层 (4A 事件时间累计曲线 + 4B MELD梯度)
import os, csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "Figure4"; os.makedirs(OUT, exist_ok=True)
A = os.path.expanduser("~/nurse_faers/analysis")  # adjust to your environment
BLUE, ORANGE, GRAY = "#0072B2", "#D55E00", "#666666"

hist, meld, tim = [], [], {}
for r in csv.DictReader(open(f"{A}/insp_timing_hist.csv")):
    hist.append((r["liver"].upper() == "TRUE", int(r["day_bucket"]), int(r["N"])))
for r in csv.DictReader(open(f"{A}/insp_meld.csv")):
    meld.append((r["liver"].upper() == "TRUE", r["meld_grp"], int(r["n"]), int(r["ev"])))
for r in csv.DictReader(open(f"{A}/insp_timing.csv")):
    tim[r["liver"].upper() == "TRUE"] = (int(r["n"]), float(r["median_h"]),
                                 int(r["p_le24h"].rstrip("%")),
                                 int(r["p_le7d"].rstrip("%")))

fig = plt.figure(figsize=(7.0, 2.9))
gs = fig.add_gridspec(1, 2, width_ratios=[1.55, 1], wspace=0.34,
                      left=0.09, right=0.985, top=0.82, bottom=0.24)

## ---- 4A 累计曲线 ----
ax = fig.add_subplot(gs[0])
def cdf(liver_flag, color, lw, ls, label):
    rows = sorted([(d, n) for lv, d, n in hist if lv == liver_flag])
    tot = sum(n for _, n in rows); cum = 0; xs, ys = [0], [0]
    for d, n in rows:
        cum += n; xs.append(d + 1); ys.append(100 * cum / tot)
    ax.step(xs, ys, where="post", color=color, lw=lw, ls=ls, label=label)
    return xs, ys
cdf(False, BLUE, 1.6, "-", f"Non-liver (n={tim[False][0]})")
cdf(True, GRAY, 1.0, "--", f"Liver disease (n={tim[True][0]})")
ax.axvline(1, color=GRAY, lw=0.5, ls=":")
ax.axvline(7, color=GRAY, lw=0.5, ls=":")
ax.text(1, 96, "33% ≤24 h", ha="center", fontsize=6, color="#444444")
ax.text(7, 88, "70% ≤7 d", ha="center", fontsize=6, color="#444444")
ax.scatter([2.5], [0], marker="D", s=14, color=ORANGE, zorder=5, clip_on=False)
ax.annotate("FAERS median 2.5 d", xy=(2.5, 0), xytext=(4.2, 16),
            fontsize=6, color=ORANGE, ha="left",
            arrowprops=dict(arrowstyle="-", color=ORANGE, lw=0.6))
ax.set_xlim(0, 31); ax.set_ylim(0, 100)
ax.set_xlabel("Days from first postoperative opioid dose", fontsize=7)
ax.set_ylabel("Cumulative proportion of events (%)", fontsize=7)
ax.tick_params(labelsize=6.5)
ax.legend(fontsize=6, frameon=False, loc="lower right")
ax.set_title("A  Time to ward GCS ≤12 event (INSPIRE)", loc="left", fontsize=8)

## ---- 4B MELD 梯度 ----
ax2 = fig.add_subplot(gs[1])
order = ["<10", "10-14", "15-19", ">=20"]
lab = {"<10": "<10", "10-14": "10–14", "15-19": "15–19", ">=20": "≥20"}
mm = {g: (n, ev) for lv, g, n, ev in meld if lv is False}
xs = range(len(order))
vals = [100 * mm[g][1] / mm[g][0] for g in order]
ax2.bar(xs, vals, width=0.62, color=BLUE, edgecolor="white", linewidth=0.5)
for i, g in enumerate(order):
    ax2.text(i, vals[i] + 0.35, f"{vals[i]:.1f}%", ha="center", fontsize=6.5,
             color="#222222")
ax2.axhline(3.4, color=ORANGE, lw=0.8, ls="--")
ax2.text(3.42, 3.9, "overall 3.4%", fontsize=6, color=ORANGE, ha="right")
ax2.set_xticks(list(xs))
ax2.set_xticklabels([f"{lab[g]}\n(n={mm[g][0]})" for g in order], fontsize=6.3)
ax2.set_ylim(0, 12); ax2.set_xlim(-0.55, 3.55)
ax2.set_xlabel("MELD-Na stratum", fontsize=7)
ax2.set_ylabel("GCS ≤12 event rate (%)", fontsize=7)
ax2.tick_params(labelsize=6.5)
ax2.set_title("B  Rate by liver function (non-liver patients)",
              loc="left", fontsize=8)

fig.savefig(f"{OUT}/Figure4_validation.pdf", bbox_inches="tight")
import shutil
shutil.copy(f"{OUT}/Figure4_validation.pdf", f"{OUT}/Figure4_overview_preview.pdf")
print("done", os.listdir(OUT))
