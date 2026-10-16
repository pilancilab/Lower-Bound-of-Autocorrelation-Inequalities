"""Plot rigorous lower-bound progress on C_{6.2} across versions."""
import matplotlib.pyplot as plt
import pandas as pd

TEAL = "#1f9d8a"
TEAL_LIGHT = "#8ed1c7"
RED = "#d62728"

df = pd.read_csv("objectives.csv")

general_valid = df["status"].isin(["VALID", "VALID-vacuous", "VALID-empty", "VALID-tautology"])
bad = df["status"].isin(["INVALID", "INVALID-tightening"])

g = df[general_valid & df["Omega"].notna()].sort_values("version_num").copy()
g["running_best"] = g["Omega"].cummax()
b = df[bad].copy()

fig, ax = plt.subplots(figsize=(11, 5.5))

ax.plot(
    g["version_num"], g["running_best"],
    color=TEAL, linewidth=2.2, label="Rigorous SOTA (running best)", zorder=2,
)
ax.scatter(
    g["version_num"], g["Omega"],
    color=TEAL, s=50, edgecolor="white", linewidth=0.8,
    label=r"Valid $\Omega^{(k)}$", zorder=3,
)

ymin = g["Omega"].min()
if len(b):
    ax.scatter(
        b["version_num"], [ymin - 0.02] * len(b),
        marker="x", color=RED, s=75, linewidth=2.2,
        label="Invalid", zorder=3,
    )

best = g.loc[g["Omega"].idxmax()]
ax.scatter(
    [best.version_num], [best.Omega],
    marker="*", s=280, facecolor="orange", edgecolor="orange", linewidth=0,
    label=f"Best: v{int(best.version_num)} = {best.Omega:.4f}", zorder=4,
)
ax.annotate(
    f"v{int(best.version_num)}: {best.Omega:.4f}",
    xy=(best.version_num, best.Omega),
    xytext=(-55, -35), textcoords="offset points",
    ha="center", va="top",
    fontsize=10, color=TEAL, fontweight="bold",
    arrowprops=dict(arrowstyle="->", color=TEAL, lw=1.2),
)

ax.axhline(1.2802, color="black", linestyle=":", linewidth=1.4, label="Previous SOTA = 1.2802", zorder=1)

ax.set_xlabel("Iterations", fontsize=12)
ax.set_ylabel(r"Primal objective $\Omega^{(k)}$  (lower bound on $C_{6.2}$)", fontsize=12)
ax.set_title(r"Rigorous lower bounds on $C_{6.2}$ across iterations", fontsize=13)
ax.grid(True, alpha=0.3)
ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0, fontsize=9, frameon=True)
ax.set_xlim(-1, df["version_num"].max() + 2)

fig.tight_layout()
fig.savefig("objectives.png", dpi=160)
fig.savefig("objectives.pdf")
print("Saved objectives.png and objectives.pdf")
