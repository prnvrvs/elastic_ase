"""Energy-strain curves with their fits."""
import numpy as np


def plot_fits(fits, coefficients, filename="energy_strain.png", show=False):
    """One panel per deformation: (E - E0)/V0 in GPa against the strain in %, points and the
    fitted parabola. ``fits`` and ``coefficients`` as in ElasticResults. Returns the file name."""
    import matplotlib
    if not show:
        matplotlib.use("Agg", force=False)
    import matplotlib.pyplot as plt

    names = list(fits)
    ncols = min(3, len(names))
    nrows = int(np.ceil(len(names) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.2 * ncols, 3.4 * nrows), squeeze=False)
    for ax, name in zip(axes.flat, names):
        fit = fits[name]
        x = np.array(fit["strains"])
        dense = np.linspace(x.min(), x.max(), 200)
        ax.plot(dense * 100, np.polyval(fit["coefficients_GPa"], dense), "-", color="C0", lw=1.5)
        ax.plot(x * 100, fit["y_GPa"], "o", color="C3", ms=5)
        ax.set_title(f"{name}\na = {coefficients[name]:.2f} GPa,  R² = {fit['r2']:.5f}", fontsize=9)
        ax.set_xlabel("strain δ (%)")
        ax.set_ylabel("(E - E0)/V0 (GPa)")
        ax.grid(alpha=0.3)
    for ax in list(axes.flat)[len(names):]:
        ax.set_visible(False)
    fig.tight_layout()
    if filename:
        fig.savefig(filename, dpi=150)
    if show:
        plt.show()
    plt.close(fig)
    return filename
