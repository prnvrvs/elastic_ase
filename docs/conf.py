# Configuration file for the Sphinx documentation builder.
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import elastic_ase  # noqa: E402

# -- Project information -------------------------------------------------------

project = "elastic_ase"
copyright = "2026, Pranav Kumar"
author = "Pranav Kumar"
release = elastic_ase.__version__
version = release

# -- General configuration -----------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx_copybutton",
]
# copy button: leave out prompts such as ">>> " and "$ " when copying
copybutton_prompt_text = r">>> |\.\.\. |\$ "
copybutton_prompt_is_regexp = True
# links to the documentation of types used in signatures (ase.Atoms, numpy.ndarray, ...)
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "ase": ("https://docs.ase-lib.org", None),
}
# numpy-style type phrases that are not objects
nitpick_ignore = [("py:class", "optional"), ("py:class", "ASE optimizer class")]
templates_path = ["_templates"]
exclude_patterns = ["_build", "_generated", "Thumbs.db", ".DS_Store"]  # _generated is included, not a page
default_role = "code"
autodoc_default_options = {"members": True, "member-order": "bysource"}
autodoc_mock_imports = ["matplotlib", "elastic_vis", "tensorpotential"]

# -- HTML output ---------------------------------------------------------------

html_theme = "sphinx_book_theme"
html_title = "elastic_ase {}".format(release)
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_theme_options = {
    "repository_url": "https://github.com/prnvrvs/elastic_ase",
    "use_repository_button": True,
    "use_issues_button": True,
    "path_to_docs": "docs",
    "show_toc_level": 2,
}


# -- Tables generated from the code ---------------------------------------------
# The deformation tables come straight from elastic_ase.systems, so the
# documentation cannot drift from what the calculations actually do.

def _cij(name):
    """'C11' -> 'C_{11}'."""
    return "C_{%s}" % name[1:]


def _frac(value):
    """LaTeX for a small rational coefficient (1/2, 3/2, 2, ...)."""
    from fractions import Fraction

    f = Fraction(value).limit_denominator(12)
    if f.denominator == 1:
        return str(f.numerator)
    return r"\tfrac{%d}{%d}" % (f.numerator, f.denominator)


def _combination(terms):
    """LaTeX for sum(coefficient * constant); factors out a common coefficient."""
    terms = [(c, n) for c, n in terms if abs(c) > 1e-9]
    if not terms:
        return "0"
    common = abs(terms[0][0])
    if len(terms) > 1 and abs(common - 1) > 1e-9 and all(abs(abs(c) - common) < 1e-9 for c, _ in terms):
        inner = _combination([(c / common, n) for c, n in terms])
        return r"%s\,(%s)" % (_frac(common), inner)
    out = ""
    for k, (c, n) in enumerate(terms):
        mag = abs(c)
        body = _cij(n) if abs(mag - 1) < 1e-9 else r"%s\,%s" % (_frac(mag), _cij(n))
        if k == 0:
            out = ("-" if c < 0 else "") + body
        else:
            out += (" - " if c < 0 else " + ") + body
    return out


def _F(d):
    """LaTeX pmatrix of the deformation F = I + eps applied to the cell, in terms of delta.

    Diagonal: 1 + e_i delta (+ the volume correction); off-diagonal: e_k delta / 2 for the
    engineering shear strains (e4 -> yz, e5 -> xz, e6 -> xy), exactly as in ElasticASE.runs().
    """
    shear = {(0, 1): 5, (0, 2): 4, (1, 2): 3}
    rows = []
    for i in range(3):
        cells = []
        for j in range(3):
            if i == j:
                text = "1"
                c = d.strain[i]
                if c:
                    text += r" %s \delta" % ("+" if c > 0 else "-")
                if d.correction is not None and d.correction[0] == i:
                    text += r" + \frac{\delta^2}{%d - \delta^2}" % d.correction[1]
            else:
                c = d.strain[shear[(min(i, j), max(i, j))]]
                text = "0" if not c else (r"\tfrac{\delta}{2}" if c > 0 else r"-\tfrac{\delta}{2}")
            cells.append(text)
        rows.append(" & ".join(cells))
    return r"\begin{pmatrix}" + r" \\ ".join(rows) + r"\end{pmatrix}"


def _pattern(spec):
    """LaTeX pmatrix of the stiffness-matrix layout in terms of the independent constants."""
    import numpy as np

    n = len(spec.constants)
    units = [spec.matrix(np.eye(n)[j]) for j in range(n)]
    rows = []
    for i in range(6):
        cells = []
        for j in range(6):
            cells.append(_combination([(u[i, j], c) for u, c in zip(units, spec.constants)]))
        rows.append(" & ".join(cells))
    return "\\begin{pmatrix}\n" + " \\\\\n".join(rows) + "\n\\end{pmatrix}"


def _generate_system_tables(app):
    from elastic_ase.systems import SYSTEMS

    out = ROOT / "docs" / "_generated"
    out.mkdir(exist_ok=True)
    titles = {
        "cubic": "Cubic", "hexagonal": "Hexagonal", "tetragonal_1": "Tetragonal I",
        "tetragonal_2": "Tetragonal II", "trigonal_1": "Trigonal I", "trigonal_2": "Trigonal II",
        "orthogonal": "Orthorhombic", "monoclinic": "Monoclinic", "triclinic": "Triclinic",
    }
    order = ["cubic", "hexagonal", "tetragonal_1", "tetragonal_2", "trigonal_1",
             "trigonal_2", "orthogonal", "monoclinic", "triclinic"]
    parts = []
    for name in order:
        spec = SYSTEMS[name]
        design = spec.design_matrix()
        title = titles[name]
        parts += [title, "~" * len(title), "",
                  '``symmetry="{}"`` – {} independent constants: :math:`{}`'.format(
                      name, len(spec.constants), r",\ ".join(_cij(c) for c in spec.constants)), ""]
        if name != "triclinic":  # the triclinic matrix is just all 21 entries
            parts += [".. math::", "", "   C = " + _pattern(spec).replace("\n", "\n   "), ""]
        parts += [".. list-table::", "   :header-rows: 1", "   :widths: 20 48 32", "",
                  "   * - Run folder",
                  "     - Deformation :math:`F = I + \\varepsilon` applied to the cell",
                  "     - Fitted :math:`a` in :math:`(E - E_0)/V_0 = a\\,\\delta^2 + \\dots`"]
        for k, d in enumerate(spec.deformations):
            a = _combination(list(zip(design[k], spec.constants)))
            parts += ["   * - ``{}``".format(d.folder),
                      "     - :math:`F = {}`".format(_F(d)),
                      "     - :math:`a = {}`".format(a)]
        parts.append("")
    (out / "system_tables.rst").write_text("\n".join(parts) + "\n")


def setup(app):
    app.connect("builder-inited", _generate_system_tables)

# The ASME and APS DOIs are valid but the publishers answer automated requests with 403.
linkcheck_ignore = [r"https://doi\.org/10\.1115/1\.4055097", r"https://doi\.org/10\.1103/"]  # APS too
