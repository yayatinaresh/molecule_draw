#!/usr/bin/env python
"""
Qt app: draw molecules from SMILES / PubChem names as aligned ACS 1996 SVGs.

Type or load molecules as "SMILES  name" lines, see live ACS 1996 previews
aligned on a common core, click a tile to copy its SVG (paste straight into
Inkscape), or export everything to a folder.

Self-contained: needs only rdkit and PySide6.

Usage:
    python scheme_gui.py [molecules.smi]
"""
import json
import os
import sys
import urllib.error
import urllib.parse
import urllib.request

from rdkit import Chem, RDLogger
from rdkit.Chem import rdDepictor, rdFMCS, rdCIPLabeler
from rdkit.Geometry import Point3D
from rdkit.Chem.Draw import rdMolDraw2D

import re
from PySide6.QtGui import QShortcut, QKeySequence
from PySide6.QtCore import Qt, QTimer, QByteArray, QMimeData, QSize, QThread, Signal
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QSplitter, QVBoxLayout, QHBoxLayout,
    QFormLayout, QPlainTextEdit, QLineEdit, QComboBox, QPushButton, QLabel,
    QScrollArea, QGridLayout, QFrame, QFileDialog, QDoubleSpinBox, QSpinBox,
    QMessageBox, QCheckBox,
)
from PySide6.QtSvgWidgets import QSvgWidget


def orient(mol, match, direction):
    """Rotate the 2D drawing so the core's heteroatom points straight
    down (or up) from the core centroid, i.e. the ring sits upright."""
    import math
    conf = mol.GetConformer()
    pts = [conf.GetAtomPosition(i) for i in match]
    cx = sum(p.x for p in pts) / len(pts)
    cy = sum(p.y for p in pts) / len(pts)
    het = [i for i in match if mol.GetAtomWithIdx(i).GetAtomicNum() != 6]
    if not het:
        return
    h = conf.GetAtomPosition(het[0])
    target = -math.pi / 2 if direction == "down" else math.pi / 2
    t = target - math.atan2(h.y - cy, h.x - cx)
    c, s = math.cos(t), math.sin(t)
    for i in range(mol.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        x, y = p.x - cx, p.y - cy
        conf.SetAtomPosition(i, Point3D(c * x - s * y + cx, s * x + c * y + cy, 0.0))

rdDepictor.SetPreferCoordGen(True)
RDLogger.DisableLog("rdApp.*")  # bad SMILES are reported in the tile instead

EXAMPLE = """# SMILES  name  [(degrees)]   e.g.  CCO ethanol (90)  rotates 90 deg counterclockwise
# lines starting with # are ignored
OCc1ccc(C=O)o1         HMF
OCc1ccc(CO)o1          BHMF
O=Cc1ccc(C=O)o1        DFF
OC(=O)c1ccc(C(=O)O)o1  FDCA
Cc1ccc(C)o1            DMF
"""


# ---------------------------------------------------------------- chemistry
_ROT = re.compile(r"^(.*?)\s*\(\s*([-+]?\d+(?:\.\d+)?)\s*(?:deg)?\s*\)\s*$", re.I)


def parse(text):
    """Return list of (name, smiles, angle) from 'SMILES  name  (degrees)'
    lines. The trailing (degrees) is optional: counterclockwise rotation of
    the drawing, e.g.  CCO ethanol (90)."""
    out = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(None, 1)
        smi = parts[0]
        rest = parts[1].strip() if len(parts) > 1 else ""
        angle = 0.0
        mt = _ROT.match(rest)
        if mt:
            rest, angle = mt.group(1).strip(), float(mt.group(2))
        out.append((rest or "mol%d" % len(out), smi, angle))
    return out


def make_template(core):
    core = core.strip()
    if not core or core.lower() == "auto":
        return None
    t = Chem.MolFromSmiles(core) or Chem.MolFromSmarts(core)
    if t is None:
        raise ValueError("Bad core: " + core)
    rdDepictor.Compute2DCoords(t)
    return t


def _is_linear(t):
    """True for an acyclic, unbranched, connected chain of >= 3 atoms."""
    try:
        nring = t.GetRingInfo().NumRings()
    except RuntimeError:  # SMARTS-derived: ring info not computed yet
        Chem.FastFindRings(t)
        nring = t.GetRingInfo().NumRings()
    if nring or t.GetNumAtoms() < 3:
        return False
    if any(a.GetDegree() > 2 for a in t.GetAtoms()):
        return False
    return len(Chem.GetMolFrags(t)) == 1


def _has_ring(t):
    try:
        return t.GetRingInfo().NumRings() > 0
    except RuntimeError:
        Chem.FastFindRings(t)
        return t.GetRingInfo().NumRings() > 0


def longest_chain(mol):
    """Atom indices of the longest heavy-atom chain of an acyclic molecule
    (ties: most carbons). None if not applicable (rings, <3 atoms, triple
    bonds or specified C=C stereo on the chain, where a plain zigzag could
    misdraw the geometry)."""
    import numpy as np
    if _has_ring(mol):
        return None
    dm = Chem.GetDistanceMatrix(mol)
    n = mol.GetNumAtoms()
    finite = np.where(dm < n, dm, 0)
    maxd = int(finite.max()) if n else 0
    if maxd < 2:
        return None
    best, best_c = None, -1
    for i, j in np.argwhere(np.triu(finite == maxd)):
        path = list(Chem.GetShortestPath(mol, int(i), int(j)))
        c = sum(1 for k in path if mol.GetAtomWithIdx(k).GetAtomicNum() == 6)
        if c > best_c:
            best, best_c = path, c
    prev_double = False
    for x, y in zip(best, best[1:]):
        bd = mol.GetBondBetweenAtoms(x, y)
        bt = bd.GetBondType()
        if bt == Chem.BondType.TRIPLE:
            return None
        is_d = bt == Chem.BondType.DOUBLE
        if is_d and (prev_double or bd.GetStereo() not in (
                Chem.BondStereo.STEREONONE, Chem.BondStereo.STEREOANY)):
            return None
        prev_double = is_d
    return best


def chain_template(mol):
    """(reference, pattern): a straight zigzag of this molecule's longest
    chain (including its terminal heteroatoms), or None."""
    path = longest_chain(mol)
    if path is None:
        return None
    rw = Chem.RWMol()
    for i in path:
        rw.AddAtom(Chem.Atom(mol.GetAtomWithIdx(i).GetAtomicNum()))
    for k in range(len(path) - 1):
        rw.AddBond(k, k + 1, Chem.BondType.SINGLE)
    ref = rw.GetMol()
    ref.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(ref)
    rdDepictor.Compute2DCoords(ref)
    patt = Chem.MolFromSmarts(
        "~".join("[#%d]" % mol.GetAtomWithIdx(i).GetAtomicNum() for i in path))
    return ref, patt


def auto_template(smiles_list, hetero="up"):
    """Largest common skeleton (elements kept, bond orders ignored) of all
    valid molecules, drawn on the first molecule that contains it.
    Returns (reference_mol, pattern) or None."""
    mols = [m for m in (Chem.MolFromSmiles(x) for x in smiles_list) if m is not None]
    if len(mols) < 2:
        return None
    res = rdFMCS.FindMCS(
        mols, timeout=5,
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareAny,
        ringMatchesRingOnly=True, completeRingsOnly=True)
    if res.numAtoms < 3:
        return None
    patt = Chem.MolFromSmarts(res.smartsString)
    ref = next((m for m in mols if m.HasSubstructMatch(patt)), None)
    if ref is None:
        return None
    ct = chain_template(ref)  # straighten the reference's own long chain
    if ct is not None:
        rdDepictor.GenerateDepictionMatching2DStructure(ref, ct[0], refPatt=ct[1])
    else:
        rdDepictor.Compute2DCoords(ref)
        if hetero != "free" and _has_ring(patt):
            orient(ref, ref.GetSubstructMatch(patt), hetero)
    return ref, patt


def drawing_ok(orig, m):
    """Does the 2D drawing still encode the SMILES' stereochemistry?
    (Catches a template forcing e.g. a Z double bond to be drawn trans.)"""
    try:
        mm = Chem.Mol(m)
        Chem.WedgeMolBonds(mm, mm.GetConformer())
        back = Chem.MolFromMolBlock(Chem.MolToMolBlock(mm))
        if back is None:
            return True
        for b in orig.GetBonds():  # unspecified C=C: ignore the geometry drawn
            if b.GetBondType() == Chem.BondType.DOUBLE and b.GetStereo() in (
                    Chem.BondStereo.STEREONONE, Chem.BondStereo.STEREOANY):
                back.GetBondWithIdx(b.GetIdx()).SetStereo(Chem.BondStereo.STEREONONE)
        return Chem.MolToSmiles(back) == Chem.MolToSmiles(orig)
    except Exception:  # noqa: BLE001
        return True


def level_ends(mol):
    """Rotate slightly so the first and last atom of the aligned chain sit at
    the same height (an even-length zigzag has its ends half a step apart)."""
    import math
    if not mol.HasProp("chain_idx") or mol.GetNumConformers() == 0:
        return
    idx = [int(x) for x in mol.GetProp("chain_idx").split(",")]
    conf = mol.GetConformer()
    p0, p1 = conf.GetAtomPosition(idx[0]), conf.GetAtomPosition(idx[-1])
    ang = math.atan2(p1.y - p0.y, p1.x - p0.x)
    if ang > math.pi / 2:
        ang -= math.pi
    elif ang < -math.pi / 2:
        ang += math.pi
    c, s = math.cos(-ang), math.sin(-ang)
    for i in range(mol.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, Point3D(c * p.x - s * p.y, s * p.x + c * p.y, 0.0))


def build(smi, template, hetero):
    """template: None | "chain" | (ref, patt) | core Mol (typed by the user)."""
    m = Chem.MolFromSmiles(smi)
    if m is None:
        raise ValueError("Invalid SMILES")
    orig = Chem.Mol(m)
    aligned, chain = False, None
    if isinstance(template, str):  # "chain": this molecule's own longest chain
        ct = chain_template(m)
        if ct is not None and m.HasSubstructMatch(ct[1]):
            rdDepictor.GenerateDepictionMatching2DStructure(m, ct[0], refPatt=ct[1])
            aligned, chain = True, m.GetSubstructMatch(ct[1])
        else:
            rdDepictor.Compute2DCoords(m)
    elif isinstance(template, tuple):  # common skeleton of the whole set
        ref, patt = template
        if m.HasSubstructMatch(patt):
            rdDepictor.GenerateDepictionMatching2DStructure(m, ref, refPatt=patt)
            aligned = True
            if _is_linear(patt):
                chain = m.GetSubstructMatch(patt)
        else:
            rdDepictor.Compute2DCoords(m)
    elif template is not None and m.HasSubstructMatch(template):  # typed core
        rdDepictor.GenerateDepictionMatching2DStructure(m, template)
        aligned = True
        if _has_ring(template):
            if hetero != "free":  # up/down rule only makes sense for rings
                orient(m, m.GetSubstructMatch(template), hetero)
        elif _is_linear(template):
            chain = m.GetSubstructMatch(template)
    else:
        rdDepictor.Compute2DCoords(m)
    if chain:
        m.SetProp("chain_idx", ",".join(str(i) for i in chain))
    if aligned and not drawing_ok(orig, m):
        m = Chem.Mol(orig)  # alignment would misdraw E/Z: fall back, and say so
        rdDepictor.Compute2DCoords(m)
        m.SetBoolProp("unaligned", True)
    return m


def flip_coords(mol, horizontal=False, vertical=False):
    """Mirror the 2D drawing. Wedges are recomputed from the chiral tags at
    draw time, so the molecule itself (its stereochemistry) is unchanged."""
    if not (horizontal or vertical) or mol.GetNumConformers() == 0:
        return
    conf = mol.GetConformer()
    for i in range(mol.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, Point3D(-p.x if horizontal else p.x,
                                        -p.y if vertical else p.y, 0.0))


def rotate_coords(mol, degrees):
    """Rotate the 2D drawing counterclockwise. Atom labels stay upright
    because they are drawn at render time, not rotated with the structure."""
    import math
    if not degrees or mol.GetNumConformers() == 0:
        return
    t = math.radians(degrees)
    c, s = math.cos(t), math.sin(t)
    conf = mol.GetConformer()
    for i in range(mol.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, Point3D(c * p.x - s * p.y, s * p.x + c * p.y, 0.0))


def stereo_badge(mol):
    """None if the molecule has no possible stereo elements, else
    (text, is_warning, tooltip). Only a warning, never blocks drawing."""
    try:
        infos = Chem.FindPotentialStereo(mol)
    except Exception:  # noqa: BLE001
        return None
    kinds = {Chem.StereoType.Atom_Tetrahedral: "stereocentres",
             Chem.StereoType.Bond_Double: "double bonds"}
    tot, dfn = {}, {}
    for si in infos:
        k = kinds.get(si.type)
        if k is None:
            continue
        tot[k] = tot.get(k, 0) + 1
        if si.specified == Chem.StereoSpecified.Specified:
            dfn[k] = dfn.get(k, 0) + 1
    n, d = sum(tot.values()), sum(dfn.values())
    if n == 0:
        return None
    detail = "; ".join("%s: %d of %d defined" % (k, dfn.get(k, 0), v)
                       for k, v in tot.items())
    tip = detail
    if d < n:
        tip += ("\nUndefined stereo is drawn without wedges. That is fine for a "
                "racemate or an unspecified isomer; otherwise use an isomeric "
                "SMILES (e.g. from PubChem).")
    return "stereo: %d of %d defined" % (d, n), d < n, tip


def annotate_cip(mol):
    """Write R/S and E/Z labels into the drawing (as atom/bond notes)."""
    rdCIPLabeler.AssignCIPLabels(mol)
    for a in mol.GetAtoms():
        if a.HasProp("_CIPCode"):
            a.SetProp("atomNote", a.GetProp("_CIPCode"))
    for b in mol.GetBonds():
        if b.HasProp("_CIPCode"):
            b.SetProp("bondNote", b.GetProp("_CIPCode"))


PUBCHEM = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/%s"
           "/property/IsomericSMILES/JSON")


def pubchem_smiles(name, timeout=10):
    """Isomeric SMILES for a compound name or CAS number via PubChem PUG REST."""
    url = PUBCHEM % urllib.parse.quote(name, safe="")
    req = urllib.request.Request(url, headers={"User-Agent": "scheme-gui"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            data = json.load(r)
    except urllib.error.HTTPError as e:
        raise ValueError("not found on PubChem" if e.code == 404 else "HTTP %d" % e.code)
    except urllib.error.URLError as e:
        raise ValueError("no connection (%s)" % e.reason)
    props = data["PropertyTable"]["Properties"][0]
    # PubChem renamed the field (IsomericSMILES -> SMILES); accept either
    for key in ("SMILES", "IsomericSMILES"):
        if props.get(key):
            return props[key]
    for k, v in props.items():
        if k.endswith("SMILES") and not k.startswith(("Connectivity", "Canonical")):
            return v
    raise ValueError("no SMILES in the response")


class LookupThread(QThread):
    """Runs the PubChem requests off the GUI thread."""
    got = Signal(list, list)  # [(name, smiles)], [error strings]

    def __init__(self, names):
        super().__init__()
        self.names = names

    def run(self):
        ok, bad = [], []
        for n in self.names:
            try:
                ok.append((n, pubchem_smiles(n)))
            except Exception as e:  # noqa: BLE001
                bad.append("%s: %s" % (n, e))
        self.got.emit(ok, bad)


def svg_text(mol):
    d = rdMolDraw2D.MolDraw2DSVG(-1, -1)
    rdMolDraw2D.DrawMoleculeACS1996(d, mol)
    d.FinishDrawing()
    txt = d.GetDrawingText()
    # drop RDKit's white background rectangle so SVGs are transparent
    return re.sub(r"<rect[^>]*fill:#FFFFFF[^>]*>\s*</rect>", "", txt, count=1)


def combine_svgs(svgs, gap=30.0):
    """Place several RDKit SVGs side by side in one SVG (plain <g> groups)."""
    parts, x, height = [], 0.0, 0.0
    for svg in svgs:
        head = re.search(r"<svg[^>]*>", svg, re.S).group(0)
        w = float(re.search(r'width=.([0-9.]+)', head).group(1))
        h = float(re.search(r'height=.([0-9.]+)', head).group(1))
        body = svg[svg.index(">", svg.index("<svg")) + 1: svg.rindex("</svg>")]
        parts.append('<g transform="translate(%.1f,0)">%s</g>' % (x, body))
        x += w + gap
        height = max(height, h)
    width = max(x - gap, 1.0)
    return ('<?xml version="1.0" encoding="UTF-8"?>\n'
            '<svg xmlns="http://www.w3.org/2000/svg" width="%.1fpx" height="%.1fpx" '
            'viewBox="0 0 %.1f %.1f">\n%s\n</svg>\n'
            % (width, height, width, height, "\n".join(parts)))


# ---------------------------------------------------------------- widgets
class Tile(QFrame):
    """One molecule preview. Click copies its SVG to the clipboard."""

    def __init__(self, name, svg=None, error=None, zoom=1.0, on_select=None,
                 on_copy=None, badge=None):
        super().__init__()
        self.name, self.svg = name, svg
        self.on_select, self.on_copy = on_select, on_copy
        self.setFrameShape(QFrame.StyledPanel)
        self.setStyleSheet("Tile { background: white; border: 1px solid #ccc; }"
                           "QLabel { color: #222; background: transparent; }")
        lay = QVBoxLayout(self)
        lay.setContentsMargins(8, 8, 8, 6)

        if svg is not None:
            w = QSvgWidget()
            w.load(QByteArray(svg.encode()))
            size = w.renderer().defaultSize()
            w.setFixedSize(QSize(int(size.width() * zoom), int(size.height() * zoom)))
            lay.addWidget(w, alignment=Qt.AlignCenter)
            self.setCursor(Qt.PointingHandCursor)
            self.setToolTip("Click to select, double-click to copy SVG")
        else:
            err = QLabel(error or "error")
            err.setStyleSheet("color: #b00020;")
            err.setWordWrap(True)
            lay.addWidget(err, alignment=Qt.AlignCenter)

        lbl = QLabel(name)
        lbl.setAlignment(Qt.AlignCenter)
        lay.addWidget(lbl)
        if badge:
            text, warn, tip = badge
            b = QLabel(text)
            b.setAlignment(Qt.AlignCenter)
            b.setToolTip(tip)
            b.setStyleSheet("color: %s; font-size: 11px;" % ("#b26a00" if warn else "#777"))
            lay.addWidget(b)

    def set_selected(self, on):
        border = "#2a7de1" if on else "#ccc"
        width = 2 if on else 1
        self.setStyleSheet("Tile { background: white; border: %dpx solid %s; }"
                           "QLabel { color: #222; background: transparent; }" % (width, border))

    def mousePressEvent(self, ev):
        if self.svg is not None and ev.button() == Qt.LeftButton and self.on_select:
            self.on_select(self)

    def mouseDoubleClickEvent(self, ev):
        if self.svg is not None and ev.button() == Qt.LeftButton and self.on_copy:
            self.on_copy(self)


class MainWindow(QMainWindow):
    def __init__(self, smi_path=None):
        super().__init__()
        self.setWindowTitle("Scheme Molecules")
        self.resize(1150, 680)
        self.results = []  # (name, svg or None, error or None)
        self.tiles = []
        self.selected = None

        # ---- left panel: input + options
        self.editor = QPlainTextEdit()
        self.editor.setPlaceholderText("SMILES  name  (optional rotation in degrees), one per line")
        font = self.editor.font()
        font.setFamily("monospace")
        self.editor.setFont(font)

        self.core = QLineEdit("auto")
        self.core.setPlaceholderText("auto (default), a SMILES/SMARTS, or none")
        self.core.setToolTip(
            "auto (or empty): align on the largest skeleton shared by all molecules;\n"
            "  one molecule or nothing shared: straighten its longest chain\n"
            "typed SMILES/SMARTS (e.g. OCCCCCCO, c1ccoc1): always used instead\n"
            "none: no alignment")
        self.hetero = QComboBox()
        self.hetero.addItems(["up", "down", "free"])
        self.zoom = QDoubleSpinBox()
        self.zoom.setRange(0.5, 4.0)
        self.zoom.setSingleStep(0.25)
        self.zoom.setValue(1.5)
        self.zoom.setSuffix(" x")
        self.cols = QSpinBox()
        self.cols.setRange(1, 10)
        self.cols.setValue(3)

        self.level = QCheckBox("Level chain ends (small tilt)")
        self.level.setToolTip("An even-length zigzag has its two ends half a step apart;\n"
                              "this tilts the drawing slightly so they sit at the same height")
        self.flip_v = QCheckBox("Flip vertical")
        self.flip_h = QCheckBox("Flip horizontal")
        self.cip = QCheckBox("Show CIP labels (R/S, E/Z)")
        self.lookup_edit = QLineEdit()
        self.lookup_edit.setPlaceholderText("PubChem name or CAS; separate several with ;")
        self.btn_lookup = QPushButton("Add from PubChem")
        self.btn_lookup.clicked.connect(self.lookup)
        self.lookup_edit.returnPressed.connect(self.lookup)
        lrow = QHBoxLayout()
        lrow.addWidget(self.lookup_edit, 1)
        lrow.addWidget(self.btn_lookup)
        self._thread = None

        form = QFormLayout()
        form.addRow("Align on core", self.core)
        form.addRow("Heteroatom", self.hetero)
        form.addRow("Preview zoom", self.zoom)
        form.addRow("Columns", self.cols)
        form.addRow(self.level)
        form.addRow(self.flip_v)
        form.addRow(self.flip_h)
        form.addRow(self.cip)

        btn_open = QPushButton("Open .smi")
        btn_save = QPushButton("Save .smi")
        btn_export = QPushButton("Export all SVGs...")
        self.live = QCheckBox("Live update while typing")
        btn_run = QPushButton("Run  (Ctrl+Enter)")
        btn_run.setDefault(True)
        btn_copy = QPushButton("Copy SVG (selected / only one)")
        btn_copy_all = QPushButton("Copy all as one SVG")
        btn_run.clicked.connect(self.render)
        btn_copy.clicked.connect(self.copy_selected)
        btn_copy_all.clicked.connect(self.copy_all)
        QShortcut(QKeySequence("Ctrl+Return"), self, self.render)
        QShortcut(QKeySequence("Ctrl+Enter"), self, self.render)
        btn_open.clicked.connect(self.open_smi)
        btn_save.clicked.connect(self.save_smi)
        btn_export.clicked.connect(self.export_all)

        row = QHBoxLayout()
        row.addWidget(btn_open)
        row.addWidget(btn_save)

        left = QWidget()
        lv = QVBoxLayout(left)
        lv.addWidget(QLabel("Molecules"))
        lv.addWidget(self.editor, 1)
        lv.addLayout(lrow)
        lv.addLayout(form)
        lv.addWidget(self.live)
        lv.addWidget(btn_run)
        lv.addWidget(btn_copy)
        lv.addWidget(btn_copy_all)
        lv.addLayout(row)
        lv.addWidget(btn_export)

        # ---- right panel: preview grid
        self.grid_host = QWidget()
        self.grid = QGridLayout(self.grid_host)
        self.grid.setAlignment(Qt.AlignTop | Qt.AlignLeft)
        self.grid.setSpacing(12)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(self.grid_host)

        split = QSplitter()
        split.addWidget(left)
        split.addWidget(scroll)
        split.setSizes([360, 790])
        self.setCentralWidget(split)

        # ---- live update with a short debounce
        self.timer = QTimer(self, singleShot=True, interval=350)
        self.timer.timeout.connect(self.render)
        self.editor.textChanged.connect(self._maybe_live)
        self.core.textChanged.connect(self._maybe_live)
        self.hetero.currentIndexChanged.connect(self.timer.start)
        self.zoom.valueChanged.connect(self.timer.start)
        self.cols.valueChanged.connect(self.timer.start)
        for cb in (self.level, self.flip_v, self.flip_h, self.cip):
            cb.toggled.connect(self.timer.start)

        if smi_path:
            self.load_file(smi_path)
        else:
            self.editor.setPlainText(EXAMPLE)
        self.render()

    # ---- rendering
    def render(self):
        while self.grid.count():
            w = self.grid.takeAt(0).widget()
            if w:
                w.setParent(None)  # hide immediately, not at next event loop
                w.deleteLater()

        self.results, self.tiles, self.selected = [], [], None
        entries = parse(self.editor.toPlainText())
        hetero = self.hetero.currentText()
        txt = self.core.text().strip()
        try:
            if txt.lower() in ("", "auto"):  # typed core always wins over this
                template = auto_template([e[1] for e in entries], hetero) or "chain"
            elif txt.lower() == "none":
                template = None
            else:
                template = make_template(txt)
            core_ok = True
        except ValueError as e:
            template, core_ok = None, False
            self.statusBar().showMessage(str(e))

        for name, smi, angle in entries:
            try:
                m = build(smi, template, hetero)
                if self.level.isChecked():
                    level_ends(m)
                flip_coords(m, self.flip_h.isChecked(), self.flip_v.isChecked())
                rotate_coords(m, angle)  # last, so (90) is the final on-screen angle
                badge = stereo_badge(m)
                if m.HasProp("unaligned"):
                    note = "not aligned: the shared layout would draw its E/Z geometry wrongly"
                    badge = (("%s | not aligned" % badge[0]) if badge else "not aligned (E/Z)",
                             True, (badge[2] + "\n" if badge else "") + note)
                if self.cip.isChecked():
                    annotate_cip(m)
                self.results.append((name, svg_text(m), None, badge))
            except Exception as e:  # noqa: BLE001
                self.results.append((name, None, "%s\n%s" % (e, smi), None))

        ncol = self.cols.value()
        for i, (name, svg, err, badge) in enumerate(self.results):
            tile = Tile(name, svg, err, self.zoom.value(), self.select_tile,
                        self.copy_svg, badge)
            self.tiles.append(tile)
            self.grid.addWidget(tile, i // ncol, i % ncol)

        if core_ok:
            nbad = sum(1 for r in self.results if r[1] is None)
            msg = "%d molecules" % len(self.results)
            if nbad:
                msg += ", %d with errors" % nbad
            self.statusBar().showMessage(msg)

    # ---- actions
    def lookup(self):
        names = [n.strip() for n in self.lookup_edit.text().split(";") if n.strip()]
        if not names or (self._thread is not None and self._thread.isRunning()):
            return
        self.btn_lookup.setEnabled(False)
        self.statusBar().showMessage("Looking up %s on PubChem ..." % "; ".join(names))
        self._thread = LookupThread(names)
        self._thread.got.connect(self._lookup_done)
        self._thread.start()

    def _lookup_done(self, ok, bad):
        self.btn_lookup.setEnabled(True)
        if ok:
            text = self.editor.toPlainText()
            if text and not text.endswith("\n"):
                text += "\n"
            text += "".join("%s  %s\n" % (smi, name) for name, smi in ok)
            self.editor.setPlainText(text)
            self.lookup_edit.clear()
            self.render()
        if bad:
            self.statusBar().showMessage("PubChem lookup failed - " + " | ".join(bad), 12000)
        elif ok:
            self.statusBar().showMessage("Added %d from PubChem" % len(ok), 4000)

    def _maybe_live(self):
        if self.live.isChecked():
            self.timer.start()

    def select_tile(self, tile):
        if self.selected is not None:
            self.selected.set_selected(False)
        self.selected = tile
        tile.set_selected(True)

    def _to_clipboard(self, svg, label):
        mime = QMimeData()
        mime.setData("image/svg+xml", QByteArray(svg.encode()))
        mime.setText(svg)
        QApplication.clipboard().setMimeData(mime)
        self.statusBar().showMessage("Copied %s as SVG to clipboard" % label, 4000)

    def copy_svg(self, tile):
        self.select_tile(tile)
        self._to_clipboard(tile.svg, tile.name)

    def copy_selected(self):
        good = [t for t in self.tiles if t.svg is not None]
        if self.selected is not None:
            self.copy_svg(self.selected)
        elif len(good) == 1:
            self.copy_svg(good[0])
        else:
            self.statusBar().showMessage("Click a tile to select it first, or use Copy all", 4000)

    def copy_all(self):
        good = [t.svg for t in self.tiles if t.svg is not None]
        if good:
            self._to_clipboard(combine_svgs(good), "%d molecules" % len(good))

    def load_file(self, path):
        with open(path) as fh:
            self.editor.setPlainText(fh.read())
        self.setWindowTitle("Scheme Molecules - " + os.path.basename(path))

    def open_smi(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Open SMILES file", "", "SMILES (*.smi *.txt);;All files (*)")
        if path:
            self.load_file(path)

    def save_smi(self):
        path, _ = QFileDialog.getSaveFileName(
            self, "Save SMILES file", "molecules.smi", "SMILES (*.smi)")
        if path:
            with open(path, "w") as fh:
                fh.write(self.editor.toPlainText())
            self.statusBar().showMessage("Saved " + path, 3000)

    def export_all(self):
        good = [(r[0], r[1]) for r in self.results if r[1] is not None]
        if not good:
            QMessageBox.information(self, "Export", "Nothing to export.")
            return
        outdir = QFileDialog.getExistingDirectory(self, "Export SVGs to folder")
        if not outdir:
            return
        for name, svg in good:
            safe = "".join(c if c.isalnum() or c in "-_." else "_" for c in name)
            with open(os.path.join(outdir, safe + ".svg"), "w") as fh:
                fh.write(svg)
        self.statusBar().showMessage(
            "Exported %d SVGs to %s" % (len(good), outdir), 5000)


def selftest():
    """Headless smoke test for a frozen build: python/exe --selftest"""
    app = QApplication([])
    win = MainWindow()
    win.editor.setPlainText("OC[C@@H](O)[C@@H](O)[C@H](O)[C@@H](O)CO sorbitol\nCC(=O)O acetic_acid\n")
    win.render()
    ok = len(win.tiles) == 2 and all(r[1] and "<svg" in r[1] for r in win.results)
    ok = ok and win.results[0][3] is not None and pubchem_ok_parse()
    print("selftest:", "OK" if ok else "FAILED")
    return 0 if ok else 1


def pubchem_ok_parse():
    """Offline check that the lookup code path imports and parses."""
    import io
    class R(io.BytesIO):
        def __enter__(self): return self
        def __exit__(self, *a): return False
    real = urllib.request.urlopen
    urllib.request.urlopen = lambda req, timeout=0: R(
        b'{"PropertyTable":{"Properties":[{"CID":1,"SMILES":"CCO"}]}}')
    try:
        return pubchem_smiles("ethanol") == "CCO"
    finally:
        urllib.request.urlopen = real


def main():
    if "--selftest" in sys.argv:
        sys.exit(selftest())
    app = QApplication(sys.argv)
    win = MainWindow(sys.argv[1] if len(sys.argv) > 1 else None)
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
