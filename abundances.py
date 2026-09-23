"""
abundances.py
=============
Convert [X/Fe] abundance columns to [X/H] for the X/H re-run of the
phylogenetic benchmark, using [X/H] = [X/Fe] + [Fe/H].

This lives in the data-prep layer (alongside the notebook cleaning), NOT in
phylograph.py: the tree/graph/spectral core is representation-agnostic and
must stay identical between the X/Fe and X/H experiments. The ONLY thing that
changes is the abundance matrix handed to RunConfig.

Two experiments, selected by `include_feh`:

  include_feh=False : re-reference the SAME 7 elements to H.
      [Mg/Fe]..[Ni/Fe]  ->  [Mg/H]..[Ni/H]     (7 features)
      Tests whether changing the reference element alone changes the result.
      The [Fe/H] information is smeared into all 7 features but is not present
      as its own axis.

  include_feh=True  : re-reference AND add [Fe/H] as its own feature.
      [Mg/H]..[Ni/H] + [Fe/H]                    (8 features)
      This is the experiment that re-introduces the overall-metallicity axis
      that [X/Fe] removes. It is the one whose separation you predicted to be
      >= the [X/Fe] separation. Use this as the primary X/H run.

IMPORTANT (dropna trap): the conversion is done on the full cleaned frame with
[Fe/H] present, and the returned frame contains only the new feature columns.
Because RunConfig.frames() calls df[elements].dropna(), any row missing [Fe/H]
would be dropped from the X/H run but kept in the X/Fe run, changing the star
set per population and reintroducing the leaf-count asymmetry that caused the
earlier bug. To keep the two experiments on the SAME stars, this function drops
rows lacking [Fe/H] (or any used element) HERE, up front, identically for every
population, and reports how many rows each population loses.
"""

from __future__ import annotations
import pandas as pd


# The 7 elements, in the two column-naming conventions actually used in the
# notebook (APOGEE upper-underscore vs Nissen slash notation).
_APOGEE_XFE = ["MG_FE", "AL_FE", "SI_FE", "CA_FE", "TI_FE", "CR_FE", "NI_FE"]
_NISSEN_XFE = ["Mg/Fe", "Al/Fe", "Si/Fe", "Ca/Fe", "Ti/Fe", "Cr/Fe", "Ni/Fe"]


def _xh_name(xfe_col: str) -> str:
    """MG_FE -> MG_H ; 'Mg/Fe' -> 'Mg/H'. Purely a rename of the reference."""
    if xfe_col.endswith("_FE"):
        return xfe_col[:-3] + "_H"
    if xfe_col.endswith("/Fe"):
        return xfe_col[:-3] + "/H"
    raise ValueError(f"unrecognised [X/Fe] column name: {xfe_col!r}")


def to_xh(df: pd.DataFrame,
          xfe_cols: list[str],
          feh_col: str,
          *,
          include_feh: bool = True,
          label: str = "") -> tuple[pd.DataFrame, list[str]]:
    """Return (frame_with_XH_columns, xh_element_list).

    df       : a cleaned frame that still contains xfe_cols AND feh_col.
    xfe_cols : the 7 [X/Fe] column names for this frame's naming scheme.
    feh_col  : the [Fe/H] column name ('FE_H' for APOGEE, 'Fe/H' for Nissen).
    include_feh : if True, [Fe/H] is appended as its own feature (8 features);
                  if False, only the 7 re-referenced [X/H] columns are returned.
    label    : optional name for the printed row-loss report.

    The returned element list is what you pass to RunConfig as `elements`.
    Only the returned feature columns survive, and rows missing ANY used
    column (including feh_col) are dropped here so the star set is fixed
    before RunConfig.frames() ever calls dropna.
    """
    missing = [c for c in (*xfe_cols, feh_col) if c not in df.columns]
    if missing:
        raise KeyError(f"[{label}] frame is missing columns: {missing}")

    # Drop rows lacking any input we need, UP FRONT and identically for every
    # population, so X/Fe-run and X/H-run see the same stars.
    needed = [*xfe_cols, feh_col]
    n_before = len(df)
    df = df.dropna(subset=needed).reset_index(drop=True)
    n_after = len(df)
    if n_after < n_before:
        print(f"[{label}] dropped {n_before - n_after} rows lacking "
              f"[Fe/H] or an element ({n_after} remain)")

    feh = df[feh_col]
    out = pd.DataFrame(index=df.index)
    xh_cols = []
    for xfe in xfe_cols:
        xh = _xh_name(xfe)
        out[xh] = df[xfe] + feh          # [X/H] = [X/Fe] + [Fe/H]
        xh_cols.append(xh)

    if include_feh:
        feh_feature = _xh_name(feh_col)  # FE_H -> FE_H (unchanged) / 'Fe/H' -> 'Fe/H'
        # feh_col already IS [Fe/H]; keep its original name as the feature.
        out[feh_col] = feh
        xh_cols.append(feh_col)

    return out, xh_cols


# Convenience wrappers that know the two naming schemes, so the notebook call
# site stays a one-liner per frame.

def apogee_to_xh(df, *, include_feh=True, label=""):
    return to_xh(df, _APOGEE_XFE, "FE_H", include_feh=include_feh, label=label)


def nissen_to_xh(df, *, include_feh=True, label=""):
    return to_xh(df, _NISSEN_XFE, "Fe/H", include_feh=include_feh, label=label)
