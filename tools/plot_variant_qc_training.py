#!/usr/bin/env python3
"""Plot H1 fit rounds from a variant-qc training trace (NumPy, SciPy, Matplotlib)."""
import argparse
import gzip
from pathlib import Path
import sys
import zlib

import numpy as np
from scipy.stats import beta as beta_distribution
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

SCHEMA = "mitoquest.variant_qc.training.v1"
COLUMNS = (
    "iteration", "stage", "alpha", "beta", "training_count", "evaluable_alt_count",
    "selected_alt_count", "changed_alt_count", "call_change_fraction",
    "max_pp_change", "max_log_parameter_change",
)
INTEGER_META = {"max_iter", "dp_threshold", "hq_threshold", "max_alt_alleles", "bins"}
FLOAT_META = {"convergence_eps", "threshold", "pi"}
STOP_REASONS = {
    "converged", "maximum fitting rounds reached",
    "fewer than 3 selected ALT observations", "no evaluable ALT observations",
}
DEFAULT_TITLE = "variant-qc H1 training"

PLOT_THEME = {
    "pdf.fonttype": 42, # 确保 pdf 字体正常矢量
    "font.family": "sans-serif",
    "font.size": 10,
    "figure.figsize": (12.8, 9.2),
    "figure.facecolor": "white",
    "figure.titlesize": 17,
    "figure.titleweight": "bold",
    "axes.facecolor": "white",
    "axes.titlesize": 12,
    "axes.labelsize": 10,
    "axes.edgecolor": "#444444",
    "axes.labelcolor": "#333333",
    "axes.linewidth": 0.8,
    "axes.grid": True,
    "axes.axisbelow": True,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "grid.color": "#b0b0b0",
    "grid.alpha": 0.18,
    "grid.linestyle": "-",
    "text.color": "#333333",
    "xtick.color": "#444444",
    "ytick.color": "#444444",
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "lines.linewidth": 1.5,
    "lines.markersize": 4,
    "legend.fontsize": 8,
    "legend.frameon": True,
    "legend.framealpha": 0.9,
}

PLOT_COLORS = {
    "primary": "#2166ac",
    "secondary": "#c65d28",
    "tertiary": "#26816a",
    "reference": "#666666",
    "muted": "#555555",
    "status_text": "#444444",
    "status_background": "#f2f4f6",
    "warning_text": "#a32323",
    "warning_background": "#ffeded",
    "fallback_text": "#a34d13",
    "density_cmap": "viridis",
}


def _number(text, label, integer=False, positive=False):
    try:
        value = float(text)
        if integer and not (text.isascii() and text.isdecimal()):
            raise ValueError
        if not np.isfinite(value) or value < 0 or (positive and value == 0):
            raise ValueError
        return int(text) if integer else value
    except (ValueError, OverflowError):
        kind = "positive" if positive else "nonnegative"
        raise ValueError(f"{label} must be a finite {kind} "
                         f"{'integer' if integer else 'number'}; got {text!r}") from None


def _validate_metadata(meta):
    if meta.get("schema") != SCHEMA:
        raise ValueError(f"expected #schema={SCHEMA}")
    for key in sorted(INTEGER_META | FLOAT_META):
        if key not in meta:
            raise ValueError(f"missing metadata #{key}")
        meta[key] = _number(meta[key], f"#{key}", integer=key in INTEGER_META,
                            positive=key not in {"dp_threshold", "hq_threshold", "threshold"})
    if meta["convergence_eps"] > 1 or meta["threshold"] > 1 or meta["pi"] >= 1:
        raise ValueError("convergence_eps must be in (0, 1], threshold in [0, 1], pi in (0, 1)")


def _parse_row(fields, rows, meta, line_number):
    label = f"line {line_number}"
    if len(fields) != len(COLUMNS):
        raise ValueError(f"{label}: expected {len(COLUMNS)} tab-separated columns")
    row = dict(zip(COLUMNS, fields))
    for key in COLUMNS:
        if key == "stage":
            continue
        if key in COLUMNS[7:] and not rows:
            if row[key] != ".":
                raise ValueError(f"{label}: initial {key} must be '.' (not a measured change)")
            row[key] = np.nan
        else:
            row[key] = _number(row[key], f"{label}: {key}",
                               integer=key == "iteration" or key.endswith("_count"),
                               positive=key in {"iteration", "alpha", "beta"})
    if row["iteration"] != len(rows) + 1 or row["iteration"] > meta["max_iter"]:
        raise ValueError(f"{label}: iterations must start at 1, be consecutive and not exceed max_iter")
    if row["stage"] != ("refit" if rows else "initial"):
        raise ValueError(f"{label}: stage must be initial for round 1, then refit")
    evaluable = row["evaluable_alt_count"]
    if row["selected_alt_count"] > evaluable:
        raise ValueError(f"{label}: selected_alt_count exceeds evaluable_alt_count")
    if row["training_count"] == 0 and (row["alpha"], row["beta"]) != (0.5, 0.2):
        raise ValueError(f"{label}: zero training observations require fallback alpha=0.5, beta=0.2")
    if rows:
        if evaluable != rows[0]["evaluable_alt_count"]:
            raise ValueError(f"{label}: evaluable_alt_count must remain fixed")
        if row["training_count"] != rows[-1]["selected_alt_count"] or row["training_count"] < 3:
            raise ValueError(f"{label}: refit training_count must equal the previous selected count (at least 3)")
        if row["changed_alt_count"] > evaluable:
            raise ValueError(f"{label}: changed_alt_count exceeds evaluable_alt_count")
        previous = rows[-1]["selected_alt_count"]
        selected, changed = row["selected_alt_count"], row["changed_alt_count"]
        # Changed calls are the symmetric difference of the two selected sets.
        upper = min(previous + selected, 2 * evaluable - previous - selected)
        if not abs(previous - selected) <= changed <= upper or (previous + selected - changed) % 2:
            raise ValueError(f"{label}: changed_alt_count is inconsistent with selected counts")
        if row["call_change_fraction"] > 1 or row["max_pp_change"] > 1:
            raise ValueError(f"{label}: call_change_fraction and max_pp_change must be in [0, 1]")
        expected = row["changed_alt_count"] / evaluable if evaluable else 0.0
        if not np.isclose(row["call_change_fraction"], expected, rtol=1e-12, atol=1e-15):
            raise ValueError(f"{label}: call_change_fraction disagrees with changed/evaluable counts")
    return row


def load_training(tsv_path):
    """Return validated metadata and row dictionaries; initial changes remain NaN."""
    try:
        with open(tsv_path, "rb") as probe:
            compressed = probe.read(2) == b"\x1f\x8b"
        opener = gzip.open if compressed else open
        # gzip also reads concatenated BGZF members, irrespective of the filename.
        with opener(tsv_path, "rt", encoding="utf-8") as stream:
            lines = list(stream)
    except (OSError, EOFError, UnicodeError, zlib.error) as error:
        raise ValueError(f"cannot read training TSV {tsv_path}: {error}") from None

    meta, rows, header_seen = {}, [], False
    for line_number, line in enumerate(lines, 1):
        line = line.rstrip("\r\n")
        if not line:
            continue
        if "stop_reason" in meta:
            raise ValueError(f"line {line_number}: content after stop_reason footer")
        if line.startswith("#"):
            key, separator, value = line[1:].partition("=")
            if not separator or key in meta:
                raise ValueError(f"line {line_number}: malformed or duplicate metadata")
            if key == "stop_reason":
                if not rows or value not in STOP_REASONS:
                    raise ValueError(f"line {line_number}: invalid stop_reason footer")
            elif header_seen or key not in INTEGER_META | FLOAT_META | {"schema"}:
                raise ValueError(f"line {line_number}: unexpected metadata #{key}")
            meta[key] = value
        elif not header_seen:
            if tuple(line.split("\t")) != COLUMNS:
                raise ValueError("missing or invalid TSV columns; expected: " + "\t".join(COLUMNS))
            _validate_metadata(meta)
            header_seen = True
        else:
            rows.append(_parse_row(line.split("\t"), rows, meta, line_number))
    if not header_seen or not rows:
        raise ValueError("training TSV must contain the schema, column header and at least one fit row")
    meta.setdefault("stop_reason", None)
    return meta, rows


@matplotlib.rc_context(PLOT_THEME)
def make_plot(meta, rows, title=DEFAULT_TITLE):
    """Build the four-panel figure; the caller owns closing/saving it."""
    data = {key: np.array([row[key] for row in rows]) for key in COLUMNS if key != "stage"}
    rounds = data["iteration"]
    count = len(rows)
    complete = meta["stop_reason"] is not None
    last_label = "final" if complete else "latest"
    font_size = matplotlib.rcParams["font.size"]
    line_width = matplotlib.rcParams["lines.linewidth"]
    fig, axes = plt.subplots(2, 2)
    parameters, changes, counts, density = axes.flat
    try:
        for ax in (parameters, changes, counts):
            ax.set_xlabel("H1 fit round")
            ax.set_xlim(0.7, max(count + 0.3, 1.3))
            if count == 1:
                ax.set_xticks([1])
            else:
                ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=7))

        for key, color in (("alpha", PLOT_COLORS["primary"]), ("beta", PLOT_COLORS["secondary"])):
            values = data[key]
            parameters.plot(rounds, values, "o-", color=color, linewidth=1.2 * line_width,
                            label=f"{key}: initial {values[0]:.4g}; {last_label} {values[-1]:.4g}")
        parameters.set(title="H1 shape parameters", ylabel="Beta shape parameter")
        parameters.legend()

        eps = meta["convergence_eps"]
        if count > 1:
            for key, label, color, marker in (
                ("call_change_fraction", "ALT call-change fraction", PLOT_COLORS["primary"], "o"),
                ("max_pp_change", "Maximum posterior change", PLOT_COLORS["secondary"], "s"),
                ("max_log_parameter_change", "Maximum log-parameter change", PLOT_COLORS["tertiary"], "^"),
            ):
                changes.plot(rounds[1:], data[key][1:], marker=marker, color=color, label=label)
        else:
            changes.text(0.5, 0.48, "No refit\nInitial changes are not defined", ha="center",
                         va="center", transform=changes.transAxes, color=PLOT_COLORS["muted"])
        changes.axhline(eps, color=PLOT_COLORS["reference"], linestyle="--",
                        linewidth=2 * line_width / 3, label=f"eps = {eps:g}")
        changes.set_yscale("symlog", linthresh=eps, linscale=1)
        changes.set_ylim(bottom=0)
        changes.set(title="Between-fit convergence", ylabel="Change (symlog; zero included)")
        changes.legend(loc="best")

        counts.plot(rounds, data["training_count"], "o-", color=PLOT_COLORS["primary"],
                    label="Training observations (input to this fit)")
        counts.plot(rounds, data["selected_alt_count"], "s--", color=PLOT_COLORS["secondary"],
                    label="Selected ALT candidates (this model's PP > threshold)")
        counts.axhline(data["evaluable_alt_count"][0], color=PLOT_COLORS["reference"], linestyle=":",
                       label=f"Fixed evaluable ALT candidates: {data['evaluable_alt_count'][0]}")
        counts.set(title="Training input versus model selections", ylabel="Observation / ALT candidate count")
        counts.set_ylim(0, max(1, counts.get_ylim()[1]))
        counts.yaxis.set_major_locator(MaxNLocator(integer=True))
        counts.legend()

        edge = np.geomspace(1e-5, 0.02, 180)
        x = np.unique(np.concatenate((edge, np.linspace(0.02, 0.98, 650), 1 - edge)))
        cmap = plt.get_cmap(PLOT_COLORS["density_cmap"])
        for index, row in enumerate(rows):
            initial, final = index == 0, index == count - 1
            label = (f"Initial / {last_label} (round 1)" if count == 1 else
                     "Initial (round 1)" if initial else
                     f"{last_label.capitalize()} (round {count})" if final else None)
            if label and row["training_count"] == 0:
                label += " — fallback"
            with np.errstate(over="ignore", invalid="ignore", under="ignore"):
                y = np.exp(beta_distribution.logpdf(x, row["alpha"], row["beta"]))
            if not np.all(np.isfinite(y)):
                raise ValueError(f"round {row['iteration']}: Beta PDF exceeds numerical plotting limits")
            density.plot(x, y, color=cmap(index / max(count - 1, 1)),
                         linestyle="--" if initial else "-", linewidth=line_width * (1.6 if final else 1),
                         alpha=1 if initial or final else 0.35, label=label)
        density.set_yscale("symlog", linthresh=1)
        density.set(xlim=(0, 1), ylim=(0, None), title="H1 Beta distributions (not observed VAFs)",
                    xlabel="Latent allele frequency", ylabel="Latent AF density (Beta PDF; symlog)")
        density.legend()
        if count > 2:
            scalar = matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(1, count), cmap=cmap)
            fig.colorbar(scalar, ax=density, pad=0.025, fraction=0.045,
                         ticks=sorted({1, (count + 1) // 2, count}), label="H1 fit round")

        fig.suptitle(title, y=0.987)
        status = (f"Stop: {meta['stop_reason']}" if complete else
                  "INCOMPLETE TRACE — missing stop_reason footer; convergence is unknown")
        fig.text(0.5, 0.946, status, ha="center", fontsize=1.05 * font_size,
                 color=PLOT_COLORS["status_text" if complete else "warning_text"],
                 weight="normal" if complete else "bold",
                 bbox=dict(facecolor=PLOT_COLORS["status_background" if complete else "warning_background"],
                           edgecolor="none", pad=5))
        fig.text(0.5, 0.916, f"max_iter = {meta['max_iter']}   |   eps = {eps:g}   |   "
                 f"posterior threshold = {meta['threshold']:g}", ha="center")
        fig.text(0.5, 0.045, "Rows represent completed H1 fits, not optimizer internal steps. "
                 "Beta PDFs describe latent AF, not a Beta-Binomial PMF or an empirical histogram.",
                 ha="center", fontsize=0.8 * font_size, color=PLOT_COLORS["muted"])
        if np.any(data["training_count"] == 0):
            fig.text(0.5, 0.022, "Zero training observations: fallback H1 parameters "
                     "(alpha = 0.5, beta = 0.2); this is not a data-fitted distribution.",
                     ha="center", fontsize=0.9 * font_size, color=PLOT_COLORS["fallback_text"])
        fig.tight_layout(rect=(0.015, 0.075, 0.99, 0.89), h_pad=2.5, w_pad=2)
        return fig
    except Exception:
        plt.close(fig)
        raise


def _positive_dpi(text):
    try:
        return _number(text, "dpi", integer=True, positive=True)
    except ValueError as error:
        raise argparse.ArgumentTypeError(str(error)) from None


def _validate_paths(input_path, output_path):
    if output_path.suffix.lower() not in {".png", ".pdf", ".svg"}:
        raise ValueError("output format must be PNG, PDF or SVG (.png/.pdf/.svg)")
    if not input_path.is_file():
        raise ValueError(f"input is not a readable file: {input_path}")
    if input_path.resolve() == output_path.resolve() or (
            output_path.exists() and input_path.samefile(output_path)):
        raise ValueError("output must not overwrite input (including symlinks or hardlinks)")
    if not output_path.parent.is_dir():
        raise ValueError(f"output directory does not exist: {output_path.parent}")
    if output_path.exists() and not output_path.is_file():
        raise ValueError(f"output is not a file: {output_path}")


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Plot completed H1 fits from a variant-qc training TSV (plain, gzip or BGZF).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--input", required=True, help="Training trace TSV path.")
    parser.add_argument("-o", "--output", required=True, help="Output figure (.png/.pdf/.svg).")
    parser.add_argument("--title", default=DEFAULT_TITLE, help="Figure title.")
    parser.add_argument("--dpi", type=_positive_dpi, default=180, help="Positive raster resolution.")
    args = parser.parse_args(argv)
    try:
        input_path, output_path = Path(args.input), Path(args.output)
        _validate_paths(input_path, output_path)
        meta, rows = load_training(input_path)
        if meta["stop_reason"] is None:
            print("[plot_variant_qc_training] WARNING: INCOMPLETE TRACE: missing stop_reason footer; "
                  "plotting available fits, not claiming convergence.", file=sys.stderr)
        fig = make_plot(meta, rows, args.title)
        try:
            fig.savefig(output_path, dpi=args.dpi)
        finally:
            plt.close(fig)
    except (OSError, ValueError, OverflowError, RuntimeError) as error:
        print(f"[plot_variant_qc_training] error: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
