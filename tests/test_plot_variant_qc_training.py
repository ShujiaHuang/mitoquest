"""Self-contained training-trace parser, CLI and headless rendering tests."""
import contextlib
import gzip
import importlib.util
import io
import os
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import zlib

import numpy as np

SCRIPT = Path(__file__).resolve().parents[1] / "tools" / "plot_variant_qc_training.py"
SPEC = importlib.util.spec_from_file_location("plot_variant_qc_training", SCRIPT)
plotter = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(plotter)

HEADER = ("iteration\tstage\talpha\tbeta\ttraining_count\tevaluable_alt_count\t"
          "selected_alt_count\tchanged_alt_count\tcall_change_fraction\t"
          "max_pp_change\tmax_log_parameter_change")
INITIAL = "1\tinitial\t0.5\t0.20000000000000001\t12\t10\t6\t.\t.\t.\t."
REFIT = "2\trefit\t0.75\t0.25\t6\t10\t5\t1\t0.10000000000000001\t0.12\t0.40546510810816438"
FINAL = "3\trefit\t0.75\t0.25\t5\t10\t5\t0\t0\t0\t0"
ZERO = "1\tinitial\t0.5\t0.2\t0\t0\t0\t.\t.\t.\t."


def trace(rows=None, stop="converged", max_iter=20):
    metadata = (f"#schema=mitoquest.variant_qc.training.v1\n#max_iter={max_iter}\n"
                "#convergence_eps=0.001\n#threshold=0.9\n#pi=0.00082845\n"
                "#dp_threshold=100\n#hq_threshold=20\n#max_alt_alleles=2\n#bins=100\n")
    result = metadata + HEADER + "\n" + "\n".join(
        [INITIAL, REFIT, FINAL] if rows is None else rows) + "\n"
    return result if stop is None else result + f"#stop_reason={stop}\n"


def replace_field(row, index, value):
    fields = row.split("\t")
    fields[index] = value
    return "\t".join(fields)


def bgzf_member(data):
    compressor = zlib.compressobj(wbits=-15)
    payload = compressor.compress(data) + compressor.flush()
    header = bytes.fromhex("1f8b08040000000000ff060042430200")
    return (header + struct.pack("<H", len(payload) + 25) + payload
            + struct.pack("<II", zlib.crc32(data), len(data)))


class TrainingPlotTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.root = Path(self.directory.name)
        self.input = self.root / "training.tsv"

    def tearDown(self):
        plotter.plt.close("all")
        self.directory.cleanup()

    def load(self, content):
        self.input.write_text(content, encoding="utf-8")
        return plotter.load_training(self.input)

    def invoke(self, args):
        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr):
            try:
                code = plotter.main(args)
            except SystemExit as error:
                code = error.code
        return code, stderr.getvalue()

    def test_valid_multiple_rounds_preserve_counts_precision_and_missing_changes(self):
        meta, rows = self.load(trace())
        self.assertEqual(meta["stop_reason"], "converged")
        self.assertEqual(meta["max_iter"], 20)
        self.assertEqual([row["iteration"] for row in rows], [1, 2, 3])
        self.assertEqual(rows[0]["training_count"], 12)
        self.assertEqual(rows[0]["selected_alt_count"], 6)
        self.assertEqual(rows[1]["training_count"], 6)
        self.assertEqual(rows[1]["selected_alt_count"], 5)
        self.assertEqual(rows[1]["max_log_parameter_change"], float("0.40546510810816438"))
        self.assertTrue(all(np.isnan(rows[0][key]) for key in plotter.COLUMNS[7:]))
        self.assertTrue(all(rows[2][key] == 0 for key in plotter.COLUMNS[7:]))

    def test_single_round_and_all_stop_reasons(self):
        for reason, row, limit in (
            ("maximum fitting rounds reached", INITIAL, 1),
            ("fewer than 3 selected ALT observations", replace_field(INITIAL, 6, "2"), 20),
            ("no evaluable ALT observations", ZERO, 20),
        ):
            with self.subTest(reason=reason):
                meta, rows = self.load(trace([row], reason, limit))
                self.assertEqual(meta["stop_reason"], reason)
                self.assertEqual(len(rows), 1)

    def test_zero_training_is_fallback_even_with_evaluable_candidates(self):
        zero_candidates = replace_field(replace_field(ZERO, 5, "5"), 6, "2")
        meta, rows = self.load(trace([zero_candidates], "fewer than 3 selected ALT observations"))
        self.assertEqual((rows[0]["alpha"], rows[0]["beta"], rows[0]["training_count"]), (0.5, 0.2, 0))
        fig = plotter.make_plot(meta, rows)
        self.assertTrue(any("Zero training observations: fallback" in text.get_text() for text in fig.texts))
        self.assertIn("fallback", fig.axes[3].lines[0].get_label())

    def test_magic_detected_gzip_bgzf_and_misleading_suffix(self):
        content = trace().encode()
        for name, data in (
            ("compressed.tsv", gzip.compress(content)),
            ("blocks.unknown", bgzf_member(content[:111]) + bgzf_member(content[111:]) + bgzf_member(b"")),
            ("plain.tsv.gz", content),
        ):
            with self.subTest(name=name):
                path = self.root / name
                path.write_bytes(data)
                meta, rows = plotter.load_training(path)
                self.assertEqual(meta["schema"], plotter.SCHEMA)
                self.assertEqual(len(rows), 3)

    def test_incomplete_trace_warns_without_claiming_final_or_converged(self):
        meta, rows = self.load(trace(stop=None))
        self.assertIsNone(meta["stop_reason"])
        fig = plotter.make_plot(meta, rows)
        status = " ".join(text.get_text() for text in fig.texts)
        self.assertIn("INCOMPLETE TRACE", status)
        self.assertNotIn("Stop: converged", status)
        self.assertIn("Latest", fig.axes[3].lines[-1].get_label())
        code, stderr = self.invoke(["-i", str(self.input), "-o", str(self.root / "partial.png"), "--dpi", "40"])
        self.assertEqual(code, 0, stderr)
        self.assertIn("WARNING: INCOMPLETE TRACE", stderr)
        self.assertGreater((self.root / "partial.png").stat().st_size, 1000)

    def test_invalid_schema_headers_and_footer(self):
        valid = trace()
        for content in (
            "", valid.replace("#schema=mitoquest.variant_qc.training.v1\n", ""),
            valid.replace("training.v1", "training.v2"),
            valid.replace(HEADER, HEADER.replace("alpha", "shape")),
            valid.replace(HEADER, HEADER + "\textra"),
            valid.replace(HEADER, HEADER.replace("\tstage", "")),
            valid.replace("#max_iter=20\n", ""), valid.replace("#bins=100", "#bins=100\n#bins=100"),
            valid.replace("#bins=100", "#unknown=100"), valid.replace("#bins=100", "#bins"),
            valid.replace("#stop_reason=converged", "#stop_reason=success"),
            trace([], stop=None), trace([]), valid + FINAL + "\n",
            valid.replace(FINAL, FINAL + "\textra"), valid.replace(REFIT, "\t".join(REFIT.split("\t")[:-1])),
        ):
            with self.subTest(content=content):
                with self.assertRaises(ValueError):
                    self.load(content)

    def test_invalid_numeric_metadata(self):
        for key, value in (
            ("max_iter", "0"), ("max_iter", "2"), ("max_iter", "1.5"),
            ("convergence_eps", "0"), ("convergence_eps", "1.01"), ("convergence_eps", "nan"),
            ("threshold", "-0.1"), ("threshold", "1.01"), ("pi", "0"), ("pi", "1"),
            ("dp_threshold", "-1"), ("hq_threshold", "inf"), ("bins", "0"), ("max_alt_alleles", "0"),
        ):
            with self.subTest(key=key, value=value):
                lines = trace().splitlines()
                content = "\n".join(f"#{key}={value}" if line.startswith(f"#{key}=") else line for line in lines)
                with self.assertRaises(ValueError):
                    self.load(content)

    def test_valid_metadata_boundaries(self):
        for threshold in ("0", "1"):
            content = trace().replace("#threshold=0.9", f"#threshold={threshold}")
            content = content.replace("#convergence_eps=0.001", "#convergence_eps=1")
            content = content.replace("#dp_threshold=100", "#dp_threshold=0").replace("#hq_threshold=20", "#hq_threshold=0")
            self.load(content)

    def test_invalid_rows(self):
        cases = [(0, 0, "0"), (0, 0, "2"), (1, 0, "3"), (2, 0, "2"),
                 (0, 1, "refit"), (1, 1, "initial"), (0, 4, "-1"), (0, 4, "1.5"),
                 (0, 5, "-1"), (0, 6, "11"), (1, 4, "5"), (1, 5, "11"),
                 (1, 7, "11"), (1, 7, "-1"), (1, 7, "1.1"), (1, 8, "0.2"),
                 (1, 8, "1.01"), (1, 9, "1.01"), (1, 10, "-0.1")]
        cases += [(0, index, value) for index in (2, 3) for value in ("0", "-1", "nan", "inf", "1e999", "oops")]
        cases += [(0, index, "0") for index in range(7, 11)]
        cases += [(1, index, ".") for index in range(7, 11)]
        cases += [(1, index, "nan") for index in range(8, 11)]
        for row_index, field_index, value in cases:
            with self.subTest(row=row_index, field=field_index, value=value):
                rows = [INITIAL, REFIT, FINAL]
                rows[row_index] = replace_field(rows[row_index], field_index, value)
                with self.assertRaises(ValueError):
                    self.load(trace(rows))

    def test_inconsistent_selected_sets_and_fallback_parameters(self):
        for selected, changed in ((4, 0), (5, 1), (10, 7)):
            row = replace_field(replace_field(FINAL, 6, str(selected)), 7, str(changed))
            row = replace_field(row, 8, str(changed / 10))
            with self.subTest(selected=selected, changed=changed):
                with self.assertRaisesRegex(ValueError, "inconsistent with selected counts"):
                    self.load(trace([INITIAL, REFIT, row]))
        with self.assertRaisesRegex(ValueError, "fallback alpha=0.5, beta=0.2"):
            self.load(trace([replace_field(ZERO, 2, "1")], "no evaluable ALT observations"))

    def test_malformed_compression_and_encoding(self):
        for content in (b"\x1f\x8bgarbage", gzip.compress(trace().encode())[:-8], b"\xff\xfe"):
            with self.subTest(content=content[:20]):
                self.input.write_bytes(content)
                with self.assertRaisesRegex(ValueError, "cannot read training TSV"):
                    plotter.load_training(self.input)

    def test_plot_curves_keep_zero_changes_and_omit_initial_changes(self):
        self.assertEqual(plotter.matplotlib.get_backend().lower(), "agg")
        meta, rows = self.load(trace())
        fig = plotter.make_plot(meta, rows, "Custom title")
        self.assertEqual(fig._suptitle.get_text(), "Custom title")
        changes = fig.axes[1]
        self.assertEqual(changes.get_yscale(), "symlog")
        for curve in changes.lines[:3]:
            np.testing.assert_array_equal(curve.get_xdata(), [2, 3])
            self.assertEqual(curve.get_ydata()[-1], 0)
        np.testing.assert_array_equal(fig.axes[2].lines[0].get_ydata(), [12, 6, 5])
        np.testing.assert_array_equal(fig.axes[2].lines[1].get_ydata(), [6, 5, 5])
        density = fig.axes[3]
        self.assertEqual(len(density.get_legend().get_texts()), 2)
        self.assertEqual(density.lines[0].get_linestyle(), "--")
        self.assertGreater(density.lines[-1].get_linewidth(), density.lines[0].get_linewidth())
        for curve in density.lines:
            self.assertTrue(np.all((curve.get_xdata() > 0) & (curve.get_xdata() < 1)))
            self.assertTrue(np.all(np.isfinite(curve.get_ydata())))
            self.assertTrue(np.all(curve.get_ydata() >= 0))

    def test_custom_theme_and_palette_preserve_curve_data(self):
        meta, rows = self.load(trace())
        baseline = plotter.make_plot(meta, rows)
        theme = {"font.family": "serif", "font.size": 9, "figure.figsize": (14, 10),
                 "figure.facecolor": "#fafafa", "figure.titlesize": 15, "axes.labelsize": 11,
                 "axes.grid": False, "axes.spines.top": True, "lines.linewidth": 2.4,
                 "lines.markersize": 6, "legend.fontsize": 7}
        palette = {"primary": "#663399", "secondary": "#a55194", "tertiary": "#008080",
                   "reference": "#777777", "density_cmap": "plasma"}
        with patch.dict(plotter.PLOT_THEME, theme), patch.dict(plotter.PLOT_COLORS, palette):
            fig = plotter.make_plot(meta, rows)
        np.testing.assert_array_equal(fig.get_size_inches(), [14, 10])
        self.assertEqual(fig.get_facecolor(), plotter.matplotlib.colors.to_rgba("#fafafa"))
        self.assertEqual(fig._suptitle.get_fontsize(), 15)
        self.assertEqual(fig._suptitle.get_fontfamily(), ["serif"])
        parameters, changes, counts, density = fig.axes[:4]
        self.assertEqual(parameters.xaxis.label.get_fontsize(), 11)
        self.assertTrue(parameters.spines["top"].get_visible())
        self.assertFalse(any(line.get_visible() for line in parameters.get_xgridlines()))
        self.assertEqual(parameters.get_legend().get_texts()[0].get_fontsize(), 7)
        self.assertAlmostEqual(parameters.lines[0].get_linewidth(), 1.2 * 2.4)
        self.assertEqual(parameters.lines[0].get_markersize(), 6)
        for ax in (parameters, changes, counts):
            self.assertEqual(ax.lines[0].get_color(), palette["primary"])
            self.assertEqual(ax.lines[1].get_color(), palette["secondary"])
        self.assertEqual(changes.lines[2].get_color(), palette["tertiary"])
        self.assertEqual(changes.lines[3].get_color(), palette["reference"])
        self.assertEqual(changes.lines[0].get_linewidth(), 2.4)
        self.assertEqual(counts.lines[2].get_color(), palette["reference"])
        cmap = plotter.plt.get_cmap("plasma")
        for index, curve in enumerate(density.lines):
            np.testing.assert_array_equal(curve.get_color(), cmap(index / 2))
        self.assertEqual(fig.axes[4].collections[-1].get_cmap().name, "plasma")
        for original, styled in zip(baseline.axes[:4], fig.axes[:4]):
            self.assertEqual(len(original.lines), len(styled.lines))
            for before, after in zip(original.lines, styled.lines):
                np.testing.assert_array_equal(before.get_xdata(), after.get_xdata())
                np.testing.assert_array_equal(before.get_ydata(), after.get_ydata())
        fig.savefig(self.root / "styled.png", dpi=60)
        self.assertEqual(plotter.plt.imread(self.root / "styled.png").shape[:2], (600, 840))

    def test_theme_restores_global_settings_on_success_and_failure(self):
        meta, rows = self.load(trace())
        with plotter.matplotlib.rc_context({"font.size": 6, "axes.grid": False}):
            original = plotter.matplotlib.rcParams.copy()
            for _ in range(2):
                fig = plotter.make_plot(meta, rows)
                self.assertEqual(dict(plotter.matplotlib.rcParams), dict(original))
                plotter.plt.close(fig)
            figures = plotter.plt.get_fignums()
            with patch.dict(plotter.PLOT_COLORS, {"density_cmap": "not-a-colormap"}):
                with self.assertRaises(ValueError):
                    plotter.make_plot(meta, rows)
            self.assertEqual(plotter.plt.get_fignums(), figures)
            self.assertEqual(dict(plotter.matplotlib.rcParams), dict(original))

    def test_status_warning_and_fallback_use_configured_palette(self):
        palette = {"status_text": "#123456", "status_background": "#eeeeee",
                   "warning_text": "#990000", "warning_background": "#ffeedd",
                   "fallback_text": "#884400", "muted": "#555566"}
        with patch.dict(plotter.PLOT_COLORS, palette):
            for stop, prefix in (("no evaluable ALT observations", "status"), (None, "warning")):
                with self.subTest(stop=stop):
                    meta, rows = self.load(trace([ZERO], stop))
                    fig = plotter.make_plot(meta, rows)
                    status = next(text for text in fig.texts if text.get_bbox_patch() is not None)
                    self.assertEqual(status.get_color(), palette[f"{prefix}_text"])
                    self.assertEqual(status.get_bbox_patch().get_facecolor(),
                                     plotter.matplotlib.colors.to_rgba(palette[f"{prefix}_background"]))
                    fallback = next(text for text in fig.texts if "Zero training" in text.get_text())
                    self.assertEqual(fallback.get_color(), palette["fallback_text"])
                    self.assertEqual(fig.axes[1].texts[0].get_color(), palette["muted"])
                    plotter.plt.close(fig)

    def test_zero_candidate_plot_uses_integer_count_ticks(self):
        meta, rows = self.load(trace([ZERO], "no evaluable ALT observations"))
        fig = plotter.make_plot(meta, rows)
        ticks = fig.axes[2].get_yticks()
        np.testing.assert_array_equal(ticks, np.round(ticks))
        self.assertEqual(fig.axes[2].get_ylim(), (0, 1))

    def test_initial_only_no_refit_and_twenty_round_legend(self):
        meta, rows = self.load(trace([INITIAL], "maximum fitting rounds reached", 1))
        fig = plotter.make_plot(meta, rows)
        self.assertTrue(any("No refit" in text.get_text() for text in fig.axes[1].texts))
        self.assertEqual(len(fig.axes[1].lines), 1)
        self.assertIn("Initial / final", fig.axes[3].lines[0].get_label())
        for ax in fig.axes[:3]:
            np.testing.assert_array_equal(ax.get_xticks(), [1])
        more_rows = [INITIAL, REFIT] + [replace_field(FINAL, 0, str(index)) for index in range(3, 21)]
        meta, rows = self.load(trace(more_rows, "maximum fitting rounds reached"))
        fig = plotter.make_plot(meta, rows)
        self.assertEqual(len(fig.axes[3].lines), 20)
        self.assertEqual(len(fig.axes[3].get_legend().get_texts()), 2)
        self.assertEqual(fig.axes[4].get_ylabel(), "H1 fit round")

    def test_cli_bad_paths_output_format_and_dpi(self):
        self.load(trace())
        base = ["-i", str(self.input), "-o", str(self.root / "out.png")]
        cases = [(base + ["--dpi", value], "dpi") for value in ("0", "-2", "nan", "1.5")]
        cases += [([], "required"), (["-i", str(self.input)], "required"),
                  (["-i", str(self.root / "absent"), "-o", str(self.root / "out.gif")], "output format"),
                  (["-i", str(self.root / "absent"), "-o", str(self.root / "out.png")], "input"),
                  (["-i", str(self.root), "-o", str(self.root / "out.png")], "input"),
                  (["-i", str(self.input), "-o", str(self.root / "absent" / "out.png")], "output directory")]
        directory = self.root / "directory.png"
        directory.mkdir()
        cases.append((["-i", str(self.input), "-o", str(directory)], "not a file"))
        for args, message in cases:
            with self.subTest(args=args):
                code, stderr = self.invoke(args)
                self.assertNotEqual(code, 0)
                self.assertIn(message, stderr)
                self.assertNotIn("Traceback", stderr)
        self.assertFalse((self.root / "out.png").exists())

    def test_cli_rejects_same_path_symlink_and_hardlink_without_overwriting(self):
        source = self.root / "source.png"
        source.write_text(trace(), encoding="utf-8")
        symlink, hardlink = self.root / "symlink.png", self.root / "hardlink.png"
        symlink.symlink_to(source)
        os.link(source, hardlink)
        original = source.read_bytes()
        for output in (source, symlink, hardlink):
            code, stderr = self.invoke(["-i", str(source), "-o", str(output)])
            self.assertNotEqual(code, 0)
            self.assertIn("overwrite input", stderr)
            self.assertEqual(source.read_bytes(), original)

    def test_cli_malformed_input_is_friendly_and_does_not_write(self):
        self.input.write_bytes(b"\x1f\x8bgarbage")
        output = self.root / "out.png"
        code, stderr = self.invoke(["-i", str(self.input), "-o", str(output)])
        self.assertNotEqual(code, 0)
        self.assertIn("cannot read training TSV", stderr)
        self.assertNotIn("Traceback", stderr)
        self.assertFalse(output.exists())

    def test_cli_renders_png_pdf_svg_with_agg(self):
        fixtures = {"png": trace(), "pdf": trace([INITIAL], "maximum fitting rounds reached", 1),
                    "svg": trace([ZERO], "no evaluable ALT observations")}
        for extension, content in fixtures.items():
            with self.subTest(extension=extension):
                self.load(content)
                output = self.root / f"render.{extension}"
                result = subprocess.run([sys.executable, "-B", str(SCRIPT), "--input", str(self.input),
                                         "--output", str(output), "--title", "CLI training test", "--dpi", "60"],
                                        capture_output=True, text=True, timeout=60)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(result.stderr, "")
                self.assertGreater(output.stat().st_size, 1000)
                if extension == "png":
                    image = plotter.plt.imread(output)
                    self.assertEqual(image.shape[:2], (552, 768))
                elif extension == "pdf":
                    self.assertTrue(output.read_bytes().startswith(b"%PDF"))
                else:
                    self.assertIn("<svg", output.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
