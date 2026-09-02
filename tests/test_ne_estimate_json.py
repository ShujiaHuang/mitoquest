#!/usr/bin/env python3
"""Regression test for strict parsing of ne-estimate's JSON object."""

import json
import pathlib
import subprocess
import sys
import tempfile


def main():
    if len(sys.argv) != 2:
        raise SystemExit("usage: test_ne_estimate_json.py <mitoquest>")

    mitoquest = sys.argv[1]
    with tempfile.TemporaryDirectory(prefix="mitoquest-ne-json-") as directory:
        directory = pathlib.Path(directory)
        input_tsv = directory / "pairs.tsv"
        output_json = directory / "result.json"
        family_id = 'F"\\family'
        mother_id = 'M"\\mother'
        input_tsv.write_text(
            "CHROM\tPOS\tREF\tALT\tFAM_ID\tMOTHER_ID\tCHILD_ID\t"
            "MOTHER_DP\tMOTHER_AD_ALT\tMOTHER_VAF\tCHILD_DP\t"
            "CHILD_AD_ALT\tQC\n"
            f"chrM\t100\tA\tG\t{family_id}\t{mother_id}\tC1\t"
            "100\t30\t0.3\t100\t20\tPASS\n"
        )
        subprocess.run(
            [
                mitoquest,
                "ne-estimate",
                "-i",
                str(input_tsv),
                "--model",
                "continuous",
                "--min-vaf",
                "0",
                "--max-vaf",
                "1",
                "--min-ne",
                "2",
                "--max-ne",
                "2",
                "--per-family",
                "--min-family-sites",
                "1",
                "--cross-check",
                "kimura",
                "--kimura-bootstrap",
                "10",
                "-o",
                str(output_json),
            ],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        lines = output_json.read_text().splitlines()
        assert lines[0].startswith("#mitoquest_ne_estimate_command=")
        payload = json.loads("\n".join(lines[1:]))
        family = payload["Per_Family_Estimates"][0]
        assert family["FAM_ID"] == family_id
        assert family["Mother_ID"] == mother_id
        assert family["Kimura_Cross_Check"]["N_Bootstrap"] == 10
        assert "Ne_Kimura_CI_95_Low" in family["Kimura_Cross_Check"]

        base_command = [
            mitoquest,
            "ne-estimate",
            "-i",
            str(input_tsv),
            "--min-ne",
            "2",
            "--max-ne",
            "2",
        ]
        for option in ("--min-vaf", "--max-vaf", "--kimura-trim", "--ne-profile-step"):
            rejected = subprocess.run(
                [*base_command, option, "nan"],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            assert rejected.returncode != 0, option
            assert "finite" in rejected.stderr.lower(), rejected.stderr

        impossible = subprocess.run(
            [
                mitoquest,
                "ne-estimate",
                "-i",
                str(input_tsv),
                "--min-vaf",
                "0",
                "--max-vaf",
                "1",
                "--min-ne",
                "1",
                "--max-ne",
                "1",
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        assert impossible.returncode != 0
        assert "no finite continuous-model likelihood" in impossible.stderr.lower()

        mixed_model = subprocess.run(
            [*base_command, "--model", "discrete", "--per-family"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        assert mixed_model.returncode != 0
        assert "per-family" in mixed_model.stderr.lower()


if __name__ == "__main__":
    main()