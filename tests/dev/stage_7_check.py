#!/usr/bin/env python3
"""Stage-7 developer self-check (handoff/0007-words-cli-mcp/brief.md section 7).

Checks:
(a) Every word of the catalogue parses and prints back to the same canonical text;
(b) Five malformed or unsupported inputs are refused with the expected messages
    (unknown operator, rank 3, mixed basis without transform, U(t) without time,
    unknown key in the run file);
(c) runs/goldens/g1.json and runs/goldens/g5.json through the CLI give files
    byte-identical to G1 and G5 (cmp);
(d) The MCP server started in a subprocess answers words_describe and inspect
    with the same JSON as the CLI.
"""

from __future__ import annotations

import asyncio
import filecmp
import glob
import json
import os
import subprocess
import sys
import tempfile

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if os.path.join(REPO, "python") not in sys.path:
    sys.path.insert(0, os.path.join(REPO, "python"))

from lsquant.runfile import RunFile  # noqa: E402
from lsquant.words import CATALOGUE, compile as compile_word, describe, parse  # noqa: E402


def check_a_catalogue_roundtrip() -> None:
    print("== (a) Catalogue canonical parse/print round-trip ==")
    words_to_test = [entry.canonical_example for entry in CATALOGUE] + [
        "Tr[sea[VX delta(E-H) VY dG-(E)/dE]]",
        "Tr[P:PX U(t)^+ SX U(t) delta(E-H) P:PX]",
        "Tr[B^+ VX@bloch B delta(E-H)]",
    ]
    for w_text in words_to_test:
        parsed = parse(w_text)
        canonical = str(parsed)
        if canonical != w_text:
            raise AssertionError(
                f"Round-trip mismatch: input '{w_text}' -> '{canonical}'"
            )
        info = describe(w_text)
        assert info["word"]["canonical"] == w_text
        print(
            f"  PASS  {w_text:<42} -> shape={info['catalogue']['shape_id']:<20} "
            f"plan_builder={info['catalogue']['plan_builder']}"
        )


def check_b_refusals() -> None:
    print("\n== (b) Refusal of 5 malformed / unsupported inputs ==")

    cases = [
        (
            "unknown_operator",
            lambda: compile_word(
                "Tr[FOO delta(E-H)]",
                basis={"name": "real", "operators": ["VX", "VY"]},
                params={"moments": 64, "states": 1},
            ),
            "Unknown operator 'FOO'",
        ),
        (
            "rank_3",
            lambda: compile_word(
                "Tr[VX delta(E-H) VY delta(E-H) VZ delta(E-H)]",
                basis="real",
                params={"moments": 64, "states": 1},
            ),
            "rank 3",
        ),
        (
            "mixed_basis_without_transform",
            lambda: compile_word(
                "Tr[VX@bloch delta(E-H) VX@real delta(E-H)]",
                basis="real",
                params={"moments": 64, "states": 1},
            ),
            "Mixed basis without explicit transform",
        ),
        (
            "ut_without_time",
            lambda: compile_word(
                "Tr[A delta(E-H) U(t)^+ B U(t)]",
                basis="real",
                params={"moments": 64, "states": 1},
            ),
            "without 'time' parameters",
        ),
        (
            "unknown_key_in_runfile",
            lambda: RunFile.from_dict(
                {
                    "model": {
                        "kind": "csr",
                        "dir": "goldens/stage0/chain1d/operators",
                        "label": "chain1d",
                    },
                    "word": "Tr[delta(E-H)]",
                    "expansion": {"moments": 64},
                    "unknown_extra_key": 123,
                },
                repo_root=REPO,
            ),
            "Unknown key 'unknown_extra_key' in run file",
        ),
    ]

    for label, fn, expected_substr in cases:
        try:
            fn()
        except ValueError as exc:
            msg = str(exc)
            if expected_substr not in msg:
                raise AssertionError(
                    f"Case '{label}' raised ValueError('{msg}'), "
                    f"missing expected substring '{expected_substr}'"
                ) from exc
            print(f"  PASS  [{label}] refused with: {msg}")
        else:
            raise AssertionError(
                f"Case '{label}' did NOT raise ValueError (expected '{expected_substr}')"
            )

    # Also verify all committed run files under runs/goldens/*.json validate cleanly
    rf_paths = sorted(glob.glob(os.path.join(REPO, "runs", "goldens", "*.json")))
    assert rf_paths, "No run files found in runs/goldens/"
    for p in rf_paths:
        RunFile.from_file(p, repo_root=REPO)
    print(f"  PASS  [runfile_schema] all {len(rf_paths)} files in runs/goldens/*.json validated")


def _run_cli(args: list[str]) -> subprocess.CompletedProcess[str]:
    env = dict(os.environ, KPM_SEED="12345")
    return subprocess.run(
        [sys.executable, "-m", "lsquant", *args],
        cwd=REPO,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )


def _cmp_file(new_path: str, gold_path: str, tag: str) -> None:
    if not os.path.isfile(new_path):
        raise AssertionError(f"DIFFERS {tag}: missing produced file {new_path}")
    if not os.path.isfile(gold_path):
        raise AssertionError(f"DIFFERS {tag}: missing golden file {gold_path}")
    r = subprocess.run(["cmp", "-s", new_path, gold_path], check=False)
    ok = (r.returncode == 0) and filecmp.cmp(new_path, gold_path, shallow=False)
    if not ok:
        raise AssertionError(f"DIFFERS {tag}: {new_path} != {gold_path}")
    print(f"  PASS  cmp IDENTICAL: {tag}")


def check_c_cli_goldens(tmp_root: str) -> str:
    print("\n== (c) CLI reproduction of G1 (runs/goldens/g1.json) and G5 (runs/goldens/g5.json) ==")
    out_g1 = os.path.join(tmp_root, "g1")
    r1 = _run_cli(["run", "runs/goldens/g1.json", "--out", out_g1])
    if r1.returncode != 0:
        raise RuntimeError(f"CLI g1.json failed ({r1.returncode}): {r1.stderr}")

    g1_name = "SpectralOp1chain1dKPM_M64_stateexact.chebmom1D"
    g1_new = os.path.join(out_g1, g1_name)
    g1_gold = os.path.join(REPO, "goldens", "stage0", "chain1d", g1_name)
    _cmp_file(g1_new, g1_gold, f"G1 ({g1_name})")

    out_g5 = os.path.join(tmp_root, "g5")
    r5 = _run_cli(["run", "runs/goldens/g5.json", "--out", out_g5])
    if r5.returncode != 0:
        raise RuntimeError(f"CLI g5.json failed ({r5.returncode}): {r5.stderr}")

    g5_files = [
        "NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D",
        "NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom2D",
        "KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat",
        "KuboBastinSea_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat",
        "KuboBastinSurf_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat",
        "KuboGreenWood_NonEqOpVX-VXhaldaneKPM_M128x128_stateexactJACKSON.dat",
    ]
    for fname in g5_files:
        _cmp_file(
            os.path.join(out_g5, fname),
            os.path.join(REPO, "goldens", "stage0", "haldane", fname),
            f"G5 ({fname})",
        )
    return g1_new


async def _run_mcp_check_async(g1_gold_path: str) -> None:
    from mcp import ClientSession, StdioServerParameters
    from mcp.client.stdio import stdio_client

    word_query = "Tr[VX delta(E-H) VY dG-(E)/dE]"

    cli_desc_proc = _run_cli(["words", "describe", word_query])
    if cli_desc_proc.returncode != 0:
        raise RuntimeError(f"CLI words describe failed: {cli_desc_proc.stderr}")
    cli_desc_json = json.loads(cli_desc_proc.stdout)

    cli_insp_proc = _run_cli(["inspect", g1_gold_path])
    if cli_insp_proc.returncode != 0:
        raise RuntimeError(f"CLI inspect failed: {cli_insp_proc.stderr}")
    cli_insp_json = json.loads(cli_insp_proc.stdout)

    server_params = StdioServerParameters(
        command=sys.executable,
        args=["-m", "lsquant.mcp_server", "--root", REPO],
        env=dict(os.environ, LSQUANT_MCP_ROOT=REPO),
    )

    async with stdio_client(server_params) as (read_stream, write_stream):
        async with ClientSession(read_stream, write_stream) as session:
            await session.initialize()
            tools_resp = await session.list_tools()
            tool_names = [t.name for t in tools_resp.tools]
            print(f"  MCP tools registered: {tool_names}")

            # 1. Call words_describe
            req_1 = {"tool": "words_describe", "arguments": {"word": word_query}}
            print(f"  MCP -> call_tool: {json.dumps(req_1)}")
            res_1 = await session.call_tool("words_describe", {"word": word_query})
            mcp_desc_text = res_1.content[0].text
            mcp_desc_json = json.loads(mcp_desc_text)
            print("  MCP <- words_describe response:")
            for line in mcp_desc_text.splitlines():
                print(f"    {line}")
            if mcp_desc_json != cli_desc_json:
                raise AssertionError(
                    "MCP words_describe JSON differs from CLI words describe JSON"
                )
            print("  PASS  MCP words_describe JSON == CLI words describe JSON")

            # 2. Call inspect
            rel_g1 = os.path.relpath(g1_gold_path, REPO)
            req_2 = {"tool": "inspect", "arguments": {"path": rel_g1}}
            print(f"  MCP -> call_tool: {json.dumps(req_2)}")
            res_2 = await session.call_tool("inspect", {"path": rel_g1})
            mcp_insp_text = res_2.content[0].text
            mcp_insp_json = json.loads(mcp_insp_text)
            print("  MCP <- inspect response:")
            for line in mcp_insp_text.splitlines():
                print(f"    {line}")
            if mcp_insp_json != cli_insp_json:
                raise AssertionError(
                    "MCP inspect JSON differs from CLI inspect JSON"
                )
            print("  PASS  MCP inspect JSON == CLI inspect JSON")


def check_d_mcp_subprocess() -> None:
    print("\n== (d) MCP server subprocess vs CLI JSON equivalence ==")
    g1_gold = os.path.join(
        REPO,
        "goldens",
        "stage0",
        "chain1d",
        "SpectralOp1chain1dKPM_M64_stateexact.chebmom1D",
    )
    asyncio.run(_run_mcp_check_async(g1_gold))


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="stage7_check_") as tmp_root:
        check_a_catalogue_roundtrip()
        check_b_refusals()
        check_c_cli_goldens(tmp_root)
        check_d_mcp_subprocess()
    print("\nSTAGE 7 SELF-CHECK PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
