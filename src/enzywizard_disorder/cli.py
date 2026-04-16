from __future__ import annotations

import argparse

from .commands.disorder import add_disorder_parser


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="enzywizard-disorder",
        description="EnzyWizard-Disorder: Predict intrinsically disordered regions from a cleaned protein structure and generating a detailed JSON report."
    )
    add_disorder_parser(parser)
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)