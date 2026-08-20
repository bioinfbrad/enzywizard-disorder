from __future__ import annotations
from argparse import Namespace, ArgumentParser
import sys
from ..services.disorder_service import run_disorder_service


def add_disorder_parser(parser:ArgumentParser) -> None:
    parser.add_argument("-i", "--input_path",required=True,help="Path to the input cleaned protein structure file in CIF or PDB format.")
    parser.add_argument("-o", "--output_dir",required=True,help="Path to the output directory for saving the JSON report.")
    parser.add_argument("--window_size",type=int,default=11,help="Sliding window size for disorder score calculation (default: 11). Larger values produce smoother regional trends, while smaller values are more sensitive to local fluctuations.")
    parser.add_argument("--min_region_length",type=int,default=5,help="Minimum number of consecutive residues required to define a disordered region (default: 5). Shorter predicted segments will be ignored.")

    parser.set_defaults(func=run_disorder)


def run_disorder(args: Namespace) -> None:
    success = run_disorder_service(input_path=args.input_path,output_dir=args.output_dir,window_size=args.window_size,min_region_length=args.min_region_length)
    if not success:
        sys.exit(1)
