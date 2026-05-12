from __future__ import annotations
from ..utils.logging_utils import Logger
from typing import Any, Dict, List

def moving_average(values: List[float], window_size: int, logger: Logger) -> List[float] | None:
    if window_size <= 0:
        logger.print("[ERROR] window_size must be positive")
        return None

    n = len(values)
    if n == 0:
        logger.print("[ERROR] Empty input for moving_average")
        return None

    half = window_size // 2
    result: List[float] = []

    for i in range(n):
        left = max(0, i - half)
        right = min(n, i + half + 1)
        sub = values[left:right]
        result.append(sum(sub) / len(sub))

    return result



def postprocess_disorder_report_to_schema(
    raw_report: Dict[str, Any],
    logger: Logger,
) -> Dict[str, Any] | None:
    """
    Convert the raw EnzyWizard-Disorder report into the schema-compliant
    EnzyWizard-Disorder JSON report format.

    """

    if not isinstance(raw_report, dict):
        logger.print("[ERROR] raw_report must be a dictionary.")
        return None

    output_type = raw_report.get("output_type")
    if output_type != "enzywizard_disorder":
        logger.print(
            f"[ERROR] Invalid output_type in disorder report: {output_type}. "
            "Expected 'enzywizard_disorder'."
        )
        return None

    raw_statistics = raw_report.get("disorder_region_statistics")
    if not isinstance(raw_statistics, dict):
        logger.print("[ERROR] Missing or invalid disorder_region_statistics in disorder report.")
        return None

    region_num = raw_statistics.get("region_num")
    max_region_length = raw_statistics.get("max_region_length")
    total_region_length = raw_statistics.get("total_region_length")

    if not isinstance(region_num, int):
        logger.print("[ERROR] disorder_region_statistics.region_num must be an integer.")
        return None

    if not isinstance(max_region_length, int):
        logger.print("[ERROR] disorder_region_statistics.max_region_length must be an integer.")
        return None

    if not isinstance(total_region_length, int):
        logger.print("[ERROR] disorder_region_statistics.total_region_length must be an integer.")
        return None

    if region_num < 0:
        logger.print("[ERROR] disorder_region_statistics.region_num must be non-negative.")
        return None

    if max_region_length < 0:
        logger.print("[ERROR] disorder_region_statistics.max_region_length must be non-negative.")
        return None

    if total_region_length < 0:
        logger.print("[ERROR] disorder_region_statistics.total_region_length must be non-negative.")
        return None

    raw_regions = raw_report.get("disorder_regions")
    if not isinstance(raw_regions, list):
        logger.print("[ERROR] Missing or invalid disorder_regions in disorder report.")
        return None

    disordered_regions: List[Dict[str, Any]] = []

    for region_index, raw_region in enumerate(raw_regions):
        if not isinstance(raw_region, dict):
            logger.print(f"[ERROR] Invalid disorder region at index {region_index}.")
            return None

        region_length = raw_region.get("length")
        if not isinstance(region_length, int):
            logger.print(
                f"[ERROR] disorder_regions[{region_index}].length must be an integer."
            )
            return None

        if region_length < 0:
            logger.print(
                f"[ERROR] disorder_regions[{region_index}].length must be non-negative."
            )
            return None

        raw_residues = raw_region.get("residues")
        if not isinstance(raw_residues, list):
            logger.print(
                f"[ERROR] disorder_regions[{region_index}].residues must be a list."
            )
            return None

        residues: List[Dict[str, Any]] = []

        for residue_index_in_region, raw_residue in enumerate(raw_residues):
            if not isinstance(raw_residue, dict):
                logger.print(
                    f"[ERROR] Invalid residue at disorder_regions[{region_index}]"
                    f".residues[{residue_index_in_region}]."
                )
                return None

            aa_id = raw_residue.get("aa_id")
            aa_name = raw_residue.get("aa_name")

            if not isinstance(aa_id, int):
                logger.print(
                    f"[ERROR] disorder_regions[{region_index}]"
                    f".residues[{residue_index_in_region}].aa_id must be an integer."
                )
                return None

            if not isinstance(aa_name, str):
                logger.print(
                    f"[ERROR] disorder_regions[{region_index}]"
                    f".residues[{residue_index_in_region}].aa_name must be a string."
                )
                return None

            if len(aa_name) != 1 or aa_name not in "ACDEFGHIKLMNPQRSTVWY":
                logger.print(
                    f"[ERROR] disorder_regions[{region_index}]"
                    f".residues[{residue_index_in_region}].aa_name must be a valid "
                    "one-letter amino acid code."
                )
                return None

            residues.append({
                "residue_index": aa_id,
                "residue_name": aa_name,
            })

        disordered_regions.append({
            "disordered_region_length": region_length,
            "residues": residues,
        })

    schema_report: Dict[str, Any] = {
        "report_type": "enzywizard_disorder",
        "disordered_region_statistics": {
            "disordered_region_count": region_num,
            "max_disordered_region_length": max_region_length,
            "total_disordered_region_length": total_region_length,
        },
        "disordered_regions": disordered_regions,
    }

    return schema_report