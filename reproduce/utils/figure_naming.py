from __future__ import annotations

from pathlib import Path

# Mapping from internal publish source names to paper figure labels.
# Keys are grouped by publish section to avoid name collisions.
_SECTION6_STEM_MAP: dict[str, str] = {
    # Figure 5: k-regular graphs
    "3-regular-graph_time_comparison": "Figure_5_a",
    "4-regular-graph_time_comparison": "Figure_5_b",
    "5-regular-graph_time_comparison": "Figure_5_c",
    # Figure 6: k-regular l-colored graphs
    "3-regular-2-colored-graph_time_comparison": "Figure_6_a",
    "3-regular-3-colored-graph_time_comparison": "Figure_6_b",
    "3-regular-4-colored-graph_time_comparison": "Figure_6_c",
    "4-regular-2-colored-graph_time_comparison": "Figure_6_d",
    "4-regular-3-colored-graph_time_comparison": "Figure_6_e",
    "5-regular-2-colored-graph_time_comparison": "Figure_6_f",
    # Figure 7: directed graphs
    "2-regular-directed-graph_time_comparison": "Figure_7_a",
    "3-regular-directed-graph_time_comparison": "Figure_7_b",
    # Figure 8: BA graphs
    "BA_CC_time_comparison": "Figure_8_a",
    "BA_time_comparison": "Figure_8_b",
    # Figure 9: modulo-counting graph families
    "0mod2-regular-graph_time_comparison": "Figure_9_a",
    "1mod2-regular-graph_time_comparison": "Figure_9_b",
    "2mod4-regular-graph_time_comparison": "Figure_9_c",
    # Figure 10: odd-degree runtime
    "m-odd-degree-graph-sc2_m2_time_comparison": "Figure_10_a",
    "m-odd-degree-graph-sc2_m4_time_comparison": "Figure_10_b",
    "m-odd-degree-graph-sc2_m6_time_comparison": "Figure_10_c",
}

_APPENDIX_B1_STEM_MAP: dict[str, str] = {
    # Figure 11: correctness on k-regular graphs
    "2-regular-graph_correctness": "Figure_11_a",
    "3-regular-graph_correctness": "Figure_11_b",
    "4-regular-graph_correctness": "Figure_11_c",
    # Figure 12: correctness on modulo-counting regular graphs
    "0mod2-regular-graph_correctness": "Figure_12_a",
    "1mod2-regular-graph_correctness": "Figure_12_b",
    "2mod4-regular-graph_correctness": "Figure_12_c",
    # Figure 13: correctness on odd-degree graphs
    "odd_degree_n": "Figure_13_a",
    "odd_degree_k": "Figure_13_b",
    "odd_degree_m": "Figure_13_c",
}

_APPENDIX_B2_STEM_MAP: dict[str, str] = {
    # Figure 14: memory on k-regular graphs
    "3-regular-graph_memory_comparison": "Figure_14_a",
    "4-regular-graph_memory_comparison": "Figure_14_b",
    "5-regular-graph_memory_comparison": "Figure_14_c",
    # Figure 15: memory on colored graph settings
    "3-regular-2-colored-graph_memory_comparison": "Figure_15_a",
    "3-regular-4-colored-graph_memory_comparison": "Figure_15_b",
    "5-regular-2-colored-graph_memory_comparison": "Figure_15_c",
    # Figure 16: memory on directed graphs
    "2-regular-directed-graph_memory_comparison": "Figure_16_a",
    "3-regular-directed-graph_memory_comparison": "Figure_16_b",
}

_PUBLISH_STEM_MAP_BY_SECTION: dict[str, dict[str, str]] = {
    "Section6": _SECTION6_STEM_MAP,
    "AppendixB.1": _APPENDIX_B1_STEM_MAP,
    "AppendixB.2": _APPENDIX_B2_STEM_MAP,
}


def resolve_publish_figure_stem(section: str, source_stem: str) -> str | None:
    """Return mapped paper-style figure stem for a publish target, if available."""
    section_map = _PUBLISH_STEM_MAP_BY_SECTION.get(section)
    if section_map is None:
        return None
    return section_map.get(source_stem)


def resolve_publish_figure_filename(section: str, source_filename: str) -> str | None:
    """Map source filename to paper-style publish filename.

    Returns None when no mapping exists for the given section/source pair.
    """
    source_path = Path(source_filename)
    suffix = source_path.suffix or ".pdf"

    mapped_stem = resolve_publish_figure_stem(section=section, source_stem=source_path.stem)
    if mapped_stem is None:
        return None

    return f"{mapped_stem}{suffix}"
