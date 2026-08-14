"""Unit tests for the Fig 1 pipeline architecture spec.

The diagram is a data structure — we test its shape, connectedness,
and per-phase counts so a bad edit doesn't silently produce a broken
diagram. The matplotlib rendering is exercised by the notebook.
"""

from __future__ import annotations

from hivdrwatch.paper1.fig01_pipeline import (
    PHASE_COLORS,
    Block,
    Edge,
    PipelineDiagram,
    pipeline_diagram,
    stage_summary,
)


class TestDiagramShape:
    def test_returns_populated_diagram(self) -> None:
        d = pipeline_diagram()
        assert isinstance(d, PipelineDiagram)
        assert len(d.blocks) >= 8
        assert len(d.edges) >= 8
        assert len(d.row_labels) >= 3

    def test_block_keys_are_unique(self) -> None:
        keys = [b.key for b in pipeline_diagram().blocks]
        assert len(keys) == len(set(keys)), f"duplicate block keys: {keys}"

    def test_every_edge_references_real_blocks(self) -> None:
        d = pipeline_diagram()
        keys = {b.key for b in d.blocks}
        for e in d.edges:
            assert e.src in keys, f"edge from unknown block: {e.src}"
            assert e.dst in keys, f"edge to unknown block: {e.dst}"

    def test_every_block_has_valid_phase(self) -> None:
        for b in pipeline_diagram().blocks:
            assert b.phase in PHASE_COLORS, f"unknown phase: {b.phase}"

    def test_every_phase_has_a_color(self) -> None:
        for phase, color in PHASE_COLORS.items():
            assert color.startswith("#"), f"bad color for phase {phase!r}: {color}"
            assert len(color) == 7


class TestConnectivity:
    def test_source_blocks_are_at_top_row(self) -> None:
        d = pipeline_diagram()
        sources = [b for b in d.blocks if b.phase == "source"]
        assert sources
        for b in sources:
            assert b.row == 0

    def test_serve_blocks_are_at_bottom(self) -> None:
        d = pipeline_diagram()
        serve = [b for b in d.blocks if b.phase == "serve"]
        assert serve
        max_row = max(b.row for b in d.blocks)
        for b in serve:
            assert b.row == max_row

    def test_no_upward_edges(self) -> None:
        """Data flows top-to-bottom; an edge going up is a spec error."""
        d = pipeline_diagram()
        row_of = {b.key: b.row for b in d.blocks}
        for e in d.edges:
            assert row_of[e.dst] >= row_of[e.src], (
                f"edge {e.src} → {e.dst} goes upward (rows {row_of[e.src]} → {row_of[e.dst]})"
            )

    def test_ingest_is_the_only_bronze_producer(self) -> None:
        """Every bronze-tier edge must originate from `ingest`."""
        d = pipeline_diagram()
        for e in d.edges:
            if e.tier == "bronze":
                assert e.src == "ingest", f"bronze edge from non-ingest block: {e.src}"

    def test_every_non_source_block_is_reachable_from_a_source(self) -> None:
        d = pipeline_diagram()
        sources = {b.key for b in d.blocks if b.phase == "source"}
        # Simple BFS forward reachability.
        reached = set(sources)
        added = True
        while added:
            added = False
            for e in d.edges:
                if e.src in reached and e.dst not in reached:
                    reached.add(e.dst)
                    added = True
        for b in d.blocks:
            assert b.key in reached, f"unreachable block: {b.key}"


class TestStageSummary:
    def test_shape(self) -> None:
        summary = stage_summary()
        assert set(summary) == set(PHASE_COLORS)
        assert sum(summary.values()) == len(pipeline_diagram().blocks)

    def test_all_phases_populated(self) -> None:
        summary = stage_summary()
        for phase, count in summary.items():
            assert count >= 1, f"empty phase: {phase}"


class TestConstants:
    def test_phase_order_locked(self) -> None:
        expected = ("source", "ingest", "call", "curate", "analyze", "serve")
        assert tuple(PHASE_COLORS) == expected

    def test_block_is_hashable(self) -> None:
        # frozen dataclass → usable as dict key / set member
        b = Block(
            key="k",
            label="L",
            module="m",
            phase="source",
            row=0,
            col=0.0,
        )
        assert b in {b}

    def test_edge_is_hashable(self) -> None:
        e = Edge(src="a", dst="b", tier="bronze")
        assert e in {e}
