"""Guards on the Paper 1 style constants.

`resistance_cmap()` once returned unreversed viridis while its docstring
claimed the opposite direction, so Figure 2 and the dashboard rendered the
same data with inverted colour meanings. Nothing failed, because nothing
asserted the direction. These tests do.
"""

from __future__ import annotations

from hivdrwatch.paper1.style import COLORBLIND_PALETTE, resistance_cmap


def _luminance(rgba: tuple[float, float, float, float]) -> float:
    r, g, b, _ = rgba
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def test_resistance_cmap_is_darker_at_higher_prevalence() -> None:
    cmap = resistance_cmap()
    assert _luminance(cmap(1.0)) < _luminance(cmap(0.0))


def test_resistance_cmap_is_monotonic_in_luminance() -> None:
    cmap = resistance_cmap()
    lums = [_luminance(cmap(i / 10)) for i in range(11)]
    assert lums == sorted(lums, reverse=True)


def test_palette_is_the_okabe_ito_eight() -> None:
    assert len(COLORBLIND_PALETTE) == 8
    assert COLORBLIND_PALETTE[0] == "#0072B2"
