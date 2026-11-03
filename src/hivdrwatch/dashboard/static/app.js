/* SAWA resistance summary: draws summary.json as inline SVG.
 *
 * No framework and no network beyond the one same-origin fetch. The order
 * constants below are copies of paper1/style.py and dashboard/export/aggregate.py
 * (the page cannot import Python); tests/dashboard/test_static.py keeps them in step.
 * The landscape ramp is built from the same Okabe-Ito blues as paper1/style.py.
 */

const SUMMARY_URL = "./summary.json";
const SVG_NS = "http://www.w3.org/2000/svg";

const DRUG_CLASS_ORDER = Object.freeze(["PI", "NRTI", "NNRTI", "INSTI"]);
const GENE_ORDER = Object.freeze(["PR", "RT", "IN"]);
const SDRM_GENES = Object.freeze(["PR", "RT"]);
const LEVEL_ORDER = Object.freeze(["susceptible", "potential_low", "low", "intermediate", "high"]);

const LEVEL_LABEL = Object.freeze({
  susceptible: "Susceptible",
  potential_low: "Potential low-level",
  low: "Low-level",
  intermediate: "Intermediate",
  high: "High-level",
});
const CLASS_GENE = Object.freeze({ PI: "PR", NRTI: "RT", NNRTI: "RT", INSTI: "IN" });
const SDRM_TOP_N = 12;

const NUMBER = new Intl.NumberFormat("en-US");

/* ---- small helpers ------------------------------------------------------ */

const fmtInt = (n) => NUMBER.format(n);

function fmtPct(frac) {
  if (frac > 0 && frac < 0.001) return "<0.1%";
  return `${(frac * 100).toFixed(1)}%`;
}

function fmtShare(part, whole) {
  return whole > 0 ? `${((part / whole) * 100).toFixed(1)}%` : "0.0%";
}

function svgEl(tag, attrs = {}, text = null) {
  const el = document.createElementNS(SVG_NS, tag);
  for (const [name, value] of Object.entries(attrs)) el.setAttribute(name, String(value));
  if (text !== null) el.textContent = text;
  return el;
}

function newSvg(width, height, label) {
  const root = svgEl("svg", { viewBox: `0 0 ${width} ${height}`, role: "group", "aria-label": label });
  return root;
}

function withTitle(el, text) {
  el.append(svgEl("title", {}, text));
  return el;
}

function niceMax(maxFrac) {
  const step = 0.05;
  return Math.max(step, Math.ceil(maxFrac / step - 1e-9) * step);
}

/* ---- schema gate -------------------------------------------------------- */

function schemaMismatch(expected, actual) {
  return typeof actual !== "string" || actual !== expected;
}

// A bound sentence with its numbers missing reads as a claim, so withhold it along with the charts.
function showBanner(message) {
  const banner = document.getElementById("schema-banner");
  banner.textContent = message;
  banner.hidden = false;
  for (const el of document.querySelectorAll("main > section, [data-bind]")) {
    (el.closest("p, section") ?? el).hidden = true;
  }
}

/* ---- input digests ------------------------------------------------------ */

function digestsDiffer(meta) {
  return meta.input_file_md5 !== meta.input_dvc_md5 || meta.curated_file_md5 !== meta.curated_dvc_md5;
}

function syncDigestNote(meta) {
  document.getElementById("digest-note").hidden = !digestsDiffer(meta);
}

/* ---- text bindings ------------------------------------------------------ */

function median(values) {
  if (values.length === 0) return 0;
  const sorted = [...values].sort((a, b) => a - b);
  const mid = Math.floor(sorted.length / 2);
  return sorted.length % 2 ? sorted[mid] : (sorted[mid - 1] + sorted[mid]) / 2;
}

function listNames(names) {
  return new Intl.ListFormat("en", { style: "long", type: "conjunction" }).format(names);
}

function absentCountries(meta, landscape) {
  return meta.landscape_countries.filter(
    (country) => !landscape.some((row) => row.country === country.code && row.frac_resistant !== null),
  );
}

function deriveFigures(summary) {
  const { meta, views } = summary;
  const landscape = views.resistance_landscape;
  const drawnSizes = landscape.filter((row) => row.frac_resistant !== null).map((row) => row.n_sequences);
  return {
    landscape_max_cell_n: drawnSizes.length ? Math.max(...drawnSizes) : 0,
    landscape_median_cell_n: median(drawnSizes),
    landscape_region_names: listNames(meta.landscape_countries.map((country) => country.name)),
    landscape_absent_countries: listNames(absentCountries(meta, landscape).map((country) => country.name)),
    pct_landscape_undated: fmtShare(
      meta.landscape_n_excluded_no_year,
      meta.landscape_n_accessions + meta.landscape_n_excluded_no_year,
    ),
    n_unscored: meta.n_accessions_ingested - meta.n_accessions_scored,
    n_drugs: new Set(views.levels_by_drug.map((row) => row.drug)).size,
    n_classes: new Set(views.levels_by_drug.map((row) => row.drug_class)).size,
    n_repeated_sequences: meta.n_accessions_scored - meta.n_distinct_sequences,
    pct_sdrm_multistate: fmtShare(meta.n_sdrm_observations_multistate, meta.n_sdrm_observations),
  };
}

// Bound values that may be null in a valid summary (Sierra can omit the publish date).
// Every other bound value is required, and null there is a malformed file, not a blank.
const NULLABLE_BINDINGS = Object.freeze(["meta.algorithm_publish_date"]);
const MISSING = "\u2014";

function findLandscapeCell(summary, cls, country, year) {
  return summary.views.resistance_landscape.find(
    (row) => row.drug_class === cls && row.country === country && row.collection_year === Number(year),
  );
}

// `cell.<class>.<country>.<year>.<n|pct>` quotes one heatmap cell, so a sentence that cites a
// cell cannot drift from the figure. A cell under the floor has no percentage and is an error.
function cellBinding([cls, country, year, field], summary) {
  const row = findLandscapeCell(summary, cls, country, year);
  if (row === undefined) return undefined;
  if (field === "n") return row.n_sequences;
  if (field === "pct") return row.frac_resistant === null ? null : fmtPct(row.frac_resistant);
  return undefined;
}

// `ratio.<class>.<countryA>.<yearA>.<countryB>.<yearB>`: A's share over B's, to the nearest whole.
function ratioBinding([cls, countryA, yearA, countryB, yearB], summary) {
  const a = findLandscapeCell(summary, cls, countryA, yearA)?.frac_resistant;
  const b = findLandscapeCell(summary, cls, countryB, yearB)?.frac_resistant;
  if (a === undefined || b === undefined || a === null || b === null || b === 0) return undefined;
  return Math.round(a / b);
}

function resolveBinding(path, summary, derived) {
  const [scope, ...rest] = path.split(".");
  let value;
  if (scope === "meta") value = summary.meta[rest[0]];
  else if (scope === "derived") value = derived[rest[0]];
  else if (scope === "gene") value = summary.views.gene_coverage.find((row) => row.gene === rest[0])?.n_sequences;
  else if (scope === "cell") value = cellBinding(rest, summary);
  else if (scope === "ratio") value = ratioBinding(rest, summary);
  if (value === undefined) throw new Error(`summary.json has no value for ${path}`);
  if (value === null && !NULLABLE_BINDINGS.includes(path)) {
    throw new Error(`summary.json has null for ${path}, which must always have a value`);
  }
  return value;
}

// Never `String(null)`: a missing value reads as a dash, not as the word "null".
function boundText(value, trunc) {
  if (value === null) return MISSING;
  if (typeof value === "number") return fmtInt(value);
  return trunc ? value.slice(0, Number(trunc)) : value;
}

function applyBindings(summary) {
  const derived = deriveFigures(summary);
  for (const el of document.querySelectorAll("[data-bind]")) {
    const raw = resolveBinding(el.dataset.bind, summary, derived);
    el.textContent = boundText(raw, el.dataset.trunc);
    if (el.dataset.href && raw !== null) el.setAttribute("href", el.dataset.href.replace("{value}", raw));
  }
}

/* ---- charts ------------------------------------------------------------- */

function renderPrevalence(rows, container) {
  const width = Math.max(container.clientWidth, 200);
  const rowHeight = 64;
  const barHeight = 14;
  const axisMax = niceMax(Math.max(...rows.map((row) => row.frac_resistant)));
  const svg = newSvg(width, rows.length * rowHeight + 22, "Resistance prevalence by drug class");

  rows.forEach((row, i) => {
    const top = i * rowHeight;
    const gene = CLASS_GENE[row.drug_class];
    const group = svgEl("g");
    group.append(svgEl("text", { x: 0, y: top + 14, class: "strong" }, row.drug_class));
    group.append(svgEl("text", { x: width, y: top + 14, class: "strong num", "text-anchor": "end" }, fmtPct(row.frac_resistant)));
    group.append(svgEl("rect", { x: 0, y: top + 22, width, height: barHeight, class: "track" }));
    group.append(
      svgEl("rect", {
        x: 0,
        y: top + 22,
        width: Math.max((row.frac_resistant / axisMax) * width, 1),
        height: barHeight,
        class: `bar-${row.drug_class}`,
      }),
    );
    group.append(
      svgEl(
        "text",
        { x: 0, y: top + 55, class: "muted num" },
        `${fmtInt(row.n_resistant)} of ${fmtInt(row.n_sequences)} sequences scored for ${gene}`,
      ),
    );
    withTitle(group, `${row.drug_class}: ${fmtInt(row.n_resistant)} of ${fmtInt(row.n_sequences)} sequences scored for ${gene} resistant`);
    svg.append(group);
  });

  svg.append(
    svgEl("text", { x: width, y: rows.length * rowHeight + 14, class: "muted", "text-anchor": "end" }, `Bars span 0 to ${Math.round(axisMax * 100)}%`),
  );
  return svg;
}

function groupLevels(rows) {
  const byDrug = new Map();
  for (const row of rows) {
    if (!byDrug.has(row.drug)) byDrug.set(row.drug, { drug: row.drug, drug_class: row.drug_class, levels: {} });
    byDrug.get(row.drug).levels[row.level] = row.n_sequences;
  }
  return DRUG_CLASS_ORDER.map((cls) => ({
    drug_class: cls,
    drugs: [...byDrug.values()].filter((d) => d.drug_class === cls),
  })).filter((group) => group.drugs.length > 0);
}

function renderLevels(rows, container) {
  const width = Math.max(container.clientWidth, 280);
  const labelWidth = 54;
  const valueWidth = 48;
  const barX = labelWidth;
  const barWidth = width - labelWidth - valueWidth - 8;
  const headingHeight = 34;
  const rowHeight = 21;
  const barHeight = 12;
  const stackOrder = [...LEVEL_ORDER].reverse();

  const groups = groupLevels(rows);
  const height = groups.reduce((total, g) => total + headingHeight + g.drugs.length * rowHeight + 8, 0);
  const svg = newSvg(width, height, "Resistance level by drug");

  let y = 0;
  for (const group of groups) {
    const totals = group.drugs.map((d) => LEVEL_ORDER.reduce((sum, level) => sum + (d.levels[level] ?? 0), 0));
    const low = Math.min(...totals);
    const high = Math.max(...totals);
    const range = low === high ? fmtInt(low) : `${fmtInt(low)} to ${fmtInt(high)}`;
    const gene = CLASS_GENE[group.drug_class];

    svg.append(svgEl("text", { x: 0, y: y + 22, class: "strong" }, group.drug_class));
    svg.append(svgEl("text", { x: labelWidth, y: y + 22, class: "muted num" }, `n = ${range}, scored for ${gene}`));
    y += headingHeight;

    group.drugs.forEach((drug, i) => {
      const total = totals[i];
      const row = svgEl("g");
      row.append(svgEl("text", { x: 0, y: y + 11, class: "num" }, drug.drug));
      row.append(svgEl("rect", { x: barX, y, width: barWidth, height: barHeight, class: "track" }));
      let x = barX;
      for (const level of stackOrder) {
        const count = drug.levels[level] ?? 0;
        const w = total > 0 ? (count / total) * barWidth : 0;
        if (w > 0) {
          const seg = svgEl("rect", { x, y, width: w, height: barHeight, class: `lvl-${level}` });
          withTitle(seg, `${drug.drug}, ${LEVEL_LABEL[level]}: ${fmtInt(count)} of ${fmtInt(total)} sequences (${fmtPct(count / total)})`);
          row.append(seg);
        }
        x += w;
      }
      const highCount = drug.levels.high ?? 0;
      row.append(svgEl("text", { x: width, y: y + 11, class: "num", "text-anchor": "end" }, fmtPct(total > 0 ? highCount / total : 0)));
      svg.append(row);
      y += rowHeight;
    });
    y += 8;
  }
  return svg;
}

function topSdrm(rows, gene) {
  return rows
    .filter((row) => row.panel === gene || row.panel.startsWith(`${gene}_`))
    .sort((a, b) => b.n_sequences - a.n_sequences || a.mutation.localeCompare(b.mutation))
    .slice(0, SDRM_TOP_N);
}

function renderSdrm(rows, gene, container) {
  const width = Math.max(container.clientWidth, 280);
  const labelWidth = 56;
  const valueWidth = 104;
  const barWidth = width - labelWidth - valueWidth - 8;
  const rowHeight = 22;
  const headingHeight = 30;

  const top = topSdrm(rows, gene);
  const denominator = top.length > 0 ? top[0].n_panel_denominator : 0;
  const axisMax = niceMax(Math.max(...top.map((row) => row.frac)));
  const svg = newSvg(width, headingHeight + top.length * rowHeight, `Most frequent surveillance mutations in ${gene}`);

  svg.append(svgEl("text", { x: 0, y: 16, class: "strong" }, gene));
  svg.append(svgEl("text", { x: labelWidth, y: 16, class: "muted num" }, `n = ${fmtInt(denominator)} scored; bars span 0 to ${Math.round(axisMax * 100)}%`));

  top.forEach((row, i) => {
    const y = headingHeight + i * rowHeight;
    const group = svgEl("g");
    group.append(svgEl("text", { x: 0, y: y + 11, class: "num" }, row.mutation));
    group.append(svgEl("rect", { x: labelWidth, y, width: barWidth, height: 12, class: "track" }));
    group.append(
      svgEl("rect", { x: labelWidth, y, width: Math.max((row.frac / axisMax) * barWidth, 1), height: 12, class: `bar-${row.panel}` }),
    );
    group.append(svgEl("text", { x: width, y: y + 11, class: "num", "text-anchor": "end" }, `${fmtInt(row.n_sequences)} · ${fmtPct(row.frac)}`));
    withTitle(group, `${gene} ${row.mutation}: ${fmtInt(row.n_sequences)} of ${fmtInt(row.n_panel_denominator)} sequences`);
    svg.append(group);
  });
  return svg;
}

function renderCoverage(rows, scored, container) {
  const width = Math.max(container.clientWidth, 280);
  const labelWidth = 40;
  const valueWidth = 112;
  const barWidth = width - labelWidth - valueWidth - 8;
  const rowHeight = 30;
  const svg = newSvg(width, GENE_ORDER.length * rowHeight, "Gene coverage of the scored cohort");

  GENE_ORDER.forEach((gene, i) => {
    const row = rows.find((r) => r.gene === gene);
    if (!row) throw new Error(`gene_coverage has no row for ${gene}`);
    const y = i * rowHeight + 6;
    const frac = row.n_sequences / scored;
    const group = svgEl("g");
    group.append(svgEl("text", { x: 0, y: y + 12, class: "strong" }, gene));
    group.append(svgEl("rect", { x: labelWidth, y, width: barWidth, height: 16, class: "track" }));
    group.append(svgEl("rect", { x: labelWidth, y, width: Math.max(frac * barWidth, 1), height: 16, class: "bar-gene" }));
    group.append(svgEl("text", { x: width, y: y + 12, class: "num", "text-anchor": "end" }, `${fmtInt(row.n_sequences)} · ${fmtPct(frac)}`));
    withTitle(group, `${gene}: ${fmtInt(row.n_sequences)} of ${fmtInt(scored)} scored accessions`);
    svg.append(group);
  });
  return svg;
}

/* ---- landscape heatmap -------------------------------------------------- */

// Light to dark, so darker always means a higher share. The middle two stops are Okabe-Ito sky
// blue and blue from paper1/style.py; the ends are a tint of the first and a shade of the second.
const RAMP_STOPS = Object.freeze(["#EAF4FB", "#56B4E9", "#0072B2", "#00263D"]);
const RAMP_STEPS = 20;
const LABEL_MIN_GAP = 34;
const MIN_WIDTH_FOR_N = 22;

function hexToRgb(hex) {
  const value = parseInt(hex.slice(1), 16);
  return [(value >> 16) & 255, (value >> 8) & 255, value & 255];
}

function rgbToHex(rgb) {
  return `#${rgb.map((c) => Math.round(c).toString(16).padStart(2, "0")).join("")}`;
}

function rampColor(frac) {
  const position = Math.min(1, Math.max(0, frac)) * (RAMP_STOPS.length - 1);
  const index = Math.min(Math.floor(position), RAMP_STOPS.length - 2);
  const t = position - index;
  const from = hexToRgb(RAMP_STOPS[index]);
  const to = hexToRgb(RAMP_STOPS[index + 1]);
  return rgbToHex(from.map((c, i) => c + (to[i] - c) * t));
}

function relativeLuminance(hex) {
  const [r, g, b] = hexToRgb(hex).map((c) => {
    const v = c / 255;
    return v <= 0.03928 ? v / 12.92 : ((v + 0.055) / 1.055) ** 2.4;
  });
  return 0.2126 * r + 0.7152 * g + 0.0722 * b;
}

// Whichever of white and near-black text has the higher contrast on the cell's fill.
function textTone(fillHex) {
  return relativeLuminance(fillHex) >= 0.211 ? "dark" : "light";
}

function fmtCellN(n) {
  return n >= 1000 ? `${(n / 1000).toFixed(1)}k` : String(n);
}

// A cell is coloured only at or above the floor. The export applies the floor and writes a null
// share below it; a cell where the two disagree is a malformed file, not something to draw.
function landscapeModel(rows, drugClass, countryCodes, floor) {
  const own = rows.filter((row) => row.drug_class === drugClass);
  for (const row of own) {
    if ((row.n_sequences < floor) !== (row.frac_resistant === null)) {
      throw new Error(`landscape cell ${row.country} ${row.collection_year} ${drugClass} disagrees with the floor of ${floor}`);
    }
  }
  const coloured = own.filter((row) => row.frac_resistant !== null);
  const years = [...new Set(coloured.map((row) => row.collection_year))].sort((a, b) => a - b);
  const countries = countryCodes.filter((code) => coloured.some((row) => row.country === code));
  const yearSet = new Set(years);
  const countrySet = new Set(countries);
  const cells = new Map();
  let omittedCells = 0;
  let omittedSequences = 0;
  for (const row of own) {
    if (countrySet.has(row.country) && yearSet.has(row.collection_year)) {
      cells.set(`${row.country}|${row.collection_year}`, row);
    } else {
      omittedCells += 1;
      omittedSequences += row.n_sequences;
    }
  }
  return { years, countries, cells, omittedCells, omittedSequences };
}

// Years are not spaced evenly: a run of years with no coloured cell is left out, and a narrow
// break marks the place so that two columns that are not neighbours do not look like they are.
function columnLayout(years, plotWidth, breakWidth) {
  const breaksBefore = (i) => i > 0 && years[i] - years[i - 1] > 1;
  const breaks = years.filter((_, i) => breaksBefore(i)).length;
  const cellWidth = (plotWidth - breaks * breakWidth) / Math.max(years.length, 1);
  const xs = new Map();
  let x = 0;
  years.forEach((year, i) => {
    if (breaksBefore(i)) x += breakWidth;
    xs.set(year, x);
    x += cellWidth;
  });
  return { cellWidth, xs };
}

// The first and last year always, then the multiples of five that do not crowd a label already placed.
function yearLabels(years, xs, cellWidth) {
  const centre = (year) => xs.get(year) + cellWidth / 2;
  const placed = [];
  const place = (year) => {
    if (placed.every((other) => Math.abs(centre(other) - centre(year)) >= LABEL_MIN_GAP)) placed.push(year);
  };
  if (years.length > 0) {
    place(years[0]);
    place(years[years.length - 1]);
  }
  for (const year of years) if (year % 5 === 0) place(year);
  return placed.sort((a, b) => a - b);
}

function describeLandscapeCell(row, name, drugClass, gene, floor) {
  const head = `${name} ${row.collection_year}, ${drugClass}`;
  const counts = `${fmtInt(row.n_resistant)} of ${fmtInt(row.n_sequences)} sequences scored for ${gene} are resistant`;
  if (row.frac_resistant === null) {
    return `${head}: ${counts}. Fewer than ${floor} sequences, so no percentage is given.`;
  }
  return `${head}: ${counts} (${fmtPct(row.frac_resistant)}).`;
}

function drawLandscapeCell(g, geometry, row, label) {
  const { x, y, w, h, showN } = geometry;
  if (row === undefined) {
    g.append(svgEl("rect", { x, y, width: w, height: h, class: "hm-empty" }));
    return withTitle(g, `${label}: no sequences`);
  }
  const sparse = row.frac_resistant === null;
  if (sparse) {
    g.append(svgEl("rect", { x, y, width: w, height: h, class: "hm-nodata" }));
    g.append(svgEl("line", { x1: x + 2, y1: y + h - 2, x2: x + w - 2, y2: y + 2, class: "hm-slash" }));
  } else {
    const fill = rampColor(row.frac_resistant);
    g.append(svgEl("rect", { x, y, width: w, height: h, class: "hm-fill", fill }));
  }
  if (showN) {
    const tone = sparse ? "nodata" : textTone(rampColor(row.frac_resistant));
    g.append(svgEl("text", { x: x + w / 2, y: y + h / 2 + 3.5, class: `hm-n hm-n-${tone}`, "text-anchor": "middle" }, fmtCellN(row.n_sequences)));
  }
  return g;
}

function renderLandscape(model, drugClass, ctx, container, onSelect) {
  const width = Math.max(container.clientWidth, 280);
  const labelWidth = 30;
  const cellHeight = 24;
  const pitch = 26;
  const axisHeight = 22;
  const breakWidth = 8;
  const names = new Map(ctx.countries.map((country) => [country.code, country.name]));
  const gene = CLASS_GENE[drugClass];
  const { cellWidth, xs } = columnLayout(model.years, width - labelWidth, breakWidth);
  const showN = cellWidth >= MIN_WIDTH_FOR_N;
  const height = Math.max(model.countries.length, 1) * pitch + axisHeight;
  const svg = newSvg(width, height, `Resistance prevalence for ${drugClass} by country and year`);

  if (model.countries.length === 0) {
    svg.append(svgEl("text", { x: 0, y: 18, class: "muted" }, `No ${drugClass} cell reaches ${ctx.floor} sequences.`));
    return svg;
  }

  model.countries.forEach((code, r) => {
    const y = r * pitch;
    svg.append(withTitle(svgEl("text", { x: 0, y: y + cellHeight / 2 + 4, class: "strong" }, code), names.get(code)));
    for (const year of model.years) {
      const row = model.cells.get(`${code}|${year}`);
      const label = `${names.get(code)} ${year}, ${drugClass}`;
      const g = svgEl("g", { class: "hm-cell" });
      const geometry = { x: labelWidth + xs.get(year), y, w: Math.max(cellWidth - 1, 1), h: cellHeight, showN };
      drawLandscapeCell(g, geometry, row, label);
      if (row !== undefined) {
        const text = describeLandscapeCell(row, names.get(code), drugClass, gene, ctx.floor);
        withTitle(g, text);
        g.addEventListener("click", () => onSelect(text));
      }
      svg.append(g);
    }
  });

  const axisY = model.countries.length * pitch + 14;
  for (const year of yearLabels(model.years, xs, cellWidth)) {
    const centre = labelWidth + xs.get(year) + cellWidth / 2;
    const atEdge = centre + 14 > width;
    svg.append(
      svgEl("text", { x: atEdge ? width : centre, y: axisY, class: "muted num", "text-anchor": atEdge ? "end" : "middle" }, year),
    );
  }
  return svg;
}

function renderLandscapeLegend(ctx, rows) {
  const legend = document.getElementById("legend-landscape");
  const ramp = svgEl("svg", { viewBox: `0 0 ${RAMP_STEPS * 6} 12`, width: RAMP_STEPS * 6, height: 12, role: "img", "aria-label": "Colour scale from 0% to 100%" });
  for (let i = 0; i < RAMP_STEPS; i += 1) {
    ramp.append(svgEl("rect", { x: i * 6, y: 0, width: 6, height: 12, fill: rampColor((i + 0.5) / RAMP_STEPS) }));
  }
  const scale = document.createElement("li");
  scale.append("Share resistant 0%", ramp, "100%");

  const sparse = document.createElement("li");
  const sparseSwatch = document.createElement("span");
  sparseSwatch.className = "swatch swatch-nodata";
  sparse.append(sparseSwatch, `Fewer than ${ctx.floor} sequences: not enough data`);

  const empty = document.createElement("li");
  const emptySwatch = document.createElement("span");
  emptySwatch.className = "swatch swatch-empty";
  empty.append(emptySwatch, "No sequences");

  const shown = new Set(rows.filter((row) => row.frac_resistant !== null).map((row) => row.country));
  const rowKey = document.createElement("li");
  rowKey.append(ctx.countries.filter((c) => shown.has(c.code)).map((c) => `${c.code} ${c.name}`).join(" \u00b7 "));

  legend.replaceChildren(scale, sparse, empty, rowKey);
}

function drawLandscape(summary) {
  const { views, meta } = summary;
  const ctx = { countries: meta.landscape_countries, floor: meta.landscape_min_cell_n };
  const codes = ctx.countries.map((country) => country.code);
  for (const drugClass of DRUG_CLASS_ORDER) {
    const container = document.getElementById(`chart-landscape-${drugClass}`);
    const readout = document.getElementById(`readout-landscape-${drugClass}`);
    const model = landscapeModel(views.resistance_landscape, drugClass, codes, ctx.floor);
    container.replaceChildren(renderLandscape(model, drugClass, ctx, container, (text) => {
      readout.textContent = text;
    }));
    const omitted = document.getElementById(`omitted-landscape-${drugClass}`);
    omitted.hidden = model.omittedCells === 0;
    const plural = model.omittedCells === 1 ? "cell" : "cells";
    omitted.textContent = `${fmtInt(model.omittedCells)} ${plural} under the floor, holding ${fmtInt(model.omittedSequences)} sequences, fall in years or countries with no coloured cell and are not drawn.`;
  }
}

function syncLandscapeAbsent(summary) {
  const absent = absentCountries(summary.meta, summary.views.resistance_landscape);
  document.getElementById("landscape-absent").hidden = absent.length === 0;
}

function renderLevelLegend() {
  const legend = document.getElementById("legend-levels");
  for (const level of [...LEVEL_ORDER].reverse()) {
    const item = document.createElement("li");
    const swatch = document.createElement("span");
    swatch.className = `swatch swatch-${level}`;
    item.append(swatch, LEVEL_LABEL[level]);
    legend.append(item);
  }
}

/* ---- orchestration ------------------------------------------------------ */

function drawAll(summary) {
  const { views, meta } = summary;
  const draw = (id, build) => {
    const container = document.getElementById(id);
    container.replaceChildren(build(container));
  };
  draw("chart-prevalence", (c) => renderPrevalence(views.prevalence_by_drug_class, c));
  drawLandscape(summary);
  draw("chart-levels", (c) => renderLevels(views.levels_by_drug, c));
  for (const gene of SDRM_GENES) {
    draw(`chart-sdrm-${gene.toLowerCase()}`, (c) => renderSdrm(views.sdrm_frequency, gene, c));
  }
  draw("chart-coverage", (c) => renderCoverage(views.gene_coverage, meta.n_accessions_scored, c));
}

function redrawOnResize(summary) {
  let lastWidth = document.documentElement.clientWidth;
  let timer = 0;
  new ResizeObserver(() => {
    clearTimeout(timer);
    timer = setTimeout(() => {
      const width = document.documentElement.clientWidth;
      if (width === lastWidth) return;
      lastWidth = width;
      drawAll(summary);
    }, 120);
  }).observe(document.documentElement);
}

async function main() {
  const expected = document.body.getAttribute("data-summary-schema");
  let summary;
  try {
    const response = await fetch(SUMMARY_URL);
    if (!response.ok) throw new Error(`HTTP ${response.status}`);
    summary = await response.json();
  } catch (error) {
    showBanner(`The resistance summary could not be loaded (${error.message}). Open this page from a web server, not from a file.`);
    return;
  }

  if (schemaMismatch(expected, summary?.meta?.schema_version)) {
    showBanner(
      `This page was built for summary schema ${expected}, but summary.json declares ${summary?.meta?.schema_version ?? "no version"}. No charts are drawn, because they could be wrong.`,
    );
    return;
  }

  try {
    applyBindings(summary);
    syncDigestNote(summary.meta);
    syncLandscapeAbsent(summary);
    renderLevelLegend();
    renderLandscapeLegend(
      { countries: summary.meta.landscape_countries, floor: summary.meta.landscape_min_cell_n },
      summary.views.resistance_landscape,
    );
    drawAll(summary);
  } catch (error) {
    showBanner(`The resistance summary does not have the shape this page expects (${error.message}). No charts are drawn.`);
    return;
  }
  redrawOnResize(summary);
}

main();
