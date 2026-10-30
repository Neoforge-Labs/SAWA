/* SAWA resistance summary: draws summary.json as inline SVG.
 *
 * No framework and no network beyond the one same-origin fetch. The order
 * constants below are copies of paper1/style.py and dashboard/export/aggregate.py
 * (the page cannot import Python); tests/dashboard/test_static.py keeps them in step.
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
  return meta.input_file_md5 !== meta.input_dvc_md5;
}

function syncDigestNote(meta) {
  document.getElementById("digest-note").hidden = !digestsDiffer(meta);
}

/* ---- text bindings ------------------------------------------------------ */

function deriveFigures(summary) {
  const { meta, views } = summary;
  return {
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

function resolveBinding(path, summary, derived) {
  const [scope, key] = path.split(".");
  let value;
  if (scope === "meta") value = summary.meta[key];
  else if (scope === "derived") value = derived[key];
  else if (scope === "gene") value = summary.views.gene_coverage.find((row) => row.gene === key)?.n_sequences;
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
    renderLevelLegend();
    drawAll(summary);
  } catch (error) {
    showBanner(`The resistance summary does not have the shape this page expects (${error.message}). No charts are drawn.`);
    return;
  }
  redrawOnResize(summary);
}

main();
