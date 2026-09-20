"use strict";

const state = { groups: [], selected: null, detail: null, filter: "" };

async function api(path, options) {
  const res = await fetch(path, options);
  const text = await res.text();
  const value = text ? JSON.parse(text) : null;
  if (!res.ok) throw new Error((value && value.message) || res.statusText);
  return value;
}

function el(tag, attrs, children) {
  const node = document.createElement(tag);
  if (attrs) for (const [k, v] of Object.entries(attrs)) {
    if (k === "class") node.className = v;
    else if (k.startsWith("on")) node.addEventListener(k.slice(2), v);
    else node[k] = v;
  }
  (children || []).forEach((c) => {
    if (c == null) return;
    node.appendChild(typeof c === "string" ? document.createTextNode(c) : c);
  });
  return node;
}

async function loadState() {
  const build = document.getElementById("buildFilter").value || "";
  const data = await api("/api/state" + (build ? "?build=" + encodeURIComponent(build) : ""));
  document.getElementById("rulesVersion").textContent = data.rulesVersion;
  const select = document.getElementById("buildFilter");
  if (!select.dataset.filled) {
    data.builds.forEach((b) => select.appendChild(el("option", { value: b.name }, [b.name])));
    select.dataset.filled = "1";
  }
  state.groups = data.groups;
  renderGroups(data.groups);
  renderMismatches(data.mismatches);
  document.getElementById("groupCount").textContent = "(" + data.groups.length + ")";
}

function renderGroups(groups) {
  const ul = document.getElementById("groupList");
  ul.innerHTML = "";
  groups.forEach((g) => {
    const tags = [];
    if (!g.uniquePlacement) tags.push(el("span", { class: "badge repeat" }, ["重复 " + g.equivalentAnchors.length + " 位置"]));
    if (g.frozen) tags.push(el("span", { class: "badge frozen" }, ["已冻结"]));
    if (g.markedMismatch) tags.push(el("span", { class: "badge mismatch" }, ["参考不一致"]));
    const li = el("li", {
      class: state.selected === g.canonicalId ? "active" : "",
      onclick: () => selectGroup(g.canonicalId),
    }, [
      el("div", { class: "grow-head" }, [
        el("span", null, [g.contig + ":" + g.pos + " " + g.ref + ">" + g.alt]),
        el("span", { class: "small" }, [g.build]),
      ]),
      el("div", { class: "grow-sub" }, [g.canonicalId.split("|").slice(0, 1)[0] + " · " + g.rulesVersion]),
      el("div", null, tags),
    ]);
    ul.appendChild(li);
  });
}

async function selectGroup(id) {
  state.selected = id;
  renderGroups(state.groups);
  const d = await api("/api/groups/" + encodeURIComponent(id));
  state.detail = d;
  renderDetail(d);
}

function renderDetail(d) {
  document.getElementById("detailTitle").textContent =
    d.contig + ":" + d.pos + "  " + d.ref + " > " + d.alt;
  const box = document.getElementById("detail");
  box.innerHTML = "";

  box.appendChild(kv([
    ["规范 ID", d.canonicalId],
    ["参考构建", d.build],
    ["规则版本", d.rulesVersion],
    ["物理位置", d.uniquePlacement ? "唯一" : "重复序列，存在 " + d.equivalentAnchors.length + " 个等价锚点（首选最左 " + Math.min(...d.equivalentAnchors) + "）"],
    ["等价锚点", d.equivalentAnchors.join(", ")],
    ["审阅状态", (d.frozen ? "已冻结@" + d.frozenAnchor + " " : "") + (d.markedMismatch ? "已标记参考不一致 " : "") || "开放"],
  ]));

  if (d.alignmentWindow) box.appendChild(renderWindow(d.alignmentWindow, d));
  box.appendChild(renderNormSteps(d.normSteps));
  box.appendChild(renderEvidence(d.evidence));
  d.projections.forEach((p) => box.appendChild(renderProjection(p)));
  box.appendChild(renderReviews(d));
}

function kv(rows) {
  const grid = el("div", { class: "kv" });
  rows.forEach(([k, v]) => {
    grid.appendChild(el("b", null, [k]));
    grid.appendChild(el("span", null, [String(v)]));
  });
  return grid;
}

function renderWindow(w, d) {
  const wrap = el("div", null, [el("div", { class: "section-title" }, ["对齐窗口 (" + w.contig + ":" + w.windowStart + "-" + w.windowEnd + ")"])]);
  const line = el("div", { class: "window" });
  const seq = w.referenceWindow;
  const alleleStartRel = w.alleleStart - w.windowStart;
  const alleleEndRel = w.alleleEnd - w.windowStart;
  for (let i = 0; i < seq.length; i++) {
    const gPos = w.windowStart + i;
    let cls = "";
    if (gPos >= w.alleleStart && gPos <= w.alleleEnd) cls = "anchor";
    line.appendChild(el("span", { class: cls, title: String(gPos) }, [seq[i]]));
  }
  wrap.appendChild(line);
  const altLine = el("div", { class: "window alt" }, ["ALT 锚点后插入/替换：" + d.alt]);
  wrap.appendChild(altLine);
  return wrap;
}

function renderNormSteps(steps) {
  const wrap = el("div", null, [el("div", { class: "section-title" }, ["规范化计算路径"])]);
  if (!steps || !steps.length) { wrap.appendChild(el("div", { class: "small" }, ["记录本身即为最小左对齐形式"])); return wrap; }
  const table = el("table", null, [
    el("thead", null, [el("tr", null, ["步骤", "POS", "REF", "ALT", "说明"].map((h) => el("th", null, [h])))]),
    el("tbody", null, steps.map((s) => el("tr", null, [
      el("td", null, [s.op]), el("td", { class: "mono" }, [String(s.pos)]),
      el("td", { class: "mono" }, [s.ref]), el("td", { class: "mono" }, [s.alt]),
      el("td", { class: "small" }, [s.note]),
    ]))),
  ]);
  wrap.appendChild(table);
  return wrap;
}

function renderEvidence(ev) {
  const wrap = el("div", null, [el("div", { class: "section-title" }, ["来源证据（复合等位拆分双向可追溯）"])]);
  const table = el("table", null, [
    el("thead", null, [el("tr", null,
      ["批次", "行", "contig:POS", "REF>ALT", "来源ID", "角色", "REF校验", "规范组"].map((h) => el("th", null, [h])))]),
    el("tbody", null, ev.map((o) => el("tr", null, [
      el("td", { class: "mono" }, [o.batchId]),
      el("td", null, [String(o.lineNumber)]),
      el("td", { class: "mono" }, [o.contig + ":" + o.pos]),
      el("td", { class: "mono" }, [o.ref + ">" + o.alt]),
      el("td", null, [o.sourceId]),
      el("td", null, [o.altIndex === 0
        ? "原始记录" + (o.childCount ? "（拆出 " + o.childCount + " 条）" : "")
        : "拆分子记录(父行 " + o.parentLineNumber + ", ALT#" + o.altIndex + ")"]),
      el("td", { class: o.refStatus === "MATCH" ? "ok" : "fail" }, [o.refStatus]),
      el("td", { class: "mono small" }, [o.canonicalId ? "→ " + o.canonicalId.split("|").slice(2).join(":") : "—"]),
    ]))),
  ]);
  wrap.appendChild(table);
  return wrap;
}

function renderProjection(p) {
  const wrap = el("div", null, [
    el("div", { class: "section-title" }, ["转录本投影：" + p.transcriptId + " (" + p.build + ", 规则 " + p.rulesVersion + ")"]),
  ]);
  if (p.failure) {
    wrap.appendChild(el("div", { class: "fail" }, ["失败类型：" + p.failure]));
  } else {
    wrap.appendChild(el("div", null, [
      el("span", { class: "tag unique" }, [p.consequence]),
      p.cdsStart ? el("span", { class: "small" }, ["  CDS " + p.cdsStart + (p.cdsEnd && p.cdsEnd !== p.cdsStart ? "-" + p.cdsEnd : "") +
        (p.proteinStart ? " / 蛋白 " + p.proteinStart : "")]) : null,
    ]));
    if (p.refCodon) wrap.appendChild(el("div", { class: "mono" }, [
      p.refCodon + " (" + p.refAa + ") → " + p.altCodon + " (" + p.altAa + ")",
    ]));
    if (p.affectedExons && p.affectedExons.length) {
      wrap.appendChild(el("div", { class: "small" }, ["命中外显子：" + p.affectedExons.join(", ")]));
    }
  }
  const pathList = el("div", { class: "path" }, [(p.path || []).map((s, i) => (i + 1) + ". " + s).join("\n")]);
  wrap.appendChild(pathList);
  return wrap;
}

function renderReviews(d) {
  const wrap = el("div", null, [el("div", { class: "section-title" }, ["审阅（乐观版本锁）"])]);
  const version = d.reviews.length ? d.reviews[d.reviews.length - 1].version : 0;
  const reviewerInput = el("input", { value: "reviewer@" + Math.floor(Math.random() * 1000), placeholder: "审阅者" });
  const anchorInput = el("input", { placeholder: "冻结锚点(可空=首选)", style: "width:150px" });
  const noteInput = el("input", { placeholder: "备注", style: "flex:1" });
  const bar = el("div", { class: "review-bar" }, [
    reviewerInput, anchorInput, noteInput,
    el("button", { onclick: () => submitReview(d, "FREEZE", version, reviewerInput.value, anchorInput.value, noteInput.value) }, ["冻结规范化结果"]),
    el("button", { onclick: () => submitReview(d, "MARK_MISMATCH", version, reviewerInput.value, null, noteInput.value) }, ["标记参考不一致"]),
    el("button", { onclick: () => submitReview(d, "CLEAR", version, reviewerInput.value, null, noteInput.value) }, ["清除标记"]),
  ]);
  wrap.appendChild(bar);
  wrap.appendChild(el("div", { class: "small" }, ["当前版本：" + version + "（提交时必须携带此版本号，过期提交返回 409）"]));
  const table = el("table", null, [
    el("thead", null, [el("tr", null, ["版本", "审阅者", "动作", "锚点", "备注", "时间"].map((h) => el("th", null, [h])))]),
    el("tbody", null, d.reviews.map((r) => el("tr", null, [
      el("td", null, [String(r.version)]), el("td", null, [r.reviewer]), el("td", null, [r.action]),
      el("td", null, [r.anchorOverride == null ? "—" : String(r.anchorOverride)]),
      el("td", null, [r.note || ""]), el("td", { class: "small" }, [r.submittedAt]),
    ]))),
  ]);
  wrap.appendChild(table);
  return wrap;
}

async function submitReview(d, action, expectedVersion, reviewer, anchor, note) {
  try {
    await api("/api/reviews", {
      method: "POST", headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        canonicalId: d.canonicalId, action, expectedVersion, reviewer,
        anchorOverride: anchor ? Number(anchor) : null, note: note || null,
      }),
    });
    await loadState();
    await selectGroup(d.canonicalId);
  } catch (e) { alert("提交被拒绝：" + e.message); }
}

function renderMismatches(rows) {
  const ul = document.getElementById("mismatchList");
  ul.innerHTML = "";
  if (!rows.length) { ul.appendChild(el("li", { class: "small" }, ["无"])); return; }
  rows.forEach((r) => ul.appendChild(el("li", null, [
    el("div", { class: "mono" }, [r.build + " " + r.contig + ":" + r.pos + " " + r.ref + ">" + r.alt]),
    el("div", { class: "small" }, [r.batchId + " 行" + r.lineNumber + " (" + r.sourceId + ")"],
  )])));
}

document.getElementById("refreshBtn").addEventListener("click", loadState);
document.getElementById("buildFilter").addEventListener("change", loadState);

document.getElementById("importBtn").addEventListener("click", async () => {
  const out = document.getElementById("importResult");
  try {
    const r = await api("/api/batches", {
      method: "POST", headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        batchId: document.getElementById("batchId").value,
        build: document.getElementById("batchBuild").value,
        vcf: document.getElementById("batchVcf").value,
      }),
    });
    out.textContent = JSON.stringify(r, null, 2);
    await loadState();
  } catch (e) { out.textContent = "ERROR: " + e.message; }
});

document.getElementById("loBtn").addEventListener("click", async () => {
  const out = document.getElementById("loResult");
  try {
    const r = await api("/api/liftover", {
      method: "POST", headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        sourceBuild: document.getElementById("loSrc").value,
        destBuild: document.getElementById("loDst").value,
        chainTsv: document.getElementById("loChain").value,
      }),
    });
    out.textContent = "汇总：" + JSON.stringify(r.summary) + "\n\n" +
      r.items.map((i) => i.status.padEnd(14) + " " + i.canonicalId + " -> " + (i.destContig || "-") +
        ":" + (i.destPos || "-") + "  " + i.note).join("\n");
  } catch (e) { out.textContent = "ERROR: " + e.message; }
});

loadState().catch((e) => alert("初始化失败：" + e.message));
