# 变体规范化与转录本投影审阅台（离线）

一个完全离线的服务：导入精简 VCF 与本地参考序列 fixture，将同一生物学变化的
多种写法（左对齐差异、多碱基分解、参考构建差异）归并到一个**稳定规范记录**，
再把每个等位基因投影到多个转录本。不做临床解读，不联网查询任何数据库；
但严格保留来源批次、参考构建和注释规则版本。

- Kotlin 服务（JDK 内置 HTTP Server，无 Web 框架）
- SQLite 持久状态（单文件，WAL）
- 浏览器单页审阅界面（原生 HTML/CSS/JS，无构建步骤）
- 固定、无外部依赖的 FASTA / TSV / VCF 测试数据

## 构建与演示

```bash
# 安装（跳过测试打包）
mvn -q -DskipTests package

# 演示：先跑测试，再在 5318 端口启动服务并自动载入固定 fixture
mvn -q test && mvn -q exec:java -Dexec.mainClass=app.MainKt -Dexec.args='--port 5318'
```

界面地址：<http://127.0.0.1:5318>

首次启动会自动载入 `src/main/resources/fixtures` 下的两个参考构建（b37/b38）、
两个 b37 转录本、一个种子批次 `SEED-B37-0001` 以及一条 b37→b38 映射链。
加 `--reset` 可清空数据库重新播种；加 `--no-seed` 启动空库；`--db <path>`
指定 SQLite 文件（默认 `target/variant-normalizer.db`）。

## VCF 坐标约定

遵循 VCF 4.2：

- `POS` 为 **1-based**；`REF`/`ALT` 至少各含一个碱基，包含一个共享锚点碱基。
- 替换：`len(REF)=len(ALT)=1`，变更碱基就是 `POS`。
- 删除：`len(REF)>len(ALT)`，被删参考碱基占据 `POS+1 .. POS+len(REF)-1`，
  锚点碱基位于 `POS`。
- 插入：`len(ALT)>len(REF)`，插入的 `len(ALT)-1` 个碱基附着在锚点 `POS` 之后。
- 逗号分隔的多等位 `ALT` 在导入时**拆分**为子记录；原始复合行原样保留。
- 仅支持核苷酸等位；不支持符号等位（`<...>`）与 `*`。

解析阶段只做结构校验。即使 `REF` 与参考不符（例如写成了错误的碱基），记录也会
进入管线并被标记 `MISMATCH`，而不是在解析时被丢弃——这样参考不一致才能被
审阅者看到并标记。

## 规范化（最小表示 + 左对齐）

管线顺序：

1. 去除共享右侧锚点（保留至少一个碱基）；
2. 去除共享左侧前缀；
3. 在重复序列中做最简约左对齐（vt 风格的等位滑动）；
4. 再次最小化。

每一步都以 `TRIM_RIGHT / TRIM_LEFT / LEFT_ALIGN` 记录在规范组的“计算路径”里，
界面可逐步审计。

### 等价关系

在串联重复区域，同一分子可以有多个合法 VCF 位置。系统会：

- 计算一个**首选规范表示**：可左对齐时一律取最左锚点；
- 枚举所有重新规范化后回到同一记录的**物理锚点位置**
  （`equivalentAnchors`）；
- 当存在多个等价位置时 `uniquePlacement=false`，界面标记“重复序列多位置”，
  **绝不假装只有一个物理位置**；
- 对这类未唯一左对齐的等位，转录本投影会被显式扣留
  （`NON_UNIQUE_LEFT_ALIGN`）。

注意单碱基删除与插入的区别：在纯合 run（如 8 个 A）中，

- 删除**具体某个** A 与删除其他 A 可能是不同分子，会落到不同规范组；
- 插入一个 A，无论写在 run 中哪个锚点，都是同一分子，因此该插入的
  `equivalentAnchors` 覆盖整段 run。

稳定规范 ID：`build|contig|pos|REF|ALT`。

### 多等位拆分与双向追溯

- 每个源 VCF 行只保留**一条原始证据行**（`altIndex=0`，记录拆出的子记录数）。
- 每个 ALT 产生一条子记录（`altIndex=1..n`），带 `parentLineNumber` 指回原行。
- 两条子记录即便各自独立“重放”，也只对应同一条原始证据，不会重复计数。

## 转录本投影

转录本以 TSV 提供，坐标全部 1-based、闭区间，并显式给出：

- 外显子 `exon_starts/exon_ends`；
- 编码外显子片段 `cds_exon_starts/cds_exon_ends`（用于精确处理 UTR 与
  外显子边界）；
- 链向 `+ / -`。负链转录本的 CDS 第 1 位是基因组上最高的编码碱基，密码子
  按反向互补拼接。

投影会输出后果、CDS/蛋白坐标、参考/替代密码子与氨基酸、命中外显子，以及
可读的逐步计算路径。后果类型包括：

`START_LOST`、`STOP_GAINED`、`STOP_LOST`、`MISSENSE_VARIANT`、
`SYNONYMOUS_VARIANT`、`FRAMESHIFT_VARIANT`、`INFRAME_INDEL`、
`INTRON_VARIANT`、`FIVE_PRIME_UTR_VARIANT`、`THREE_PRIME_UTR_VARIANT`、
`PROJECTION_FAILED`。

### 投影失败类型（`failure` 字段）

| 代码 | 含义 |
| --- | --- |
| `CONTIG_NOT_FOUND` | 转录本所在 contig 不在当前参考构建中 |
| `OUTSIDE_TRANSCRIPT` | 等位区间在转录本首末外显子范围之外（基因间） |
| `SPANS_INTRON` | 缺失碱基跨越外显子/内含子边界，或插入锚点正好落在外显子边界，编码影响不唯一 |
| `REFERENCE_MISMATCH` | `REF` 与当前参考构建不符，投影被扣留 |
| `NON_UNIQUE_LEFT_ALIGN` | 等位在重复序列中有多个等价物理位置，无法唯一投影 |
| `UNKNOWN_ALLELE` | 等位含歧义碱基（非 ACGT），无法可靠翻译 |

坐标转换覆盖：染色体起点、插入、删除、跨外显子边界密码子、负链转录本。

## 参考构建映射（liftover）

映射文件为 TSV 链（列：`src_contig src_start1 src_end1 dst_contig dst_start1 dst_end1 dst_strand`）。
导入是**原子操作**，并先生成覆盖源构建全部规范组的影响清单：

- `MAPPED`：唯一、同链映射（带位移说明）；
- `STRAND_FLIPPED`：唯一但目标链方向翻转，等位做反向互补；
- `GAP`：等位区间跨越相邻链段之间的组装缺口，或目标碱基与提升后的 REF 不符；
- `COLLISION`：源区间被多条链的源范围**重叠**覆盖（多对一），提升不唯一；
- `OUT_OF_MAP`：没有任何链段覆盖该源区间。

相邻（仅相接、不重叠）链段构成的边界断裂判为 `GAP` 而非 `COLLISION`；
只有源范围真正重叠才是多对一 `COLLISION`。

老结果仍按原参考构建复现：每条投影都记录 `build` 与 `rulesVersion`，
导入新链不会改写或删除 b37 的规范组与投影。注释规则语义变化时提升
`ANNOTATION_RULES_VERSION`，历史行保留其产生时的版本。

## 审阅与并发

- 审阅动作：`FREEZE`（冻结某个规范化锚点）、`MARK_MISMATCH`（标记参考不一致）、
  `CLEAR`（清除标记）、`COMMENT`。
- 每个规范组维护单调递增的 `version`；提交必须携带 `expectedVersion`，
  过期提交返回 **409**，不会覆盖别人已提交的新版本。历史审阅全部保留。
- 冻结后在规范组上记录 `frozen/frozen_anchor`；标记参考不一致记录
  `marked_mismatch`。

## 原子性与幂等

- 批次导入、参考/转录本载入、映射链导入都包在单个 SQLite 事务里：全部成功
  提交，否则整体回滚（失败批次不留任何组、证据或收据）。
- 相同 `batchId` + 相同负载重试：返回**原始收据**（`duplicate=true`）。
- 相同 `batchId` + 不同负载：返回 **409 CONFLICT**，原收据保留不变。

## HTTP 接口

| 方法 | 路径 | 说明 |
| --- | --- | --- |
| GET | `/api/state?build=` | 构建、等价组、参考不一致总览 |
| GET | `/api/groups/{id}` | 规范组详情：证据、规范化步骤、对齐窗口、投影、审阅 |
| GET | `/api/transcripts/{build}` | 转录本外显子/编码结构 |
| GET | `/api/receipts/{batchId}` | 批次原始收据 |
| POST | `/api/batches` | 原子导入 VCF（幂等） |
| POST | `/api/reviews` | 提交审阅（乐观版本锁） |
| POST | `/api/liftover` | 原子导入映射链并返回影响清单 |

## 固定 fixture

位于 `src/main/resources/fixtures`：

- `ref_b37.fa` / `ref_b38.fa`：三个 contig（chr11 正链基因、chr17 负链基因、
  chr20 含染色体起点与 poly-A/poly-T 重复）；b38 演示 +10 位移、链翻转、
  10 碱基 gap 与多对一链；
- `transcripts_b37.tsv` / `transcripts_b38.tsv`；
- `seed_b37.vcf`：覆盖起点、多等位、左对齐、参考不符、起始/终止密码子、
  跨外显子、内含子、UTR、基因间、负链、唯一 in-frame 删除及各类 liftover；
- `map_b37_to_b38.tsv`；
- `expected_norm.json` / `expected_proj.json`：由 `tools/gen_fixtures.py`
  生成的规范化/投影期望值。

`tools/gen_fixtures.py` 仅用于**生成与核对**固定数据；运行时服务不依赖 Python。
