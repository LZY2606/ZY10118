# Pair-wise GSB — 离线变异规范化与多转录本投影审阅台

一个**无联网、无外部数据库依赖**的离线应用：导入精简 VCF 与本地参考序列 fixture，
把同一生物学变化的不同写法（左对齐差异、多碱基拆分、表述差异）归入一个稳定的
**等价组**，再投影到一个或多个转录本。参考构建之间的坐标迁移（liftover）以本地
chain fixture 完成，并在导入时先生成**影响清单**。

- 语言/运行时：Kotlin（JVM 17+ 字节码）、JDK 内置 HTTP Server、SQLite（`sqlite-jdbc`，自带本地库）
- 存储：单个 SQLite 文件（默认 `.gsb-data/state.sqlite`）
- 界面：原生 HTML/CSS/JS 单页，无前端框架、无 CDN
- 规则版本固定在代码中：`annot-rules-1.0.0`（规范化与投影规则）、`tx-fixture-1`（转录本集）

## 构建与演示

```bash
# 安装/打包（跳过测试）
mvn -q -DskipTests package

# 演示：跑全部测试并启动服务
mvn -q test && mvn -q exec:java -Dexec.mainClass=app.MainKt -Dexec.args='--port 5318'
```

界面地址：<http://127.0.0.1:5318>

可选参数：

- `--port 5318`：监听端口（默认 8080，仅绑定 127.0.0.1）
- `--data 路径`：SQLite 与参考文件目录（默认 `.gsb-data`）
- `--seed`：启动时导入内置示例批次 `DEMO-BATCH-0001` 与示例 chain（已导入则仅返回原结果）

也可直接用界面上的「填入示例 VCF」「填入示例 chain」按钮加载固定测试数据。

## 目录

```
src/main/kotlin/app/
  Models.kt             领域模型与版本常量
  FastaRef.kt           参考序列访问（1-based 闭区间）
  Normalizer.kt         REF 校验、最小表示、左对齐、重复区等价位置枚举
  VcfParser.kt          精简 VCF 解析（支持多等位基因拆分）
  TranscriptProjector.kt 转录本坐标/链转换、区域分类、失败类型
  Liftover.kt           chain 解析与跨构建映射（5 种状态）
  AppService.kt         批次/等价组/证据/审阅/liftover 的事务化用例
  Db.kt                 SQLite schema 与事务封装
  WebServer.kt / Json.kt / Main.kt / Resources.kt
src/main/resources/
  fixtures/refs/        hg19_chr1.fasta、hg38_chr1.fasta（确定性生成的迷你参考）
  fixtures/transcripts/ transcripts.tsv
  fixtures/chains/      hg19_to_hg38.chain.tsv
  fixtures/sample_batch.vcf
  web/index.html        审阅界面
tools/gen_fixtures.py   重新生成全部 fixture（确定性种子，可复现）
src/test/kotlin/app/    24 个 JUnit 5 测试（无网络）
```

## VCF 坐标约定

- 全部使用 **VCF 标准的 1-based 闭区间**；`POS` 是 `REF` 第一个碱基在参考上的位置。
- `REF`/`ALT` 的第一个碱基是**锚定碱基**：
  - 纯插入 n 个碱基：`REF` 长度 1、`ALT` 长度 n+1（锚定 + 插入序列）。
  - 纯删除 n 个碱基：`REF` 长度 n+1、`ALT` 长度 1。
  - MNV/SNV 照常写等长或不等长等位。
- 处理顺序固定为：**先校验 REF，再最小表示，再左对齐**。
  - REF 校验：`POS..POS+len(REF)-1` 必须落在当前活动参考构建上且碱基完全一致；
    不一致或越界的记录**不参与规范化/投影**，单独成组并保留“参考期望碱基”，
    可由审阅者标记 `REF_FLAGGED`。
  - 最小表示：先去掉共同后缀、再去掉共同前缀，任一等位至少保留 1 个碱基。
  - 左对齐：逐次借取紧邻上游参考碱基并丢弃共同后缀，直到无法继续（bcftools/rtg 语义）。
- 只接受 DNA 字符（`ACGTN`，大小写归一）；不接受 gVCF 的 `ALT=.`。
- 多等位基因行（`ALT=A,T`）视为**复合等位基因**：原行与每个拆分等位都落库并可双向追溯。

## 等价关系

- 一个**等价组**由 `(build, chrom, 规范化 POS, 规范化 REF, 规范化 ALT)` 唯一确定。
- 规范化键与构建绑定：同一变异在不同构建上是不同记录；老构建结果不会被新构建覆盖。
- **简单重复中的物理位置歧义被显式保留**：对于 homopolymer、二核苷酸/三核苷酸串联重复，
  同一插入/删除存在多个同样合法的最小 VCF 表述。系统枚举**全部**等价 `Placement`，
  用「最小表示后最靠左」作为首选（canonical），并在结果中记录
  `preferenceRule = leftmost placement after minimal representation …`。
  界面等价组表格显示等价位置数量与完整列表，绝不假装只有一个物理位置。
- 来自不同批次/管线的同一生物学变化合并到同一组；每条来源证据都保留
  （批次、来源 `SRC`、行号、原始行、ALT 序号）。
- 复合行的证据角色：
  - `ORIGINAL`：原始复合记录（每个组对该 raw record 至多 1 行，部分唯一索引保证
    “两条拆分记录分别重放不会产生两份原始证据”）。
  - `DERIVED`：该组对应的具体拆分等位（带 `altIndex`）。

## 转录本投影

转录本 fixture（TSV）：

```
txId	gene	chrom	strand	exons	cdsStart	cdsEnd	build
NM_DEMO_P	DEMO1	chr1	+	201-289,370-459	209	456	hg19
NM_DEMO_N	DEMO2	chr1	-	1501-1580,1661-1740	1512	1732	hg19
```

- 坐标转换覆盖：染色体起始（SNV）、插入、删除、**跨外显子边界**、**负链转录本**。
- 正链：基因组坐标按外显子长度累加得到 cDNA 坐标，等位按原样使用。
- 负链：cDNA 坐标反向计数，`REF`/`ALT` 做**反向互补**，计算路径中逐步展示。
- 每个后果都带一条**计算路径**（外显子索引、受影响区间、区域分类、链处理、cDNA 等位、HGVS 风格记法）。

### 投影失败类型（README 约定，状态码同时入库与返回）

| 状态码 | 触发条件 |
| --- | --- |
| `WRONG_BUILD` | 转录本标注的构建与变异所在构建不同（老结果仍按原构建可复现） |
| `UNKNOWN_CHROM` | 变异 contig 与转录本 contig 不同 |
| `OUT_OF_TRANSCRIPT` | REF 区间完全在转录本起始终止范围之外 |
| `EXON_STRADDLE` | REF 区间跨越外显子—内含子边界：同时命中两个外显子且跨过中间内含子，或与单个外显子部分重叠并延伸进相邻内含子；cDNA 层面不是一次连续替换 |
| `CDS_INCONSISTENT` | 转录本 CDS 注释与外显子结构导出的 cDNA CDS 区间矛盾 |

非失败但不产生连续 cDNA 等位的情形：`INTRONIC`（深内含子）、`SPLICE_REGION`
（距外显子/内含子边界 ≤3 bp）、`UTR5`、`UTR3`。

## 跨构建映射（liftover）

chain fixture 为 9 列 TSV（全部 1-based 闭区间）：

```
fromBuild toBuild chrom srcStart srcEnd strand dstChrom dstStart dstEnd
```

每个等价组的 canonical placement 在导入 chain 后计算映射，状态为：

| 状态 | 含义 |
| --- | --- |
| `UNIQUE` | 单段、正链、目标坐标唯一 |
| `FLIPPED` | 命中的段方向为 `-`；坐标反向、等位反向互补，界面明确标记方向翻转 |
| `MANY_TO_ONE` | 多个源段把同一区间/坐标映射到不同目标（方向冲突或目标区间非包含式重叠），坐标不唯一 |
| `INTERRUPTED` | chain 在等位覆盖区间内断裂：段接触了区间但无法完整跨越（缺口/在段边界中断且无续接） |
| `UNMAPPED` | 附近没有任何段覆盖该坐标 |

- 导入 chain 是**原子事务**：先解析校验（源/目标长度必须一致），删除旧段与旧结果、
  插入新段、重算所有等价组映射并返回 `ImportImpact`（UNIQUE/FLIPPED/MANY_TO_ONE/
  INTERRUPTED/UNMAPPED 计数与告警）。
- 链段只决定“向目标构建的投影”；**hg19 上的等价组、后果与证据永不被改写**，
  老版本结果按原参考可复现。

## 审阅与并发

- 审阅状态：`NONE` / `FROZEN`（冻结某个规范化结果）/ `REF_FLAGGED`（标记参考不一致）。
- 提交采用**乐观锁**：请求携带所基于的 `expectedVersion`；他人已提交新版本时
  返回 400 并要求重新加载，**新审阅不能覆盖别人已提交的版本**。
- 批次导入幂等：相同 `batchId` 重试只返回**原收据**（记录数、等位组数、时间戳一致，
  `replay=true`），不新增 raw 记录与证据。同一 `batchId` 换到不同构建重放会被拒绝。
- 批次导入与参考映射导入都包在单事务中，中途任何异常整体回滚。

## HTTP API（均返回 JSON）

- `GET /api/health` 规则/转录本版本、已加载构建
- `GET /api/transcripts`、`GET /api/batches`
- `GET /api/groups?build=hg19`、`GET /api/groups/{id}`（含 placements、normSteps、
  evidence、consequences[].path、liftovers、reviewVersion）
- `POST /api/import-batch` `{batchId, build, vcf}`
- `POST /api/reviews` `{groupId, expectedVersion, reviewer, state, note}`
- `POST /api/import-chain` `{chainTsv}` → 影响清单
- `GET /api/sample-vcf`、`GET /api/sample-chain`
- `GET /` 审阅界面

## 重新生成固定数据

```bash
python3 tools/gen_fixtures.py
```

生成器使用固定随机种子写入迷你参考，并把关键 motif 放到确定坐标；样例 VCF 的
`REF` 始终从生成后的参考中读取，避免手抄漂移。生成结果已随仓库提交。
