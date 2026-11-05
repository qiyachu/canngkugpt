# canngkugpt

## RASC学习曲线分析

按48例单术者数据生成手术时间趋势、传统CUSUM、病例构成调整CUSUM及敏感性分析图。
原始临床数据应保存在本地，不提交到GitHub。此仓库只保存代码与依赖声明。

### 安装与运行

需要Python 3.12。

```bash
python -m venv .venv --system-site-packages
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python analysis.py /absolute/path/to/input.xlsx outputs
```

使用现有隔离环境和检出目录，无需另建Git worktree。Excel需包含Sheet1，48例按实际手术顺序排列，
列名与分析脚本一致。脚本检查队列结构、必需字段和研究基准值，因此不适合直接应用于其他队列。

### 分析约定

- 采用总手术时间_min列，调整BMI、同期子宫切除、后壁或会阴修补、抗尿失禁手术。
- 两组各10,000次置换，固定随机种子，使用最大绝对累计偏差和加一校正P值。
- 灰带是逐点95%置换参考区间；残差置换的交换性是假设，结果应审慎解释。
- 敏感性分析包括逐例留一、POP-Q分期调整及完整5例居中移动平均后重新中心化。
- HC3回归主结果采用t分布；另提供正态近似结果，论文应明确推断约定。
- 同日病例保留Excel行顺序，发表前需核实真实手术先后。

### 输出

outputs内生成3幅图（PNG、TIFF、PDF、SVG）、结果Excel、中文核查报告、英文图注及summary.json。
输出文件包含研究结果，默认被Git忽略。

### 验证

```bash
.venv/bin/python -m pip check
.venv/bin/python analysis.py --help
```

功能验证需要用户提供的本地工作簿；执行完整命令后应生成3幅图、48行Figure_data，
两条主CUSUM终点应接近0。仅运行帮助命令不代表分析已验证。
