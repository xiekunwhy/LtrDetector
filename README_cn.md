# LtrDetector 修改版说明

[English README](README.md)

本仓库基于 [BioinformaticsToolsmith/LtrDetector](https://github.com/BioinformaticsToolsmith/LtrDetector)，保留原项目提交历史、检测算法和许可说明，主要改进输入方式、输出文件安全性及 LTR_retriever 格式转换。

## 主要修改

### 1. 直接接收 FASTA 文件

新增 `-seqFile`，可直接扫描一个 FASTA 文件，也支持文件内包含多条序列，不再必须先建立输入目录。

以下三种输入方式互斥，必须且只能指定一种：

- `-seqFile <FASTA文件>`：直接扫描一个 FASTA 文件。
- `-fasta <目录>`：保留原代码的目录输入方式，扫描目录下的 `.fa` 文件。
- `-chromDir <目录>`：作为 `-fasta` 的兼容别名，与原 README 中的用法一致。

### 2. 自定义输出前缀

新增 `-prefix <前缀>`。例如使用 `-prefix sample1` 并开启相应输出选项时，生成：

```text
sample1.Detector.bed
sample1.RawScores.csv
sample1.CleanScores.csv
sample1.NestedDetector.bed
```

不指定 `-prefix` 时，仍保留原来的按输入 FASTA 文件名命名的行为。

### 3. 消除输出目录误删风险

原代码启动时会删除 `-destDir` 下所有以 `Detector.bed` 结尾的文件，即使它们属于其他输入或其他前缀。

本版本已取消该启动清理逻辑：

- 启动时不再删除输出目录中的旧文件。
- 不同 `-prefix` 的结果可以安全写入同一个目录。
- 普通文件和其他前缀的 BED、CSV 均不受影响。
- 同名输出在本次运行第一次写入时覆盖，避免同一前缀重跑时重复追加旧结果。
- 同一 FASTA 中的后续序列在本次运行中继续追加，避免只保留最后一条序列的 BED 结果。

建议共享输出目录时为每个输入使用唯一前缀；同名文件仍会被覆盖。

### 4. 明确的 bug 修复

- 修复可选参数不存在时仍递增尾迭代器的未定义行为。
- 增加 Windows 路径分隔符和 MinGW filesystem 路径转换支持。
- 增加 FASTA 文件无法打开、空文件、缺少标题及标题前有序列内容时的检查。
- 跳过长度小于 `k` 的区段，避免 k-mer 评分越界。
- 修复 k-mer 动态数组误用 `delete` 而非 `delete[]` 的释放错误。
- 修复奇数长度区间中位数取值错误。
- 修复 LTR 边界扩展分支中的赋值及扩展长度错误。
- 修复边界过滤后只剩一个候选时，该候选被去重流程丢弃的问题。

本次检查及限时测试中未发现主流程可复现的死循环。较长的运行时间主要来自局部比对的二次方复杂度，其规模受 `-maxLenLTR` 等参数影响；为保持最小改动，没有重写检测算法。

## Linux 静态编译

需要 GCC、Make、glibc 静态库及 OpenMP 静态运行库：

```bash
cd src
make bin  # 仅首次构建、目标目录尚不存在时需要
make clean
make tr -j4 CXXFLAGS="-std=c++11 -O3 -fopenmp -fmessage-length=0 -static"
mv ../bin/LtrDetector ../bin/LtrDetector-static
```

检查结果：

```bash
file ../bin/LtrDetector-static
ldd ../bin/LtrDetector-static
```

静态链接成功时，`ldd` 应显示 `not a dynamic executable`。

原 Makefile 的少量增量依赖不完整，修改源码后建议先运行 `make clean`，避免复用旧目标文件。

## 使用示例

```bash
./LtrDetector-static \
  -seqFile genome.fa \
  -destDir results \
  -prefix Chr01 \
  -nThreads 4
```

默认只生成主结果：

```text
results/Chr01.Detector.bed
```

其他输出需要显式开启：

```bash
-rawScores       # Chr01.RawScores.csv
-cleanedScores   # Chr01.CleanScores.csv
-nested          # Chr01.NestedDetector.bed
```

多个输入可以共享输出目录：

```bash
./LtrDetector-static -seqFile Chr01.fa -destDir results -prefix Chr01
./LtrDetector-static -seqFile Chr02.fa -destDir results -prefix Chr02
```

注意：原程序的 RawScores 和 CleanScores 在同一输入包含多条序列时会被后续序列覆盖，最终 CSV 只保留最后一条序列的分数；本版本未改变这一行为。

## 转换为 LTR_retriever 输入

仓库包含两个 Perl 转换脚本，均要求 LtrDetector 的默认 18 列输出，**不要使用 `-bedFormat`**。

### LTRharvest 风格

`convert_ltrdetector.pl` 修正了官方脚本与当前 18 列输出之间的 TSD、链方向字段映射问题：

```bash
perl convert_ltrdetector.pl genome.fa Chr01.Detector.bed > Chr01.harvest.scn

LTR_retriever \
  -genome genome.fa \
  -inharvest Chr01.harvest.scn \
  -threads 4
```

### LTR_FINDER 风格

`convert_ltrdetector_to_finder.pl` 输出 12 列 LTR_FINDER 风格文件，并将 identity 百分数转换为 0–1 的 similarity：

```bash
perl convert_ltrdetector_to_finder.pl genome.fa Chr01.Detector.bed > Chr01.finder.scn

LTR_retriever \
  -genome genome.fa \
  -infinder Chr01.finder.scn \
  -threads 4
```

转换和运行 LTR_retriever 时，应使用与 LtrDetector 结果坐标对应的同一个 FASTA。若输入是全基因组分块，合并前需要统一序列 ID 并换算成全染色体坐标。

格式转换不会改变候选或坐标，也不会使 LtrDetector 结果变成实际的 LTR_FINDER/LTRharvest 检测结果。

## 原始参数与许可

完整参数、输出列说明及原项目编译说明见 [英文 README](README.md)。原项目的学术/非盈利及商业使用许可说明保持不变。
