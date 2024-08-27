# 1. 原始数据的质控与裁剪

DamID-seq的原始数据是fastq格式，与其他二代测序数据相同，都需要事先进行质控和接头/低质量碱基裁剪后才能放心使用。在DCL data server上，我们使用fastqc和trim_galore软件进行操作。

```         
conda activate rnaseq_p3 #首先激活conda环境
# 质控
fastqc -q -t 4 -o <Output dir> <x.fq.gz y.fq.gz ...> 
# 裁剪
trim_galore -q 20 --phred33 --length 20 -e 0.1 --gzip -o ${path}/trimmed \
--paired <双端测序的第一个文件> <双端测序的第二个文件>
```

得到的fastq格式文件需要进行改名，以适应下一步的输入要求。**例如，对于对双端测序输入Example_R1.fq.gz Example_R2.fq.gz：**

1.  trim_galore得到的默认输出文件名为**Example_R1_val_1.fq.gz Example_R2_val_2.fq.gz**
2.  **将其分别改名为“Example_trimmed_1.fq.gz Example_trimmed_2.fq.gz”。**重点是两个文件除了末尾的1、2不同外，其他都相同。

# 2. 使用damidseq_pipeline软件进行上游分析

## (Optional) 2.1 参数传入

damidseq_pipeline软件可以储存传入过的index路径(bowtie2 index和Whole genome GATC coordination)。[DCL data server已经预先设置过了必要的参数，如无特殊需求，可以跳过此步骤。]{.underline}

要使用新的基因组文件组装index，使用以下命令：

```         
bowtie2-build <全基因组序列.fna> <Bowtie2 index名称>
perl ~/scripts/damidseq_pipeline/gatc.track.maker.pl <全基因组序列.fna> <output prefix>
```

-   bowtie2-build：传入fasta格式的全基因组序列，建立bowtie2 index。

-   gatc.track.maker.pl：传入fasta格式的全基因组序列，输出指定前缀的全基因组GATC坐标。

-   预先建立的bowtie index: /home/DiChen/reference_genome/bowtie2/Wbcel235_bowtie2

-   预先提取好的GTAC坐标：\~/reference_genome/damid_ref/wbcel235_gatc.GATC.gff

-   预先建立的index是使用Ensembl来源的WBcel235(ce11)线虫基因组建立的。如需进行其他分析，配套的文件为：

    1.  全基因组序列：\~/reference_genome/assembly/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa
    2.  gff3基因组注释文件：\~/reference_genome/gxf/Caenorhabditis_elegans.WBcel235.112.gff3
    3.  官方网址：<https://mart.ensembl.org/Caenorhabditis_elegans/Info/Index>

得到所有index文件后，使用以下命令重置damidseq_pipeline的参数：

```         
conda activate damid
damidseq_pipeline --gatc_frag_file=<GATC_file_path> --bowtie2_genome_dir=<Bowtie2_index_path>
```

***NOTICE: 如果没有修改的需求，请直接跳过以上步骤！***

## 2.2 使用damidseq_pipeline进行比对：

理论上，该软件可以直接输入原始测序文件进行分析。然而比对的步骤较为花费时间，不利于后期调整参数反复运行，因此我们选择先手动进行比对。

首先，需要把所有待比对的原始数据移动到一个单独的文件夹中，**将工作目录cd到这个文件夹后，确保名称符合要求(请参考第一步)。**

要使用damidseq_pipeline的比对模式，使用以下命令：

```         
conda activate damid
damidseq_pipe --just_align --paired
```

软件会自动寻找双端测序配对的两个文件，并使用bowtie2进行比对。

## 2.3 运行damidseq_pipeline

damidseq_pipeline无法处理多重复数据。对于2个及以上重复样本，需要手动对Dam_fusion组和Dam_only组进行匹配。例如，给定以下文件:

```         
D1.bam
D2.bam
D3.bam
D4.bam
D5.bam
D6.bam
```

其中D1/D2/D3为Fusion组，D4/D5/D6为control组。首先必须要在control组的开头加上"Dam\_"，否则软件将出现报错：

```         
D1.bam
D2.bam
D3.bam
Dam_D4.bam
Dam_D5.bam
Dam_D6.bam
```

对一组Fusion和control组文件进行手动配对，并输入软件开始运行：

```         
damidseq_pipeline --dam=Dam_D4.bam D1.bam
damidseq_pipeline --dam=Dam_D5.bam D2.bam
damidseq_pipeline --dam=Dam_D6.bam D3.bam
```

这里配对的方式不唯一，可以反复尝试，找到效果最好的配对方法。重复数也不一定需要三个，两个重复也足以进行分析，可以删去偏差最大的组。

## 2.4 damidseq_pipeline的输出解读

输出文件中，我们关注xx_vs_xx.kde.norm.gatc.bedgraph这个文件。如下图所示，该文件包括了四列：

[![pAkMJXR.png](https://s21.ax1x.com/2024/08/26/pAkMJXR.png)](https://imgse.com/i/pAkMJXR)

-   第一列：染色体编号
-   第二/三列：峰起始/结束位置坐标
-   第四列：峰得分，即Log2(Fusion/Dam)

对于多重复样本，我们将得到若干个bedgraph文件。使用bedgraph_avg.R合并多重复文件：

```         
Rscript /home/DiChen/scripts/damidseq_pipeline/bedgraph_avg.R <bedgraph1> <bedgraph2> <bedgraph3> <...> <outputfile.bedgraph>
```

该脚本依次读取多重复的bedgraph数据，对第四行的Log2得分进行平均，并输出平均后的bedgraph文件。**该文件即是下游分析需要的。**

# 3. Peak calling + Peak assignment

## 使用find_peaks.py进行peak calling

find_peaks软件使用bootstrap方法，对峰进行显著性计算，并输出P值校正后的显著峰。原github链接：<https://github.com/owenjm/find_peaks。该作者的原始代码存在bug，我进行了略微修改（因此直接从官网上下载可能无法跑通>）。

由于该软件python版本和damidseq_pipeline冲突，需要使用base环境运行该软件：

```         
conda activate base
python ~/scripts/damidseq_pipeline/find_peaks/find_peaks.py <bedgraph_file> --fdr=0.05 --n=100
```

-   fdr：p值矫正阈值，默认值为0.01。设置为0.05可以获得更多显著峰。
-   n：bootstrap次数，默认值为100。设置为更高的次数会使运行时间变长，结果可能变得更准确。

该软件将输出一个文件夹（以系统时间命名的，很长的名字）。其中xxx_allpeaks.gff是我们需要的结果文件，该文件包括八列：

[![pAkQKKA.png](https://s21.ax1x.com/2024/08/26/pAkQKKA.png)](https://imgse.com/i/pAkQKKA)

-   第一列：染色体编号

-   第二列：峰起始位置

-   第三列：峰结束位置

-   第四列：平均峰得分

-   第五列：峰得分

-   第六列：峰计数

-   第七列：峰宽度

-   第八列：FDR

## 使用peak2genes进行峰的分配

接下来，我们将把得到的显著峰分配给对应的基因。原github网址提供的perl脚本存在bug（或与我们的输入文件格式完全不兼容）。下面我们将使用我编写的山寨版本进行计算：

```         
/usr/bin/Rscript ./scripts/damidseq_pipeline/find_peaks/peak2genes.R \
--gxf=<Path to annotation file> \
--significant_peaks=<Path to significant peaks> \
--combine_mode=<max/随便什么> \
--anno_type=gene \
--assign_mode=<promoter/loose/center> \
--gene_pad=1000 \
--gene_up_pad=6000 \
--gene_down_pad=1000 \
```

-   gxf：基因组注释文件位置。建议与上游分析所使用的一致。

-   significant_peaks：输入的显著峰信息，即上一步得到的xxx_allpeaks.gff。

-   combine_mode：合并峰的模式。如果为max（默认为max），则对于一组重叠的峰，取最上游坐标和最下游坐标作为合并后峰的坐标。如果为其他（随便什么），则取公共重叠部分（最小部分）作为峰的新坐标。

-   anno_type：gxf文件“type”一列中，表示整个基因的字段。需要根据输入注释文件的不同修改。如使用默认gff3文件，请勿修改此参数。

-   assign_mode：分配峰的模式。可以为promoter/loose/center三种其一。

    1.  promoter：启动子模式。该模式认为基因TSS前后'gene_pad'长度为基因的调控部位，并将任意部位落在该区间的峰分配给对应基因。
    2.  loose：松弛模式。该模式认为基因TSS前'gene_pad'长度到TES后'gene_pad'长度的区域都是基因的调控部位，并将任意部位落在该区间内的峰分配给对应基因。
    3.  center: 中心模式。该模式认为基因TSS前'gene_up_pad'长度到TES后'gene_down_pad'长度的区域是基因的调控部位，并将**中心坐标**落在该区间内的峰分配给对应基因。

    如下图所示： [![pAk09HS.png](https://s21.ax1x.com/2024/08/27/pAk09HS.png)](https://imgse.com/i/pAk09HS)

# 4. 下下游分析：数据呈现

## 4.1 使用seqplots软件绘制profile图和结合热图

seqplots是一个基于R语言和shiny的交互式绘图软件。使用时，需要准备：

1.  第二步得到的bedgraph文件。
2.  第三步得到的显著峰分配的基因列表
3.  基因组注释文件

### 手动提取显著基因坐标

首先，我们需要从注释文件中提取想要呈现的显著基因的坐标：

```{r, eval=FALSE}
library(dplyr)
gxf <- rtracklayer::import("D://R-data/damid_seq/gxf/WBcel235.gff") %>% as.data.frame()
significant_gene <- c("geneA", "geneB", "geneC", "geneD")
newgxf <- gxf %>% filter(type == "gene" & Name %in% significant_gene)
rtracklayer::export.gff3(newgxf, "newgxf.gff")
```

### 上传文件

点击"Add files"按钮，进入上传文件界面。上传bedgraph文件和我们新得到的newgxf.gff即可。

[![pAk0yvt.png](https://s21.ax1x.com/2024/08/27/pAk0yvt.png)](https://imgse.com/i/pAk0yvt)

**要注意，seqplots本身存在bug，上传文件大概率出错。遇到这种情况，请反复重启软件，并反复上传，直到成功为止。**

### 绘图

点击'New plot set'按钮，进入绘图设置界面。选中我们刚才上传的track和annotation，如果是center模式，则参数按下图设置，其中upstream和downstream的缓冲区按照gene_up_pad和gene_down_pad的具体值而定。

[![pAk0WVS.png](https://s21.ax1x.com/2024/08/27/pAk0WVS.png)](https://imgse.com/i/pAk0WVS)

得到图像后，在绘图区可以调节图形参数，包括配色等。[原软件文档](https://przemol.github.io/seqplots/)有详细说明，这里不再赘述。
