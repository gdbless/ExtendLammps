# 5 理解 pair 样式

在本章中，我们将深入研究 LAMMPS 中实现的一些 pair 样式，并分析它们的源代码，以了解代码不同部分所扮演的角色。对样式根据*第 1 章，MD 理论和模拟实践*中介绍的方程式实现原子之间的对势。在本章中，我们将了解所涉及的数学运算是如何通过对式源代码执行的。

我们将涵盖以下主题：

- 查看配对样式的一般结构

- Morse 势

- Table 势

- DPD 势

在本章的最后，您将借助前面章节中介绍的材料探索几种样本对样式的内部机制，并准备对现有的对样式进行更改以在 LAMMPS 中实现自定义对势。

## 技术要求
要执行本章中的说明，您只需要一个文本编辑器（例如，**Notepad++** 或 **Gedit**）。

你可以在这里找到本章使用的完整源代码：https://github.com/PacktPublishing/Extending-and-Modifying-LAMMPS-Writing-Your-Own-Source-Code

## 回顾 pair 样式 的一般结构
如*第 3 章，源代码结构和执行阶段*中所述，每个单独的 pair 样式都继承自 `pair.cpp` 和 `pair.h`，包括 `init()` 方法。

在本节中，我们将简要介绍从 `pair.cpp` 和 `pair.h` 类继承的子对样式类中常用的方法。

父类负责验证对系数分配、混合参数、确定截止、请求邻居列表和设置计算。一些方法和变量也继承自`pair.h`。`pair_hybrid.cpp` 和 `pair_hybrid_overlay.cpp` 类提供混合样式。

子对样式类通常包含以下方法：

- `allocate()` 将内存分配给用于计算配对相互作用力的数组。

- `settings()` 读取并处理在 LAMMPS 输入脚本中的 `pair_style` 命令之后输入的全局对潜在参数。

- `coeff()` 读取并处理在 LAMMPS 输入脚本中的 `pair_coeff` 命令之后输入的本地对潜在参数。

- `init_one()` 为所有原子类型分配成对系数。

- `compute()` 计算原子类型之间的成对力和势能。

- `single()` 返回模拟运行期间其他类使用的成对势。

在以下部分中，我们将说明如何在 LAMMPS 存储库中提供的选定对样式中实现上述方法。

## 回顾Morse势

在本节中，我们分析 `pair_morse.cpp` 和 `pair_morse.h` 类实现的Morse势。

Morse势 $(V_{Morse})$ 用于表示键的共价特性，由排斥和吸引部分组成，两者之间有一个潜在的最小值。Morse势的函数形式由下式给出：

$$V_{Morse}(r)=D_0[e^{-2\alpha (r-r_0)}-2e^{-\alpha (r-r_0)}]$$

由上式可知，$D_0$ 表示Morse势阱的深度，$\alpha$ 决定了阱的曲率，$r_0$ 表示平衡距离。在 LAMMPS 输入脚本中，Morse对势是使用以下命令实现的：

```
pair_style morse GLOBAL_CUTOFF
pair_coeff TYPE1 TYPE2 D0 ALPHA R0 LOCAL_CUTOFF
```

从前面的代码中，D0、ALPHA、R0 是描述的三个Morse参数，而 `GLOBAL_CUTOFF` 和 `LOCAL_CUTOFF` 参数分别是用于所有原子对或指定原子对的截止值（在 LAMMPS 手册中有详细描述： https://lammps.sandia.gov/doc/pair_morse.html）。TYPE1 和 TYPE2 参数指定通过Morse势相互作用的原子类型对。

在接下来的部分中，我们将说明 `pair_morse.cpp` 和 `pair_morse.h` 如何促进在 MD 模拟中解析和实现Morse势所需的处理。

### PairMorse::allocate()
在实现pair style Morse时，需要的数组由`allocate()`方法创建，如下所示：

<div align=center>
<img src=./ch05/fig5_1.jpg>
</div>

图 5.1 - 来自 `pair_morse.cpp` 的代码片段显示了 `allocate()` 方法

从前面的截图可以看出，原子类型的数量存储在*第 136 行*为 n，并且使用 ``memory ->创建（）``函数。

每行和每列都从等于 0 的索引开始，并继续索引 n 以表示行和列中的所有原子类型。*第 138* 到 *141* 行插入零以填充 2D setflag 数组（参见 *Chapter 3, Source Code Structure and Stages of Execution* 中的 `pair.cpp`），它检查是否已在每对原子类型。

您可以在 *line 139* 中看到，填充数组从索引 i = 1 开始，以对应于也从索引 1 开始的原子类型（但不包括索引 i = 0）。在*第 143* 到 *150* 行创建了类似的二维数组来表示截断的平方 (`cutsq`)、截断 (`cut`)、$D_0$(`d0`)、$\alpha (`alpha `)$、$r_0$ (`r0`)、基于 $D_0(d0)$ 和 $\alpha$ (`morse1`) 的常数，以及在截止点将电位移至零的偏移量 (`offset` ) 对于每对原子类型。这些数组将由其他方法填充并用于计算配对交互。

现在让我们继续下一个方法。

### PairMorse::settings()
此方法在 `pair_style morse` 命令之后读取在 LAMMPS 输入脚本中输入的量，在这种情况下，它只是全局截止值。以下屏幕截图说明了 `settings()` 方法：

<div align=center>
<img src=./ch05/fig5_2.jpg>
</div>

图 5.2 - 来自 `pair_morse.cpp` 的代码片段显示了 `settings()` 方法

如您在前面的屏幕截图中所见，*第 159 行* 检查 LAMMPS 输入脚本中是否只有一个参数存在于 LAMMPS 手册中指定的 `pair_style` Morse命令之后，如果没有则返回错误消息。

然后，*第 161 行* 使用 ``force->numeric()`` 方法将此 arg[0] 参数作为浮点数读取，并将其存储为已在 `pair_morse.h 中声明的 `cut_global` 变量`。在*第 165* 到 *170* 行中，之前分配的所有原子类型 (`cut[i][j]`) 的截止值都设置为这个全局截止值。请注意，i 索引遍历所有原子类型 1 到 n，但 j 索引遍历从 i 到 n，这意味着二维数组切割中只有一半的数组元素被填充。其余部分将由接下来描述的 `pair_coeff()` 方法填充。

### PairMorse::coeff()
此方法在 `pair_coeff` 命令之后读取在 LAMMPS 输入脚本中输入的数量。以下屏幕截图显示了 `coeff()` 方法：

<div align=center>
<img src=./ch05/fig5_3.jpg>
</div>

图 5.3 - 来自 `pair_morse.cpp` 的代码片段显示了 `coeff()` 方法

从前面的屏幕截图中可以看到，*第 179* 到 *180* 行允许在 LAMMPS 输入脚本中的 `pair_coeff` 命令后添加 5 或 6 个参数。*第 183* 到 *185* 行使用 ``force->bounds()`` 方法读取原子类型对，此势适用，并将指定原子类型的范围分配给变量，`ilo`，` ihi`、`jlo` 和 `jhi`。

第 187 到 189 行将第三个、第四个和第五个参数分别存储为 $D_0$、$\alpha$ 和 $r_0$ 的浮点值。*第 191* 到 *192* 行将全局截止设置为默认截止，除非输入了覆盖指定原子类型对的截止的第六个参数。与这些数量相对应的二维数组 `d0`、`alpha`、`r0` 和 `cut` 在 *194* 到 *204* 行中更新，同时将这些原子的 `setflag` 元素更改为 1类型对。这些数组中只有一半被填充，其余的将被填充到 `PairMorse:init_one()` 中。计数整数用于记录分配了这些参数的原子类型对的数量，如果没有计数对，*第 206 行*返回错误。

### PairMorse::init_one()
此方法读入两个原子类型作为参数，并验证为它们定义了对系数。它还定义了存储在 morse1 数组中的数量，计算在截止距离处将势值移至零所需的偏移量，并通过从一组原子类型复制到其镜像原子类型集来填充数组。

以下屏幕截图演示了 init_one() 函数：

<div align=center>
<img src=./ch05/fig5_4.jpg>
</div>

图 5.4 - 来自 `pair_morse.cpp` 的代码片段显示了 init_one() 方法

如您所见，在*第 216 行* 验证作为参数输入的 i 和 j 原子类型是否为它们的交互分配了对系数后，*第 218 行* 将 morse1 量定义为 $2\alpha D_0$，以便计算以后就业。

第 220 到 223 行计算截止处的势能，即 $V_{Morse}(r_{cut})$ ，并将其存储为偏移量，可以添加该偏移量以将势能函数在截止点处精确地移动到零。如果 LAMMPS 输入脚本没有在势能中实现任何偏移，则 `offset` 数组用零填充。然后，在 *225* 到 *229* 行中，通过将每个原子对 (i,j) 的元素复制到其镜像对 (j,一世）。

### PairMorse::compute()
在填充了所有原子对的所有数组元素后，在此方法中计算了这些原子之间的成对力和势能。以下屏幕截图向我们展示了 `compute()` 方法：

<div align=center>
<img src=./ch05/fig5_5.jpg>
</div>

图 5.5 - 来自 `pair_morse.cpp` 的代码片段显示了 `compute()` 方法

如您所见，*第 79* 到 *89* 行循环使用局部索引 i 遍历所有中心原子，然后使用局部索引 j 遍历所有相邻原子，如 *Chapter 4, Accessing Information by Variables 中所述，数组和方法*。*第 90 行* 检查 `i` 和 `j` 原子之间的特殊键，`factor_lj` 变量说明了在 LAMMPS 输入脚本中分配给特殊键的权重。该检查是通过使用 `sbmask()` 函数的按位比较来执行的，这将在*附录 B，使用 GDB 和 Visual Studio 代码调试*中进行解释。如果指定了特殊债券，则将 `factor_lj` 设置为等于债券的权重；否则，它被设置为等于 1。

第 93 到 95 行计算从第 j 个原子到第 i 个原子的位移矢量分量 $(\Delta {x},\Delta {y},\Delta {z})$，*第 96 行* 计算距离的平方在这些原子之间， $rsq=(\Delta {x})^2+(\Delta {y})^2+(\Delta {z})^2$ 。除非需要规避与平方根计算相关的计算开销，否则不会计算距离`r`本身。在*第 99 行*，`rsq` 与这两个原子之间的截止值的平方进行比较，以确定原子是否位于截止值内。如果是这样，`r`在*第 100 行*中明确计算，然后计算两个原子之间的力。

如*第 1 章，MD 理论和模拟实践*中所述，成对力取决于成对势的梯度，在这种情况下为 $V_{Morse}(r)$。

首先，我们计算以下数量：

$$\frac{dV_{Morse}}{dr}=-2\alpha D_0[e^{-2\alpha (r-r_0)}-e^{-\alpha (r-r_0)}]$$

然后，可以计算 $(x, y, z)$ 力分量。(x) 力分量由下式给出：

$$F_x=\frac{dV}{dr}\times \frac{\Delta{x}}{r}=2\alpha D_0[e^{-2\alpha (r-r_0)}-e^{- \alpha (r-r_0)}] \frac{\Delta {x}}{r}$$

$(y)$ 力分量是这样计算的：

$$F_y=\frac{dV}{dr}\times \frac{\Delta{y}}{r}=2\alpha D_0[e^{-2\alpha (r-r_0)}-e^{- \alpha (r-r_0)}] \frac{\Delta {y}}{r}$$

$(z)$ 力分量由以下给出：

$$F_z=\frac{dV}{dr}\times \frac{\Delta{z}}{r}=2\alpha D_0[e^{-2\alpha (r-r_0)}-e^{- \alpha (r-r_0)}] \frac{\Delta {z}}{r}$$

在源代码（*图 5.5*）中，*第 101* 到 *102* 行定义了变量 $d_r=r-r_0$ 和 $dexp=e^{-\alpha(r-r_0)}$，以及 *第 103 行* 根据 `factor_lj`、`morse1` 和 `dexp` 定义了 `fpair` 变量。在这里替换 morse1（参见 `PairMorse::init_one()` 部分）后，`fpair` 结果如下：

$$fpair=factor_lj*2\alpha D_0[e^{-2\alpha (r-r_0)}-e^{-\alpha (r-r_0)}]\frac{1}{r}$$

通过将 fpair 与位移矢量分量 $(\Delta {x},\Delta {y},\Delta {z})$ 相乘，在 *105* 到 *107* 行中计算了由 `factor_lj` 缩放的相应力分量。这些力被分配给力数组元素中的 `i` 原子，`f[i][0]`、`f[i][1]` 和 `f[i][2]`（`f ` 数组已在 *line 66* 中定义）。假设牛顿对被激活，作用在相反方向的反作用力分量被分配给*第 108* 到 *112* 行中的`j`原子。

如果需要，此方法还可以计算原子之间的相互作用能量，其中包括偏移量和`factor_lj`的缩放比例（*第 114* 到 *118* 行）。

### PairMorse::single()
运行 MD 模拟时，其他类（例如`computes()`）通常需要交互能量。可以使用`single()`方法计算和访问，如下所示：

<div align=center>
<img src=./ch05/fig5_6.jpg>
</div>

图 5.6 - 来自 `pair_morse.cpp` 的代码片段显示了 `single()` 方法

如您所见，输入参数用于定义与前面描述的原子 `types`、`itype` 和 `jtype` 相同的变量、dr、dexp 和 `fforce`（相当于 `fpair`），进入。然后，在* 348* 到*349* 行，能量作为根据势函数$V_{Morse}(r)$ 及其偏移和缩放因子`factor_lj` 计算的`phi` 变量返回。

其他对样式，例如 **Lennard-Jones 势** (`pair_lj_cut.cpp`)、**Buckingham 势** (`pair_buck.cpp`) 和非长程**库仑势** (`pair_coul_cut .cpp`) 采用类似的算法，具有不同的势函数和力函数实现。在编写实现位置相关电位的自定义对样式时，建议选择现有的对电位之一，例如 `pair_morse.cpp` 并修改相关部分。

在本节中，我们了解了负责在 MD 模拟中实现Morse势的源代码。

在下一节中，我们将分析从表格中插入力和势的表格势，而不是使用函数进行计算。

## Table 势

在本节中，我们将了解实现对表势的各种方法。

`pair_table.cpp` 和 `pair_table.h` 描述的表对电位读取一个文件，其中包含与原子之间的分离相对应的数据点列表的电位和力值列表。当需要特定分离处的力或势能时，文件中最接近所需分离的一个或多个数据点用于内插力和势能。

以下 LAMMPS 输入脚本命令用于实现这种对样式：

```
pair_style table TABLE_STYLE N KEYWORD
pair_coeff TYPE1 TYPE2 FILE_NAME VARIABLE CUTOFFs
```

从前面的代码中，`TABLE_STYLE` 必须是 LAMMPS 手册中提到的样式之一（`LINEAR`、`LOOKUP`、`SPLINE` 或 `BITMAP`），`N` 是数据点的数量和KEYWORD 是手册中描述的远程求解器的可选条目。

在下一行中，TYPE1 和 TYPE2 代表与该势能交互的原子类型，`FILE_NAME` 是包含数据点的文件的名称，VARIABLE 是文件中用于标识文件开头的关键字，` CUTOFF` 已经解释过了。LAMMPS手册（https://lammps.sandia.gov/doc/pair_table.html）中提供了数据文件格式，如下图：

<div align=center>
<img src=./ch05/fig5_7.jpg>
</div>

图 5.7 - LAMMPS 手册中的表电位数据文件格式

如您所见，数据文件中的关键字必须与输入脚本中的`VARIABLE`相同，四列分别代表用于插值的索引、分隔（`r`）、势能和力。

在概述了输入脚本语法和数据文件格式之后，我们分析了 `pair_table.cpp` 中的各种方法，在接下来的部分中，处理这种配对样式，从通过 `PairTable::settings()` 解析全局系数开始.

### PairTable::settings()
此方法读取在 pair_style table 命令之后输入的全局系数。以下屏幕截图向我们展示了 `settings()` 方法：

<div align=center>
<img src=./ch05/fig5_8.jpg>
</div>

图 5.8 - 来自 `pair_table.cpp` 的代码片段显示了 `settings()` 方法

如您所见，此方法允许两个或多个全局系数（*第 216 行*）。表格样式变量 `tabstyle`（*第 220* 到 *224* 行）被分配一个整数以对应于允许的表格样式之一，这些样式在 `pair_table.h`（*第 43* 行）中枚举：

```
enum{LOOKUP,LINEAR,SPLINE,BITMAP};
```

第 226 到 227 行读取数据点的数量（`tablength`），至少需要两行数据。*如果在 LAMMPS 输入脚本中指定，*第 232* 到 *241* 行实现远程交互。

### PairTable::coeff()
此方法读取为每对原子类型输入的局部系数。下一个屏幕截图举例说明了 `coeff()` 方法：

<div align=center>
<img src=./ch05/fig5_9.jpg>
</div>

图 5.9 - 来自 `pair_table.cpp` 的代码片段显示了 `coeff()` 方法

正如你所看到的，`coeff()` 方法只允许 4 或 5 个参数（*第 265* 行）并将前两个参数读取为原子类型（*第 269* 到 *270* 行）。*第 272* 到 *279* 行使用位于同一类中的 `read_table()` 方法读取用于识别表的文件名和变量名。

`read_table()` 方法将文件名和变量名作为参数，打开此文件并循环遍历其行以定位变量名并使用 `param_extract( )`方法。完成后，根据输入参数计算或准备从文件中读取 r 变量以及能量和力值。

之后（*第 320 行*），`coeff()` 方法调用 `compute_table()` 方法计算表格内界（`rinner`）、表格外界（`cut`）、间距的平方表数据点（`delta`）和间距的倒数（`invdelta`）。

这里，`invdelta` 是 `delta` 的倒数，它被计算一次并存储为一个变量，以避免重复计算增加与除法操作相关的过载。在同一个 `compute_table()` 方法中，力和能量值分别存储在一维数组 f 和 e 中。在存储到 `f` 数组之前，为力获得的值除以 r，因为它有助于在 `compute()` 方法中计算力的 x、y、z 分量，接下来将讨论。

### PairTable::compute()
如前所述，此方法计算原子与其每个相邻原子之间的间距，并将间距的平方存储为 rsq。使用这个 rsq 值，给出最接近这个分隔而不超过它的 `r` 值的表条目索引被确定为整数，`itable`。

以下屏幕截图显示了确定 itable 的代码部分：

<div align=center>
<img src=./ch05/fig5_10.jpg>
</div>

图 5.10 - 来自 `pair_table.cpp` 的代码片段显示了 `compute()` 方法

如您所见，力的计算取决于 tabstyle 变量，该变量指示 `LOOKUP`、`LINEAR`、`SPLINE` 或 `BITMAP` 样式之一：

- 对于 `LOOKUP` 样式（*第 129 行*），直接从 `f` 数组的 `itable` 条目中读取 force，即 f[`itable`]。

- 对于 `LINEAR` 样式（*第 137* 到 *138* 行），索引处的 `r` 值之间的小数距离，`itable` 和 (`itable`+1)，存储为小数变量和用于在相同的两个索引处的`f`值之间进行线性插值，使用包含表中每两个连续`f`值之间差异的 df 数组（参见*第 676 行*）。

- 对于 `SPLINE` 样式（*第 147* 到 *151* 行），分数距离的计算与 `LINEAR` 样式相同，但被输入三次函数（*第 149* 行）以计算力。

- 对于 `BITMAP` 样式（*第 153* 到 *159* 行），执行手册中描述的算法，使用位将相关表标识为 `itable` 并进行相应的插值。

确定的力存储为 `fpair` 并乘以位移矢量分量 $(\Delta {x},\Delta {y},\Delta {z})$ 以获得力分量（*第 162* 到 *168 行*)。相应的能量在 *171* 到 *181* 行中计算：

在本节中，我们了解了用于实现对表势的方法，即 `settings()`、`coeff()` 和 `compute()` 方法。

在下一节中，我们将分析应用成对、阻力和随机力的成对样式 DPD。

## DPD 势

在本节中，我们通过 `pair_dpd.cpp` 和 `pair_dpd.h` 类查看 DPD 潜力及其实现。

**耗散粒子动力学 (DPD)** 涉及与耗散力和随机力耦合的成对保守力，这些力和随机力作用在用于表示较大分子或团簇的两个粒子上。与传统的 MD 相比，通过粗粒度消除或最小化分子或簇的原子细节，以便在更长的时间范围内进行模拟。这种技术在模拟流体时特别有用，其中粒子代表分子或流体块而不是原子。耗散力和随机力可分别用于模拟流体阻力和碰撞力。

在实现 DPD 势时，来自成对势的力占成对力的一部分，而耗散力需要使用相对粒子速度来计算，并且随机力必须遵循高斯分布。总之，作用在第 i 个粒子上的 DPD 力由三个成对加性力 $(F^C,F^D,F^R)$ 与相邻粒子 j 的总和给出，该相邻粒子 j 位于固定截止距离 $( r_C)$:

$$\overrightarrow{f} = (F_{ij}^C+F_{ij}^D+F_{ij}^R)\hat {r}_{ij}$$

单位向量，$\hat {r}_{ij}=\frac{\hat {r}_{i}-\hat {r}_{j}}{|\hat {r}_{i}- \hat {r}_{j}|}$ ，从粒子 j 指向粒子 i。代表粒子化学性质的保守力 $(F_{ij}^C)$ 提供了最大量级为 A 的软线性斥力：

$$F_{ij}^C = A(1-\frac{|\hat{r}_{ij}|}{r_c}$$

取决于两个粒子的速度差的耗散力$(F_{ij}^D)$，$\hat {v}_{ij}=\hat {v}_{i}-\hat {v} _{j}$，阻力系数 $\gamma$ 由下式给出：

$$F_{ij}^{D}=-\gamma _{ij}(1-\frac{|\hat {r}_{ij}|^2}{r_c})(\hat {r}_{ ij}\times \hat {v}_{ij})$$

耗散力和随机力$(F_{ij}^R)$的值由**涨落耗散定理**相关，统计值与系统温度(T)在高斯分布上一致分配。使用均值为 0，方差为 1 的高斯随机数 $\alpha$，$F_{ij}^R$ 由以下等式给出：

$$F_{ij}^{R}=\alpha \sqrt{2k_BT\gamma _{ij}}(1-\frac{|\hat {r}_{ij}|}{r_c})(\Delta { t})^{-1/2}$$

根据前面的等式，Δ𝑡𝑡 是时间步长，$k_B$ 是玻尔兹曼常数。在 LAMMPS 输入脚本中，以下命令实现 DPD 电位：

```
pair_style dpd TEMPERATURE GLOBAL_CUTOFF SEED
pair_coeff TYPE1 TYPE2 A GAMMA LOCAL_CUTOFF
```

从前面的代码中，TEMPERATURE 参数是不言自明的，A 和 `GAMMA` 参数分别表示前面描述的量 A 和 𝛾𝛾。`GLOBAL_CUTOFF` 和 `LOCAL_CUTOFF` 参数分别是用于所有原子对或指定原子对的截止值（在 LAMMPS 手册中有详细描述：https://lammps.sandia.gov/doc/pair_dpd.html）。`TYPE1` 和 `TYPE2` 参数指定通过 DPD 势相互作用的原子类型对。SEED 参数是一个整数，用于生成随机数的高斯分布，在 `settings()` 方法中进行了说明。

### PairDPD::settings()
此方法读取 DPD 电位的全局输入参数并初始化随机数生成器。以下屏幕截图显示了 `PairDPD::settings()` 方法：

<div align=center>
<img src=./ch05/fig5_11.jpg>
</div>

图 5.11 - 来自 `pair_dpd.cpp` 的代码片段显示了 `settings()` 方法

如您所见，*第 193* 到 *195* 行将全局参数`TEMPERATURE`和`GLOBAL_ CUTOFF`读取为浮点数，将`SEED`参数读取为整数。在 *line 201* 中，SEED 用于创建随机类（在 `pair_dpd.h` 中声明），它调用 `random_mars.cpp` 类中编码的 **Marsaglia** 随机数生成器（`random_mars.h ` 头文件已在 *line 27* 中导入）。

**Marsaglia** 随机数生成器采用整数参数并对其进行处理以生成随机数序列。`SEED` 被添加到分配给使用的多核（通过 `comm->me` 访问）中的每个内核的处理器等级中，以便将不同的整数馈入每个内核。随机类将用于生成所需分布中的随机数，如下一节所述。

### PairDPD::compute()
此方法计算 DPD 对相互作用中涉及的力。在建立每个 i 粒子的邻居之后，它会遍历这 j 个邻居并计算三个力贡献 $(F^C,F^D,F^R)$。以下屏幕截图向我们展示了 `PairDPD::compute()` 方法：

<div align=center>
<img src=./ch05/fig5_12.jpg>
</div>

图 5.12 - 来自 `pair_dpd.cpp` 的代码片段显示了 `compute()` 方法

如您所见，对于位于分界线内的邻居，计算了 r 间隔（*第 114 行*）和倒数间隔`rinv`（*第 116 行*），然后是 $(x, y, z )$ $\vec {v}_{ij}$ 的分量（*第 117* 到 *120* 行），其中 v 数组包含*第 75* 行中声明的粒子速度和变量（`vxtemp`、`vytemp `, `vztemp`) 表示粒子的速度分量，`i` (*lines 95* to *97*)。

*第 120 行* 将点积 $(\vec {r}_{ij} \times \vec {v}_{ij})$ 存储为变量 `dot`，需要除以 r 分隔为求$F^D$中所要求的数量$(\vec {r}_{ij} \times \vec {v}_{ij})$。*第 121 行* 定义了一个变量 `wd`，用作数量 $1-\frac{|\vec {r}_{ij}|}{r_c}$ 的简写符号。randnum 变量在*第 122 行* 中定义，通过使用 `random->gaussian()` 方法（定义在 `random_mars.cpp ）。

定义变量后，$F^C$ 在*第 128 行* 中计算并存储为变量 `fpair`。在*第 129 行*，在计算 $F^D$ 时，将 rinv 数量合并为 dot 除以 r，并将结果计算为 `fpair`。在计算 $F^R$ 之前，2D 数组 sigma 在 PairDPD::init_one() 内的 *line 273* 中定义以存储数量 $\sqrt {2k_BT\gamma _{ij}}$：

```
sigma[i][j] = sqrt(2.0*force->boltz*temperature*gamma[i][j]);
```
此外，在 *line 81* 上定义了 dtinvsqrt 变量来存储数量 $(\Delta {t})^(-\frac{1}{2})$：

```
双 dtinvsqrt = 1.0/sqrt(update->dt);
```

在*第 130 行*，这些数量与 randnum 和 wd 相乘以获得𝐹𝐹𝑅𝑅，并将结果计算为 fpair。最后，在*第 131 行*，包含来自 $F^C$、$F^D$ 和 $F^R$ 贡献的 `fpair` 变量按 `factor_dpd` 数量缩放以考虑任何特殊债券（*第 104 行 *) 并乘以 rinv 以方便 $(x, y, z)$ 力分量计算。在这个阶段，`fpair` 实际上等于以下内容：

$$fpari=factor_dpd*(F_{ij}^C+F_{ij}^D+F_{ij}^R)\frac{1}{r_{ij}}$$

`fpair` 的 $(x, y, z)$ 分量在 *133* 到 *139* 行中被分配给力数组 f 的适当数组元素（参见 *图 5.12*）。

DPD 势的势 $V_{DPD}(r)$ 仅根据保守力 $F^C$ 计算：

$$V_{DPD}(r)=-\int _{0}^{r} {F^C(r)\ dr}=-\int _{0}^{r} {A(1-\frac {r}{r_c})\ dr}=\frac{1}{2}Ar_c(1-\frac{r}{r_c})^2$$

这个势在 *142* 到 *148* 行和 `PairDPD::single()` 方法中计算，如果需要，可以缩放以适应特殊键。

在本节中，我们了解了 DPD 对势背后的机制，包括由 `compute()` 方法执行的力和势计算。

## 总结

在本章中，分析了一些现有的 LAMMPS 对样式，以展示各种方法在执行所需计算和交换信息中所起的作用。在 LAMMPS 中实现的势的其他变体包括多体（即非成对）**Stillinger-Weber (SW)** 势和将对系数分配给单个原子对而不是原子类型对的 Pair List 选项。

位置相关势函数使用 `compute()` 方法中的函数形式实现，或者可以通过 pair 样式读取的表格形式实现。通过定义适当的力函数，也可以同时实施非保守力。

您现在可以使用从本章中学到的课程和技能来修改成对的适当方法，从而根据您的模拟要求对其进行自定义。

在下一章中，我们将分析选定的计算样式来说明代码的内部机制，正如本章所做的那样。

## 问题

1. 哪些方法用于读取pair style中的全局和局部系数？

2.为什么`fpair`变量总是有力函数乘以分离的倒数？

3. `single()` 方法的主要目的是什么？