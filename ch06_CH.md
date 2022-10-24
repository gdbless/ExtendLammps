# 6 理解计算

在上一章中，我们探讨了一些在 LAMMPS 中实现的 pair 样式。

在本章中，我们将分析一些计算并对它们的源代码有类似的理解。计算从模拟运行中返回各种数量，如果我们希望为自定义目的编写新的计算，学习理解负责相关计算的代码是必不可少的。

我们将在本章中介绍以下主题：

- 查看compute的一般结构

- 查看compute KE 类

- 查看compute group/group类

- 查看compute RDF 类

- 查看compute heat flux类

在本章结束时，您将了解所选计算的内部机制，并能够对 LAMMPS 源代码进行进一步修改。

## 技术要求

要执行本章中的说明，您只需要一个文本编辑器（例如，**Notepad++** 或 **Gedit**）。

你可以在这里找到本章使用的完整源代码：https://github.com/PacktPublishing/Extending-and-Modifying-LAMMPS-Writing-Your-Own-Source-Code

这是下载 LAMMPS 的链接：https://lammps.sandia.gov/doc/Install.html。LAMMPS GitHub 链接是 https://github.com/lammps/lammps ，也可以在其中找到源代码。

## 查看计算的一般结构
在本节中，我们将简要介绍计算子类中最常用的一些方法。

与对样式类似，单个计算继承自`compute.cpp`和`compute.h`类描述的父计算类。这些父类从 LAMMPS 输入脚本中读取前三个参数（计算 ID、组 ID 和计算样式）。以下屏幕截图显示了 `compute.h` 中继承的一些变量和数组：

<div align=center>
<img src=ch06/fig6_1.jpg>
</div>

图 6.1 – 来自 compute.h 的代码片段

子计算类可能包含以下一种或多种方法：

- init() 方法设置类并执行初步验证检查。

- init_list() 方法设置邻居列表或指向邻居列表的指针。

- compute_scalar() 方法计算通常用作输出的标量

- compute_vector() 方法计算通常用作输出的向量。

**重要的提示：**
鉴于使用计算的应用范围很广，许多其他方法可以在单个计算中实现。

接下来，我们将查看一些精选计算来说明这些计算背后的代码。

## 查看计算 KE 类
在本节中，我们将研究在计算`KE`类中实现的方法，该类计算一组原子的*动能*。

质量为 m 和线速度 $(v_x, v_y, v_z)$ 的原子的**动能 (KE)** 使用以下等式计算：

$$KE=\frac{1}{2}m(v_x^2+v_y^2+v_z^2)$$

`compute_ke.cpp` 类可用于计算指定原子组的动能。在一个 LAMMPS 输入脚本中，对应的语法如下：

```
compute COMPUTE_ID GROUP_ID ke
```

如您所见，`COMPUTE_ID`参数是定义的compute的唯一 ID，而 GROUP_ID 参数是compute所作用的原子组的 ID。这在 LAMMPS 手册 (https://lammps.sandia.gov/doc/compute_ke.html) 中有详细描述。首先，我们来看看这个类的构造方法，`ComputeKE::ComputeKE()`。

### ComputeKE::ComputeKE()
此构造方法继承自 `Compute` 父类（*第 26 行*）并检查参数的数量，如以下屏幕截图所示：

<div align=center>
<img src=ch06/fig6_2.jpg>
</div>

图 6.2 – compute_ke.cpp 中显示构造函数的代码片段

对于这个compute，只允许三个参数（compute ID、group ID 和compute样式）（*第 28 行*），它们都在 `compute.cpp` 中读取。由于 KE 作为标量返回，因此在 *line 30* 和 *31* 中激活了 `scalar_flag` 和 `extscalar` 变量。

接下来，我们将看一下 init() 方法。

**ComputeKE::init()**
`init()` 方法用于初始化此compute `KE` 类并定义转换为正确能量单位的转换因子（`pfactor`），如 *Chapter 4, Accessing Information by Variables, Arrays, and方法*。此方法显示在以下屏幕截图中：

<div align=center>
<img src=ch06/fig6_3.jpg>
</div>

图 6.3 – compute_ke.cpp 中显示 init() 方法的代码片段

KE 计算的大部分是由 `compute_scalar()` 方法执行的，接下来将讨论。

### ComputeKE::compute_scalar()
此方法计算指定组中每个原子的 KE，并找到所有原子的总和，如以下屏幕截图所示：

<div align=center>
<img src=ch06/fig6_4.jpg>
</div>
图 6.4 – 来自 compute_ke.cpp 的代码片段显示了 compute_scalar() 方法

如您所见，*第 56 行* 检查是否激活了每个原子的质量（`rmass`），这会将单个质量分配给单个原子。代码循环遍历核心中的所有原子（*第 57 行*），识别属于指定组的原子（*第 58 行*），并使用其 rmass 值和速度计算并附加每个原子的 KE（*第 59 行*）。

如果每个原子的质量是根据其原子类型定义的，则遵循相同的程序（*第 60 行*到 *65*），除了在 KE 计算中使用原子类型的质量（*第 63 行* ）。

`ke` 变量存储核心中所有原子的组合 KE，所有核心的 KE 总和使用 `MPI_Allreduce()` 方法（*第 67 行*）计算，该方法接受以下参数：

- `&ke`：这是将核心中的局部总和作为输入传递的变量。

- `&scalar`：这是存储来自所有核心组合的全局输出的变量。

- `1`：这是输入和输出数组的大小。

- `MPI_DOUBLE`：这表示输入和输出量的数据类型。

- `MPI_SUM`：这是要执行的运算类型（即求和）。

- `world`：这表示执行操作的核心。

**Message Passing Interface (MPI)** 是编译器的一个属性（不专属于
LAMMPS），允许内核之间的信息交换。有关更多信息
MPI 可在 *附录 C，熟悉 MPI* 中找到。相关头文件
可以使用以下库导入`compute_ke.cpp`：

```
#include <mpi.h>
```

除了求和之外，`MPI_Allreduce()` 方法还允许在多个内核上计算大量归约操作，包括求数量的最大值、最小值和乘积。

最后，从 `MPI_Allreduce()` 方法获得并存储为 `scalar` 变量（继承自 compute.h，如 *图 6.1* 所示）的 KE 的全局总和通过乘法转换为适当的单位它通过转换因子 pfactor (*line 68*)，并返回它 (*line 69*)。

在本节中，我们分析了一种返回标量的相对较短的计算方式。

在下一节中，我们将研究计算两组原子之间的相互作用能和力的计算。

## 查看compute gorup/group 类

在本节中，我们将分析一个更精细的计算，它使用更多的
方法。

两组原子之间的相互作用能和力可以使用`compute group/group`获得，并由`compute_group_group`实现。cpp` 和 `compute_group_group.h` 类。实现此计算的 LAMMPS 输入脚本命令如下：

```
compute COMPUTE_ID G1 group/group G2 
```

`COMPUTE_ID` 参数是compute的唯一 ID，而 `G1` 和 `G2` 参数是计算所作用的原子组的 ID（参见手册：https://lammps.sandia.gov/doc/compute_group_group.html)。手册中列出的可选参数关键字可以在这些参数之后输入，以指定其他选项，例如相互作用类型（*对电位或静电*）和分子 ID（相同或不同）。构造函数方法容纳了这些可选参数，我们将在以下部分中看到。

### ComputeGroupGroup::ComputeGroupGroup()
此构造函数至少接受四个参数（*第 50 行*），其中三个由 `compute.cpp` 读取。标量和矢量输出类型被指定（*第 52* 到 *55* 行）以分别适应能量和力输出。

*第 57* 到 *59* 行从输入脚本命令中读取第四个参数，并将其存储为 `group2` 字符串以表示第二组原子。然后，将对应的组 ID 定位为 `jgroup` 变量（*第 61 行*），如果组不存在（*第 62* 到 *63* 行），则会生成错误。在*第 64 行*，组的按位表示生成为 `jgroupbit`，用于识别属于该组的原子（请参阅附录 B，调试程序，了解有关按位表示的更多信息）。

以下屏幕截图显示了 ComputeGroupGroup() 构造函数：

<div align=center>
<img src=ch06/fig6_5.jpg>
</div>

图 6.5 – 来自 compute_group_group.cpp 的代码片段显示了构造方法

此compute的默认设置与可选关键字有关，如 LAMMPS 手册中所述（即 `pair = yes`、`kspace = no`、`boundary = yes` 和 `molecule = off`），由 *line 66* 到 *69* 中的标志值 0 或 1 促进。*line 71* 到 *105* 包含任何可选关键字，它们循环遍历所有参数以查找允许的关键字并根据需要调整标志值。

解析完参数后，`init()` 方法会执行一些错误检查，并在需要时请求邻居列表。

### ComputeGroupGroup::init() 和 init_list()
`init()` 方法检查指定的组之间是否存在对样式和静电相互作用（如果它们被启用），并且还启动相应的对和 kspace 对象（*第 122* 行*到*144*）。如果需要compute配对样式交互，则请求邻居列表（*第 167* 到 *172* 行），如以下屏幕截图所示：

<div align=center>
<img src=ch06/fig6_6.jpg>
</div>

图 6.6 – compute_group_group.cpp 中显示 init() 方法的代码片段

`init_list()` 方法使用 `compute_group_group.h` 中定义的列表指针提供对指向邻居列表的指针的访问，如以下屏幕截图所示：

<div align=center>
<img src=ch06/fig6_7.jpg>
</div>

图 6.7 – compute_group_group.cpp 中显示 init_list() 方法的代码片段

在创建了对邻居列表的访问之后，我们现在来看看 `pair_contribution()` 方法，它计算成对的能量和力。

### ComputeGroupGroup::pair_contribution()
如果启用了配对标志并且需要计算两组之间的配对电位，则调用此方法。首先，建立邻居列表（*第 230* 到 *235* 行），然后代码循环所有属于第一组 (`G1`) 的原子`i` 及其邻居`j`。

仅选择属于第二组 (`G2`) 的邻居`j`（*第 261* 到 *271* 行）。根据分子关键字，在 *275* 到 *281* 行检查分子 ID，以确定每个原子`j`是否与原子`i`属于同一分子。以下屏幕截图显示了能量和力的计算：

<div align=center>
<img src=ch06/fig6_8.jpg>
</div>

图 6.8 – 来自 compute_group_group.cpp 的代码片段显示了 pair_contribution() 方法

`i` 和 `j` 之间的对势能相互作用的能量和力是通过调用相应对样式中的 `single()` 方法计算的（*第 290 行*），并计入一维数组的第一个元素，一个（*第 296 行*）。作用在 `i` 上的力分量在 `one` 数组的下三个元素（`one[1]`, `one[2]`, `one[3]`）中计算*第 297* 到 * 行306*。数组元素在所有核心上求和（*第 324 行*），如下所示：

```
MPI_Allreduce(one,all,4,MPI_DOUBLE,MPI_SUM,world);
```

能量总和在`标量`变量（*第 325 行*）中计算，如果需要将静电能量合并到`kspace_contribution()`方法中（*第 346 行* 和*第 363 行*），则更新能量总和。每个力分量的总和存储在 1D `vector` 数组中（*第 326 行*）。`scalar` 和 `vector` 都继承自 `compute.h` 类。

此方法执行的计算由 `compute_scalar()` 和 `compute_vector()` 方法调用，根据计算的要求。

### ComputeGroupGroup::compute_scalar() 和 compute_vector()
`compute_scalar()` 方法调用 `pair_contribution()` 和 `kspace_contribution()` 方法来获得对和静电能量。`compute_vector()` 方法调用相同的方法来获取力分量。以下屏幕截图中描述了这两种方法：

<div align=center>
<img src=ch06/fig6_9.jpg>
</div>

图 6.9 – 来自 compute_group_group.cpp 的代码片段，显示了 compute_scalar() 和 compute_vector() 方法

在本节中，我们分析了一些构成`compute group/group`并返回两组原子之间相互作用的能量和力的方法。

在下一节中，我们将了解`compute RDF`。

## 查看compute RDF 类
在本节中，我们将查看 `compute_rdf.cpp` 的源代码和
`compute_rdf.h` 类，它管理compute RDF 类。

`compute RDF` 类计算一组原子的 **径向分布函数 (RDF)** 并返回归一化的邻居数，分类到具有统一径向宽度 ($\Delta {r}$ ) 的 bin 中，并且从零延伸到指定的截止距离。输出作为一维数组返回，可用于绘制 RDF 直方图。对于每个原子，它在距离 r 处的 RDF ($g(r)$) 的值是通过取其位于距离范围 $ 内的相邻原子 ($N$) 的数量的比率来计算的$(r, r + \Delta {r})$ 从原子，到如果原子分布均匀，将位于同一范围内的原子数。

为了推导出这个函数，我们必须将中心原子的外围空间划分为宽度为 dr 的壳，每个壳占据 $dV_{shell}$ 的体积：

$$dn(r)=\rho g(r)dV_{shell}$$

从前面的等式中，$dn(r)$ 表示半径为 $r$ 和宽度为 $dr$ 的壳中的邻居数，而 $dV_{shell}$ 是该壳的体积，在 3D 中近似为以下等式：

$$dV_{shell}=\frac{4}{3}\pi (r+dr)^3-\frac{4}{3}\pi r^3\simeq 4\pi r^2dr$$

数量 $\rho$ 表示系统的原子密度，由以下等式计算：

$$\rho=\frac{Number\ of\ Neighboring\ Atoms\ in\ Simulation\ Box}{Volume\ of\ Simulation\ Box}$$

因此，可以使用以下等式求解 RDF：

$$g(r)=\frac{dn(r)}{4\pi r^2dr\rho}$$

总之，当实现离散的 bin 宽度 $\Delta {r}$ 时，3D 中的 RDF 转换为以下内容：

$$g(r)=\frac{N(r, r+\Delta {r}}{4\pi r^2(\Delta {r})\sigma}=\frac{N(r, r+ \Delta { r})}{\frac{4\pi}{3}[(r+\Delta {r})^3-r^3]\rho}$$

系统 $g(r)$ 是通过对系统中所有原子的单个 $g(r)$ 值进行平均来计算的。

在 2D 中，$g(r)$ 的等式在分母中发生变化：

$$g(r)=\frac{N(r, r+\Delta {r}}{2\pi r^2(\Delta {r})\sigma}=\frac{N(r, r+ \Delta { r})}{\pi [(r+\Delta {r})^3-r^3]\rho}$$

系统密度 $\sigma$ 是在单位面积上计算的：

$$\sigma =\frac{number\ of\ neighboring\ atoms\ in\ simulation\ box}{area\  of\  simulation\  box}$$

用于实现计算 RDF 的 LAMMPS 输入脚本命令如下：

```
compute COMPUTE_ID GROUP rdf NBIN
```

如您所见，`COMPUTE_ID`参数是compute的唯一 ID，而`GROUP`参数是compute所作用的原子组的 ID（参见手册：https://lammps.sandia.gov/doc/compute_rdf.html)。`NBIN` 整数是要compute $g(r)$ 的 bin 数量。可以在这些参数之后输入手册中列出的可选关键字，以指定原子类型对（多次）和用户定义的截止值。构造函数方法容纳了这些可选参数，我们将在以下部分中看到。

### ComputeRDF::ComputeRDF()
以下屏幕截图显示了 `ComputeRDF()` 构造方法：

<div align=center>
<img src=ch06/fig6_10.jpg>
</div>

图 6.10 – 来自 compute_rdf.cpp 的代码片段显示了构造方法

如您所见，该方法接受至少四个输入参数（*第 47 行*）并激活一个数组标志以生成输出（*第 49 行* 到 *50*）。在 *line 52* 中读取 bin 的数量并存储为整数 *nbin*。

为了计算输入的任何其他参数，必须引入 `nargpair` 变量并在第四个参数处重置为零（*第 64 行*）。在循环附加参数时（*第 66 行* 到 *74* 行），任何用户定义的截止都使用 `strcmp()` 方法（*第 67 行*）来识别，该方法定位 `cutoff` 关键字。`force->numeric()` 方法（*第 69 行*）读取此关键字后面的相应数值。如果没有提供截止值，则默认使用系统中定义的最长的成对截止值（参见 `init()` 方法）。

在 `compute_rdf.h` (*line 38*) 中引入了 `npairs` 变量，以记录可选输入参数指定的原子类型对的数量。在构造函数中，如果没有输入额外的原子类型对（*第 78 行*），则为 `npairs` 赋予默认值 `1`。对于输入的任何其他原子类型对，`npairs` 变量等于输入的对数（*第 81 行*），同时确保输入的原子类型为偶数（*第 80 行*）。

*第 84* 和 *85* 行定义了 2D 全局数组的大小（继承自 `compute.h`），它拥有 `nbin` 行数和 (`1+2*npairs`) 列以容纳 bin 坐标、$g(r)$ 值和协调数（在 LAMMPS 手册中描述）以提供输出。

*第 88* 到 *89* 行初始化了两个新结构——一个名为 `rdfpair` 的 3D 数组和一个名为 `nrdpair` 的 2D 数组。输入的任何可选原子类型对都会在 *line 98* 到 *107* 中解析。为了处理输入脚本中输入的通配符原子类型对，使用了`force->bounds`方法（*lines 101* to *102*），它可以将原子类型的下限和上限分配给` ilo`、`ihi`、`jlo` 和 `jhi` 变量。例如，如果系统中有五种原子类型，输入的pair类型是`* 4`，那么我们得到`ilo=1`、`ihi=4`和`jlo=jhi=4`。

在*第 109* 到 *120* 行中，解析的原子类型对存储在 `nrdfpair` 和 `rdfpair` 数组中，我们之前在第 88 到 89 行中声明了它们。二维的 `nrdfpair` 数组记录了原子类型的数量应该计算 RDF 的对，而 3D rdfpair 数组添加了一个额外的维度来指示原子类型对的索引。数组本身记录了对应于已解析的原子类型对的`npair`值。这两个结构有助于防止重复计算并节省 `compute_array()` 方法中的计算时间。

在构造函数方法的其余部分中，为 2D 全局数组 (`array`) 分配内存（*第 122* 到 *124* 行）以生成最终输出，并为 2D `hist` 和 `histall` 数组分配内存以容纳$g(r)$ 值，直到输出准备好。这些数组与 `compute_rdf.h` 中定义的任何其他数组一起在析构函数方法（*第 144* 到 *146* 行）中被释放。

接下来，我们将看看 `init()` 和 `init_list()` 方法。

### ComputeRDF::init() 和 init_list()
`init()` 方法计算截止值和趋肤宽度，同时根据该截止值确定 bin 宽度，如以下屏幕截图所示：

<div align=center>
<img src=ch06/fig6_11.jpg>
</div>

图 6.11 – 来自 compute_rdf.cpp 的代码片段显示了 init() 方法

如您所见，*第 162* 到 *180* 行计算在提供用户定义的截止时可应用于鬼原子的截止，然后执行检查以降低重建邻居列表的频率（*第 175* 行到*177*)。然后，通过将截止除以`nbin`来获得 bin 宽度`delr`。如果未提供用户定义的截止（*第 181 行*），则将截止设置为在所有已实现的对样式中定义的最大截止距离（通过 `force->pair->cutforce` 访问）。倒数 bin 宽度 (`delrinv`) 在*第 183 行* 中计算，以减少与除法运算符相关的任何计算开销。

在 *187* 到 *188* 行中，`array` 的第一列填充了 bin 的中点。在此方法结束时请求邻居列表（*第 205* 到 *212* 行）。`init_list()` 方法以与`compute group/group` 类中发生的方式相同的方式提供对指向邻居列表的指针的访问。

### ComputeRDF::init_norm()
`init_norm()` 方法计算 $g(r)$ 计算中使用的归一化因子。这里，确定组中存在的每种类型的原子数，并在计算 RDF 时用于计算归一化因子。以下屏幕截图显示了此方法：

<div align=center>
<img src=ch06/fig6_12.jpg>
</div>

图 6.12 – 来自 compute_rdf.cpp 的代码片段显示了 init_norm() 方法

计算组中（每种类型的）原子数（*第 235* 到 *237* 行）。在*第 243* 到 *252* 行中，确定了用于计算 RDF 的原子数（每种类型）。如果使用原子类型对的默认数量（`npairs = 1`），则在 RDF 计算中使用索引、`m = 0` 和系统中的所有原子类型。

索引 `i` (*line 245*) 和索引 `j` (*line 247*) 循环遍历所有原子类型，以计算要在计算中使用的原子数，然后将它们存储为变量；即，`icount[m]` 和 `jcount[m]`（参见构造方法中的 *lines 95* 到 *107*）。出现在两个组中的重复原子的数量也计为`duplicates[m]`。如果输入了可选的原子类型，则索引范围`i`和`j`仅涵盖指定的原子类型，并且`icount[m]`和`jcount[m]`等于这些类型的原子数只要。

`icount`、`jcount` 和 `duplicates` 数组通过 *254* 到 *261* 行中的 `MPI_Allreduce()` 方法更新为包含来自所有内核的原子。数组元素在 `compute_array()` 方法中用于计算和规范化 RDF，如下一节所述。

### ComputeRDF::compute_array()
`compute_array()` 方法计算指定组和原子类型的 RDF 和配位数。首先，邻居列表用于遍历核心中的所有原子并计算每个中心原子与其每个邻居之间的距离 r（*第 347 行*）。此 r 的 bin 索引 (`ibin`) 是通过舍入到小于 r 与 *line 348* 中的倒数 bin 宽度 (`delrinv`) 的乘积的最大整数来确定的，如下所示：

```
ibin = static_cast<int> (r*delrinv);
```

2D `hist` 数组中的相应 bin 增加 1 以计算此近邻
（*第 353 行*）：

```
hist[m][ibin] += 1.0;
```

在 `hist` 数组中，第一个索引 (`m`) 标识了计算邻居的原子类型对，而第二个索引 (`ibin`) 表示计算邻居的 bin。hist 数组对所有内核求和（*第 366 行*），结果存储在 2D `histall` 数组中。

这样， histall 数组包含了每个 bin 中的原子数，以及整个模拟框中所有中心原子的总和；也就是说， $\sum_{i=1} ^{icount[m]} {N_i(r, r+\Delta {r})}$ ，按照指定的原子类型对被排序到单独的行中。如果没有提供额外的原子类型对，`histall` 由对应于所有原子类型的单行组成。

RDF 计算的最后一部分涉及适当的归一化，如以下屏幕截图所示：

<div align=center>
<img src=ch06/fig6_13.jpg>
</div>

图 6.13 – 来自 compute_rdf.cpp 的代码片段显示了 compute_array() 方法


如您所见，3D 系统的 $g(r)$ 在 *375* 到 *394* 行中计算。数量 $\frac{4\pi}{3V}$ 在 *line 376* 中定义为常量变量，其中 V 是模拟框的体积。然后，对于每个原子类型对 (`m`)，在 *line 379* 到 *380* 中计算归一化因子 (`normfac`)，方法是从所有 bin 的邻居计数 (`jcount[m]`)。

然后代码循环遍历每个 bin（*第 382* 到 *393* 行）并计算以下数量：

- `rlower`: bin 的下限 ($r$)

- `rupper`: bin 的上限 ($r+\Delta {r}$)

- `vfrac`: bin $\frac{4\pi}{3V}[(r+\Delta {r})^3-r^2]$所占的体积分数

然后在第 387 行计算系统的 $g(r)$，如下所示：

$$g(r)=\frac{\sum_{i=1}^{icount[m]}{N_i(r,r+\Delta{r})}}{\frac{4\pi}{3V}[ (r+\Delta {r})^3-r^3]*normfac*icount[m]}$$

$\frac{normfac}{V}$ 比率表示相邻原子的密度 ($\rho$)，而 $\frac{\sum_{i=1}^{icount[m]}{N_i(r , r+\Delta {r})}}{icount[m]}$ 表示所有参与原子的平均值 $N(r, r+\Delta {r})$。

因此，这个表达式表示系统在 3D 中的平均 $g(r)$：

$$g(r)=\frac{N_{avg}(r, r+\Delta {r})}{\frac{4\pi}{3}[(r+\Delta {r})^3-r^ 3]\rho}$$

$g(r)$ 也可以在 *line 396* 到 *416* 中为 2D 系统计算。现在，`constant` 变量（*第 397 行*）被定义为 $\frac{\pi}{A}$，其中 A 是模拟框的面积，随后是 `vfrac` 变量（*第 406 行* ) 等于 $\frac{\pi}{A}[(r+\Delta {r})^2-r^2]$。这里，系统的$g(r)$如下（*第408行*）：

$$g(r)=\frac{\sum_{i=1}^{icount[m]}{N_i(r, r+\Delta {r})}}{\frac{\pi}{A}[( r+\Delta {r})^2-r^2]*normfac*iconut[m]}$$

与 3D 情况类似，$\frac{normfac}{A}$ 比率表示相邻原子的密度（$\sigma$），系统在 2D 中的平均 $g(r)$ 表达式为通过使用以下等式获得：

$$g(r)=\frac{N_{avg}(r,r+\Delta {r})}{\pi[(r+\Delta {r})^2-r^2]\sigma}$$

对于 3D 和 2D 系统，系统的平均 $g(r)$ 通过循环遍历每个 bin 来计算不同的 $r$ 值。使用 $g(r)$ 的值计算相应的配位数（*第 390* 和 *411* 行）。

结果存储在 `array` 的对应列中（*lines 391* to *392*, and *lines 412* to *413*）以生成LAMMPS手册中描述的格式的输出（即bin中点在第一列中，$g(r)$ 值在下一列中，配位数在第三列中，并且对于任何其他原子类型对，在连续列中重复 $g(r)$ 和配位数）。

应该注意的是，`compute RDF` 计算单个时间步的模拟快照的平均 $g(r)$ 和协调数。要在所需的迭代次数上找到这些量的时间平均值，需要在 LAMMPS 输入脚本中使用`fix ave/time`命令，如 LAMMPS 手册中所述。

在本节中，我们介绍了`计算 RDF`，它广泛使用 2D 数组通过数值的行和列生成 2D 输出。在下一节中，我们将查看`计算热通量`作为我们的最后一个示例。

## 查看compute heat flux类

在本节中，我们将研究计算`heat flux`类的源代码，该类包含在`compute_heat_flux.cpp`和`compute_heat_flux.h`中。

`compute heat flux`类接受每个原子的动能、每个原子的势能和每个原子的应力来计算热流 (J)。这可以计算如下：

$$J=\frac{1}{V}[\sum_{i}{e_i\vec {v_i}}-\sum_{i}{S_i\vec {v_i}}]$$

在这个方程中，$e_i$ 代表原子 i 的动能和势能之和，$\vec {v_i}$ 代表速度向量 $(v_{xi},v_{yi},v_{zi})$原子`i`，$S_i$ 表示原子`i`的应力张量，`V`表示考虑的原子所占的体积。$e_i\vec {v_i}$ 的总和是热通量的对流部分，而 $S_i\vec {v_i}$ 的总和是热通量的维里部分。

因此，该计算需要动态读取原子的动能、势能和应力，这可以通过将这些量作为其他计算输入来促进。因此，有效地，`compute heat flux`接受其他三个计算作为输入参数，每次迭代都会评估这些计算以更新所需的数量。

原子的应力张量 (`S`) 用其分量 $S_{ij}$ 表示，如下所示：

$$S=\begin{bmatrix}
S_{xx}&S_{xy}&S_{xz} \\
S_{yx}&S_{yy}&S_{yz} \\
S_{zx}&S_{zy}&S_{zz} \end{bmatrix}$$

因此，张量积 $S\vec {v}$ 计算如下：

$$
S\vec {v}=
\begin{bmatrix} 
S_{xx}&S_{xy}&S_{xz} \\
S_{yx}&S_{yy}&S_{yz} \\
S_{zx}&S_{zy}&S_{zz} \end{bmatrix}
\begin{bmatrix} 
V_x \\
V_y \\
V_z
\end{bmatrix}
=\begin{bmatrix} 
S_{xx}V_x&S_{xy}V_y&S_{xz}V_z \\
S_{yx}V_x&S_{yy}V_y&S_{yz}V_z \\
S_{zx}V_x&S_{zy}V_y&S_{zz}V_z 
\end{bmatrix}$$

张量 S 沿其对角线对称（即 $S_{xy}=S_{yx}, S_{xz}=S_{zx}, S_{yz}=S_{zy}$ `compute stress/atom`，虽然它可以是不对称的，如果由`compute centroid/stress/atom`返回。在对称`S`的情况下，有六个独特的元素作为六元素存储在LAMMPS中一维数组（$S_{xx}、S_{yy}、S_{zz}、S_{xy}、S_{xz}、S_{yz}$）。

在非对称`S`的情况下，分量存储在一个九元素一维数组中（$S_{xx}, S_{yy}, S_{zz}, S_{xy}, S_{xz}, S_{ yz}，S_{yx}，S_{zx}，S_{zy}$）。在 `commpute_heat_flux.cpp` 中，当计算热流 $J$ 中的张量积 $S\vec {v}$ 时，从这些数组中访问组件，如本节所示。

用于实现`计算热通量`的 LAMMPS 输入脚本命令是
如下：

```
compute COMPUTE_ID GROUP heat/flux KE PE STRESS
```

在上述代码中，`COMPUTE_ID` 参数是计算的唯一 ID，而 `GROUP` 参数是计算所作用的原子组的 ID（参见手册：https://lammps.sandia.gov/doc/compute_heat_flux.html)。`KE`、`PE`和`STRESS`计算分别计算用于计算热通量的每原子动能、每原子势能和每原子应力。这些计算应该适用于`GROUP`中的原子。

在头文件 `compute_heat_flux.h` 中，声明了一维字符数组以解析输入的计算名称（*第 35 行*）：

```
 char *id_ke,*id_pe,*id_stress;
```

同样，声明 `compute` 对象来评估输入的计算（*第 36 行*）：

```
 class Compute *c_ke,*c_pe,*c_stress;
```

在`compute_heat_flux.cpp`中，我们首先分析解析这些计算的构造函数。

### ComputeHeatFlux::ComputeHeatFlux()
以下屏幕截图显示了 ComputeHeatFlux() 构造函数方法：

<div align=center>
<img src=ch06/fig6_14.jpg>
</div>

图 6.14 – 来自 compute_heat_flux.cpp 的代码片段，显示了构造方法

如您所见，解析了动能（*第 47* 到 *49* 行）、势能（*第 51* 到 *53* 行）和应力（*第 55* 到 *57* 行）的三个计算名称作为字符数组；即分别为`id_ke[]`、`id_pe[]`和`id_stress[]`。对应的计算 ID 通过 `modify->find_compute()` 方法在 *line 59* 到 *61* 中找到，如果计算 ID 不存在则返回错误（*lines 62* 到 *63*）。然后，*第 64* 到*71* 行检查输入的计算是否专门计算每个原子的动能、每个原子的势能和每个原子的应力。在关闭此方法之前，会创建一个长度为 6 的输出向量（*第 73 行*）。

进一步的错误检查在 `init()` 方法中执行，如下所述。

### ComputeHeatFlux::init()
`init()` 方法，如以下屏幕截图所示，检查在构造函数中解析的计算名称是否存在（*第 92* 到 *96* 行）：

<div align=center>
<img src=ch06/fig6_15.jpg>
</div>

图 6.15 – 来自 compute_heat_flux.cpp 的代码片段显示了 init() 方法

如您所见，计算对象——即`c_ke`、`c_pe`和`c_stress`（*第 98 行*到 *100*）——被标记为计算每个原子动能的计算对象，每分别为原子势能和每原子应力。

接下来，在 `compute_vector()` 方法中计算热通量。

### ComputeHeatFlux::compute_vector()
在 `compute_vector()` 方法中，如果尚未调用计算（*第 111* 到 *122* 行），则调用计算，以及 `ke[]`、`pe[]` 和 `stress[] []` 数组用于分别提取每个原子的动能、势能和应力所需的值（*第 129 行* 到 *131* 行）。

在*第 142 行*，应力张量的对称性，`S`，由​​`pressatomflag` 标志确定，对于由`compute stress/atom`计算的对称张量，该标志设置为`1`（参见`compute_stress_atom. cpp`）或 `2` 用于由`compute centroid/stress/atom`计算的不对称张量（请参阅`compute_centroid_stress_atom.cpp`）。基于这种对称性，原子 i 的张量积 $S_i\vec {v}_i$ 的计算方式不同，如以下屏幕截图所示：

<div align=center>
<img src=ch06/fig6_16.jpg>
</div>

图 6.16 – 来自 compute_heat_flux.cpp 的代码片段显示了 compute_vector() 方法

如您所见，非对称 S 张量的情况（其中 `pressatomflag` 为 `2`）在 *143* 到 *168* 行中计算。在遍历核心中的所有原子（*第 143 行*）并选择属于指定组的原子（*第 144 行*）时，原子`i`的动能和势能之和（$e_i$）计算为eng 变量（*第 145 行*）。在 *146* 到 *148* 行，`jc[]` 数组计算对流部分的三个分量：

- 对于 x 分量，我们使用以下等式：
	$$jc[0]=\sum_{i}{e_iv_{xi}}$$

- 对于 y 分量，我们使用以下等式：
	$$jc[1]=\sum_{i}{e_iv_{yi}}$$

- 对于 z 分量，我们使用以下等式：
	$$jc[2]=\sum_{i}{e_iv_{zi}}$$

热通量的维里部分的分量被计算为`jv[]`数组，它是根据每个原子的应力计算的。在不对称`S`的情况下，`stress[][]`的数组元素访问原子`i`的应力的以下分量：

- `stress[i][0]` = $S_{xx}$
- `stress[i][1]` = $S_{yy}$
- `stress[i][2]` = $S_{zz}$
- `stress[i][3]` = $S_{xy}$
- `stress[i][4]` = $S_{xz}$
- `stress[i][5]` = $S_{yz}$
- `stress[i][6]` = $S_{yx}$
- `stress[i][7]` = $S_{zx}$
- `stress[i][8]` = $S_{zy}$

使用前面的应力分量，计算维里部分的三个分量，`jv[]`（*第 159* 到 *166* 行）：

- 对于 x 分量，我们使用以下等式：
	$$jv[0]=-\sum_{i}{s_{xxi}v_{xi}+s_{xyi}v_{yi}+s_{xzi}v_{zi}}$$
- 对于 y 分量，我们使用以下等式：
	$$jv[1]=-\sum_{i}{s_{yxi}v_{xi}+s_{yyi}v_{yi}+s_{yzi}v_{zi}}$$
- 对于 z 分量，我们使用以下等式：
	$$jv[2]=-\sum_{i}{s_{zxi}v_{xi}+s_{zyi}v_{yi}+s_{zzi}v_{zi}}$$

然后，在对称张量`S`的情况下，热通量的对流部分`jc[]`以与之前相同的方式计算（*第 172* 到 *175* 行）和维里部分, `jv[]`, 使用 `stress[i][0]` 到 `stress[i][5]` 元素计算，这些元素访问应力分量 $(S_{xx},S_{yy},S_ {zz},S_{xy},S_{xz},S_{yz})$ 以给定的顺序。因此，`jv[]` 的组件如下所示（*第 176* 到 *181* 行）：

- 对于 x 分量，我们使用以下等式：
	$$jv[0]=-\sum_{i}{s_{xxi}v_{xi}+s_{xyi}v_{yi}+s_{xzi}v_{zi}}$$

- 对于 y 分量，我们使用以下等式：
	$$jv[1]=-\sum_{i}{s_{xyi}v_{xi}+s_{yyi}v_{yi}+s_{yzi}v_{zi}}$$

- 对于 z 分量，我们使用以下等式：
	$$jv[2]=-\sum_{i}{s_{xzi}v_{xi}+s_{yzi}v_{yi}+s_{zzi}v_{zi}}$$

一旦计算了对流部分和维里部分，维里部分`jv[]`通过适当的压力单位缩放转换为与对流部分`jc[]`相同的单位（*第188行*至 *191*)。最后，使用元素创建 `data[]` 数组（*第 197 行*），如以下代码所示：

```
 double data[6] = 
{jc[0]+jv[0],jc[1]+jv[1],jc[2]+jv[2],jc[0],jc[1],jc[2]};
```

从前面的代码中，所有内核上这六个元素的总和然后通过`MPI_A llreduce()`调用（*第 198 行*）检索并作为在 `compute.h` 中声明的一维输出`vector`数组返回父类。

应该注意的是，虽然体积`V`包含在计算热通量`J`的方程中，但体积不会出现在计算`热通量`的源代码中。因此，必须单独合并体积（例如，在单独的计算或变量命令中）以计算正确的热通量。

本节到此结束，我们讨论了一个接受其他计算作为输入参数来计算热通量的计算。

## 概括

在本章中，我们检查了计算全局量、识别组之间的交互、计算范围内的距离和执行张量乘法的计算。这些计算说明了可以通过计算样式源代码中的适当实现来执行的处理和计算范围。

通过阅读本章，您将了解如果正确定义，计算可以计算并返回标量和向量，并且可以通过使用指针请求和访问相邻列表在计算中实现邻居列表。此外，`MPI_Allreduce()`方法有助于组合来自多个内核的输出。

在下一章中，我们将分析在模拟运行期间可以在模拟系统上执行各种操作的选定修复样式。

## 问题

1. 哪种方法负责解析超出 LAMMPS 输入脚本中计算命令中输入的第三个参数？

2. `MPI_Allreduce()` 方法接受哪些参数？

3. `destructor` 方法的主要目的是什么？

4. 在 `compute reduce` 中，`replace` 选项可用于查找数量的索引，以及输入数量的最大值或最小值（参见 https://lammps.sandia.gov/doc/compute_reduce.html ）。根据这些信息回答下列问题：

	在 `compute_reduce.cpp` 中，哪些行负责识别所有处理器中具有最大值的数量的索引？

	在 `compute_reduce.cpp` 中，哪些行负责识别所有处理器中具有最小值的数量的索引？