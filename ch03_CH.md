# 源代码结构和执行阶段

接续上一章，本章将通过描述源代码中的父类和子类来进一步说明源代码层次结构。还将讨论其中一些类中使用的各种方法。为了完成这幅图，将使用集成商确定的方法执行顺序来解释从顶级类到代码终止阶段的控制流。

本章旨在完成您对源代码层次结构和流程的理解，这将帮助您识别在实现自定义功能时需要修改的代码部分。

我们将涵盖以下主题：


- 源代码中的父类和子类

- 每个时间步内执行模拟的阶段

- `pointers.h` 所扮演的角色

- `input.cpp` 解析输入脚本命令

## 技术要求
要执行本章中的指令，您只需要一个文本编辑器（例如，**Notepad++** 或 **Gedit**）。

您可以在此处找到本章中使用的完整源代码：https://github.com/PacktPublishing/Extending-and-Modifying-LAMMPS-Writing-Your-Own-Source-Code

这是下载 LAMMPS 的链接：https： //lammps.sandia.gov/doc/Install.html。这是 LAMMPS GitHub 链接 (https://github.com/lammps/lammps)，也可以在其中找到源代码。

## 引入父类和子类

如*第 2 章，LAMMPS 语法和源代码层次结构*中所述，源代码中给出了可以支持子类的某些样式。这些样式作为父类，它们的子类继承它们的方法，从而确保子类的一致性，从而使分类和语法合成更加流线化。在本节中，我们将描述其中的一些父类及其子类中的一些继承方法。

### fix.cpp 和 fix.h
这些是 LAMMPS 中使用的所有修复程序的父类。除其他目的外，他们读取所有修复共有的前三个参数（修复 ID、组 ID 和修复样式）并设置能量或维里计算。以下屏幕截图显示了来自 `fix.cpp` 的代码片段，它调用 LAMMPS 实例并读取三个常用参数：

<div align=center> 
<img src=./ch03/fig3_1.jpg> 
</div>

图 3.1 - 来自的代码片段fix.cpp

所有修复程序都从 `fix.cpp` 继承这三个参数，以及设置计算的方法。从 `fix.cpp` 父类继承的典型修复子类包含以下行（从 `fix_addforce.cpp` 中提取）：

<div align=center> 
<img src=./ch03/fig3_2.jpg> 
</div>

图 3.2 – 从 fix.cpp 继承的 fixaddforce.cpp

这个子修复类从其父类继承三个参数，然后读取特定于自身的参数。同样，`fix.h` 存储了所有修复子类共享的变量列表，如下所示：

<div align=center> 
<img src=./ch03/fig3_3.jpg> 
</div>图 3.3 - fix.h 
Child

中的代码片段
修复头类通过包含 `fix.h` 继承这些变量：

``` 
#include "fix.h" 
```

这样，所有子修复共享通用语法和变量，例如 id、igroup 和样式来存储修复分别是 ID、组 ID 和修复样式。

例如，在 LAMMPS 输入脚本中输入了一个修复命令：

```
fix F1 all nvt temp 100 100 0.1 
```
数量的解析描述如下：

- `fix`：这是命令的类型，由 input.cpp 注册。

- `F1`：这是修复的名称，由`fix.cpp`解析。

- `all`：这是组名，由 f`ix.cpp` 解析。

- `nvt`：这是修复的样式，由`fix.cpp`解析。

- `temp, 100, 100, 0.1`：这些是特定于修复的关键字，由修复子类解析。


### pair.cpp 和 pair.h
这些父类通过验证是否已指定所有对系数、在需要时应用混合规则、确定截止值、请求邻居列表和设置计算来初始化所有对样式。以下屏幕截图显示了 pair.cpp 中执行这些任务的一些代码段。首先，让我们看看配对系数分配：

<div align=center> 
<img src=./ch03/fig3_4.jpg> 
</div>

图 3.4 - pair.cpp 中显示配对系数分配代码的代码片段

接下来，这里是混合系数的屏幕截图: 

<div align=center> 
<img src=./ch03/fig3_5.jpg> 
</div>

图 3.5 - 来自 pair.cpp 的代码片段显示了混合系数的代码

最后，我们将看到请求默认邻居列表的代码：

<div align=center> 
<img src=./ch03/fig3_6.jpg> 
</div>

图 3.6 - pair.cpp 中的代码片段显示了请求默认邻居列表的

代码包括pair.h，子对类可以继承用于实现原子间对势和计算相关量的变量、类和方法，包括：

- `cutforce`：这是所有原子对的最大截止距离。

- `cutsq`：这是存储每个原子对的截止平方的指针。

- `list`：这是大多数对潜在计算中使用的邻居列表的指针。

- `compute`：这是计算原子对之间力的方法。

- `single`：这是计算对相互作用能量的方法。

在输入脚本中的“pair_style”和“pair_coeff”命令之后，对系数和截止值是从子对类而不是父类中读取的。

当需要使用 `pair_style hybrid` 命令实现多种类型的势时，`pair_hybrid.cpp` 为每个指定的样式创建不同的对样式，并将每个原子对映射到正确的对样式，如下面的屏幕截图所示。首先，我们将查看显示配对样式创建的片段：

<div align=center> 
<img src=./ch03/fig3_7.jpg> 
</div>

图 3.7 - pair_hybrid.cpp 中显示配对样式创建的代码片段

接下来，我们将看到pair样式映射：

<div align=center> 
<img src=./ch03/fig3_8.jpg> 
</div>

图 3.8 - pair_hybrid.cpp 中显示对样式映射的代码片段

当使用 `pair_style hybrid/overlay` 命令叠加对电位时，`pair_hybrid_overlay。 cpp`执行类似的任务以及映射相同原子对的多对电位的选项。

### compute.cpp 和 compute.h
所有计算的父类读取前三个参数，类似于修复类，并允许子类读取剩余的参数，如下图所示：

<div align=center > 
<img src=./ch03/fig3_9.jpg> 
</div>

图 3.9 - 来自 compute.cpp 的代码片段

这样，所有子计算共享 `id`、`igroup` 和 `style` 变量来分别存储修复 ID、组 ID 和计算样式。同样，子类从 `compute.h` 继承变量，包括 `scalar` 和 `vector`，它们用作典型计算的输出。

对应于一个类的每个 `.cpp` 文件都包含几个方法，这些方法在模拟的不同阶段执行不同的功能。以正确顺序执行方法的序列由集成器执行，例如 `verlet.cpp` 中的 Verlet 类。Verlet 类是integrate.cpp 的子类，它又由顶级类`update.cpp` 发起。在下一节中，我们将描述由 `verlet.cpp` 确定的流控制。

## 执行模拟的阶段

LAMMPS 模拟通过迭代时间步长（例如，**velocity Verlet 积分**）或通过不执行时间步长的算法（例如，**最小化**）来执行。接下来，我们将描述在 verlet.cpp 中实现的 Verlet 集成方案，然后简要概述在 `min.cpp` 中实现的最小化方案。


### verlet.cpp 
`verlet.cpp` 类通过一系列按预定义顺序执行的方法实现时间步进。在时间步的开始，`verlet.cpp` 中的以下方法被调用：

- `init()`：该方法检查输入脚本中是否定义了修复并为数组设置标志，如下所示：

<div align=center> 
<img src=./ch03/fig3_10.jpg> 
</div>

图 3.10 – verlet.ccp 中显示 init() 方法的代码片段

- `force_clear()`：此方法清除所有原子上的力以在时间步长过程中存储合力：

<div align=center> 
<img src=./ch03/ fig3_11.jpg> 
</div>

图 3.11 - verlet.ccp 中显示 force_clear() 方法的代码片段

- `setup()`：该方法使用幽灵原子设置域并构建邻居列表；然后根据需要计算力以在速度 Verlet 算法中执行位置和速度更新：

<div align=center> 
<img src=./ch03/fig3_12.jpg> 
</div>

图 3.12 - verlet.ccp 中显示 setup() 方法的代码片段

时间步长中的其余步骤由 `run()` 方法执行。此方法以预定义的顺序调用其他序列控制方法，这些方法又在同一时间步内的相应点处调用修复。`run()` 方法的以下屏幕截图显示了执行顺序：

<div align=center> 
<img src=./ch03/fig3_13.jpg> 
</div>

图 3.13 - 显示来自 verlet.cpp 的 run() 的代码片段

这是一个表格，按时间顺序向我们展示了其中一些方法的列表：

<div align=center> 
<img src=./ch03/table3_1.jpg> 
</div>

表 3.1 - 显示方法列表的表格

每个修复都被分配了一个序列控制方法来确定其执行顺序，如*第 7 章，了解修复*中所述。`modify.cpp` 类存储所有修复并在时间步的适当点调用这些方法，从而促进以指定顺序执行所有修复。方法列表从 `modify.h` 导入到 `verlet.cpp`。`modify.h` 的屏幕截图按执行顺序显示了这些方法的列表：

<div align=center> 
<img src=./ch03/fig3_14.jpg> 
</div>

图 3.14 - 显示 modify.h 中的方法的代码片段
编写自定义修复程序时，可以通过在修复程序中合并上述一种或多种方法来指定执行修复程序的确切位置。

### min.cpp
`min.cpp` 类也由 update.cpp 启动，它使用非时间步长算法来执行最小化。类似于时间步长算法，执行力计算、邻居列表构建和对处理器的原子指定。但是，调用了不同的修复类序列控制方法：`min_pre_exchange()`、`min_pre_force()` 和 `min_post_force()`。可以在最小化期间​​调用包含这些方法的修复程序。

到目前为止，我们已经介绍了源代码层次结构和在时间步或非时间步过程中规定的控制流，这些过程在源代码中调用不同的类（请参阅 https://lammps.sandia.gov/doc/Developer 上的 LAMMPS 手册.html 了解更多信息）。这些类能够通过“Pointers”类（“pointers.h”）共享信息，该类包含指向 lammps.h 中列出的指针的指针，如下所述。

## 指针类

的作用 指针类的作用 `pointers.h` 文件通过创建指向 lammps.h 中列出的所有重要量的指针来促进类之间的信息传输。所有类都继承自 pointers.h 并且能够访问这些变量。以下屏幕截图说明了创建的指针：

<div align=center>
<img src=./ch03/fig3_15.jpg> 
</div>

图 3.15 - 来自 pointers.h 的代码片段显示（左）指向类的指针和（右）指向 lammps.h 中
 
的指针的指针 创建指向不同类的指针，以及指向 `lammps.h` 中列出的指针的指针，由前面的 `*&` 指示。这样，可以通过声明正确的指针直接从其他类访问来自“lammps.h”的变量，这将在后面的章节中解释。

在本章的下一节中，我们将描述如何使用 `input.cpp` 中预定义的允许命令列表来解析输入脚本命令。

## 通过 input.cpp 解析输入脚本命令

在本节中，输入脚本命令的解析被描述为由 `input.cpp` 中的 `execute_command()` 方法处理，以及每个命令之后的步骤。

input.cpp 中的 execute_command() 方法负责解析输入脚本每一行的第一个单词。此方法包含与每行的第一个单词进行比较的允许命令列表。如果没有匹配，则返回错误，并且为每个匹配调用 `input.cpp` 中的预定义方法。该方法在 `file()` 方法和 `input.cpp` 中的 `one()` 方法中调用，以方便解析和执行。

以下屏幕截图显示了 `execute_command()` 方法：

<div align=center> 
<img src=./ch03/fig3_16.jpg> 
</div>

图 3.16 – input.cpp 中的 execute_command() 方法包含允许的输入脚本命令列表

如您所见，在前面的屏幕截图中，command 变量表示正在解析的行的第一个单词，并将其与列表进行比较允许的命令。例如，如果匹配 `clear` 这个词（*第 783 行*），则会调用 `input.cpp` 中的 `clear()` 方法，如下图所示：

<div align=center> 
<img src=./ch03/fig3_17.jpg> 
</div>

图 3.17 – input.cpp 中的 clear() 方法

正如你所看到的，`clear()` 方法是自包含的，它有效地删除了以前的 LAMMPS 实例并呈现了一个全新的平台接着说。

类似地，如果匹配 `lattice` 这个词（*第 824 行*），则调用 `input.cpp` 中的 `lattice()` 方法，如下图所示：

<div align=center> 
<img src=./ch03/ fig3_18.jpg> 
</div>

图 3.18 – input.cpp 中的 lattice() 方法

如您所见，在 `lattice()` 方法中，调用了 `domain->set_lattice()` 方法，打开了 ` domain.cpp` 用于进一步执行。

为其他命令调用方法，包括 `fix` 和 `pair_style`，并且经常调用外部类来继续模拟。在输入脚本中添加新的允许命令将需要修改此方法以反映这些更改。

# 概括

在本章中，介绍了修复、对势和计算的父类，并概述了它们对子类的继承。还解释了用于控制执行流程的方法。在编写自定义功能时，这些概念将是有益的，尤其是在建立执行顺序时。

在介绍了源代码框架之后，在下一章中，我们将探讨在 MD 模拟中表示物理量的不同变量、方法和数组，例如原子的位置、速度和力，这些都需要正确实现以表示 MD 模型的物理方面。

# 进一步阅读

- LAMMPS 开发人员指南 - 时间步长如何工作：https://lammps.sandia.gov/doc/Developer_flow.html

- LAMMPS 开发人员指南 – LAMMPS 类：https://lammps.sandia.gov/doc/Classes_lammps.html 

# 问题

1. `fix.cpp`中的`id`、`igroup`和`style`变量代表什么数量? 

2. `pair.h` 中哪个方法返回原子间的对势？

3. `verlet.cpp` 中的哪些方法执行速度 Verlet 算法的两半？