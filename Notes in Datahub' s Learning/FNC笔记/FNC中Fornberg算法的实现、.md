```python
def fdweights(t, m):
    """
    fdweights(t, m)
	使用向量t中的节点值来返回函数在零处的第m阶导数的权重
    """
    # This is a compact implementation, not an efficient one.

    def weight(t, m, r, k):
        # 单重递归.
        # Input:
        '''
        t: 所有的坐标点列表。
        m: 我们要求的导数阶数。
		r: 代表“当前递归窗口中最右边那个点的索引”。
	    k: 我们正在计算哪一个点的权重（目标点的索引）。
		'''
        if (m < 0) or (m > r):  
            c = 0
        elif (m == 0) and (r == 0):  
            c = 1
        else:
            if k < r:
                denom = t[r] - t[k]
                c = (t[r] * weight(t, m, r-1, k) - m * weight(t, m-1, r-1, k)) / denom
            else:
                beta = np.prod(t[r-1] - t[:r-1]) / np.prod(t[r] - t[:r])
                c = beta * (m * weight(t, m-1, r-1, r-1) - t[r-1] * weight(t, m, r-1, r-1))
        return c

    r = len(t) - 1
    w = np.zeros(t.shape)
    return np.array([ weight(t, m, r, k) for k in range(r+1) ])
```

逐行解释：
```python
if (m < 0) or (m > r):  
    c = 0
elif (m == 0) and (r == 0):  
    c = 1
```

我们在这里处理的是“意外”的情况
if框定的条件是当导数阶数小于零或者导数的阶数大于我们用的点的数量，这是因为
- 导数阶数不能小于零
- 不能通过n-m个点算n阶导数(m>0)

elif指出的是最基础的情况：如果你只有1个点（`r=0`），且要求0阶导（即函数本身），那权重就是1 这就是递归的起点

```python

else:  # generic recursion
    if k < r:
	    denom = t[r] - t[k]
        c = (t[r] * weight(t, m, r-1, k) - m * weight(t, m-1, r-1, k)) / denom
        
        else:
            beta = np.prod(t[r-1] - t[:r-1]) / np.prod(t[r] - t[:r])
            c = beta * (m * weight(t, m-1, r-1, r-1) - t[r-1] * weight(t, m, r-1, r-1))
        
        return c
```
核心逻辑在else这里

#### if(k < r):
k < r代表k这个点的权重已经算出来过了，但是现在新加入了一个点$x_r$，所以我们需要加上$x_r$点的利用去重新计算k点的权重
对应的就是我们在刚刚推导的那个公式。当k < r时，我们认为我们是在更新旧点$x_i$的权重，所以我们应用刚刚推导出来的那个公式$$\begin{aligned} w^{(k)}(z) = \frac{u^{(k)}(z)(x_i - z) - k \cdot u^{(k-1)}(z)}{x_i - x_j} \end{aligned}$$
denom就是$x_i-x_j$， t[r]就是($x_i-z$)，因为我们求的是0处的权重，所以z=0，$x_i$就为t[r]，而weight(t, m, r-1, k)就代表着上一轮$x_i$的权重，m就是导数的阶数k，weight(t, m-1, r-1, k)就代表着在$x_i$这个点上上一阶导数的权重
这样我们就完成了旧点的新权重的更新

#### else（即k >= r）:
k >= r说明这个点完全是新的，权重从来没有计算出来过，所以我们需要完全地计算一遍
但是这段代码的编写逻辑十分地巧妙：既然$x_r$的插值多项式的分子只比$x_{r-1}$的少一项$(x-x_{r-1})$，我们就可以利用$L_{r-1}(x)$的分子去迭代$L_{r}(x)$的分子$$分子_{L_r(x)}=分子_{L_{r-1}}(x)(x-x_{r-1})$$
分子找到了关系，分母也得调整。
拉格朗日基函数的分母是常数，用来归一化（保证在基准点处值为1）。
- **旧分母**：上一个点 $x_{r-1}$ 与其他所有点的距离乘积。
- **新分母**：新点 $x_r$ 与其他所有点的距离乘积。
代码里的 `beta` 就是在做这个转换：
```python
beta = np.prod(t[r-1] - t[:r-1]) / np.prod(t[r] - t[:r])
```
$$\beta = \frac{\text{旧分母}}{\text{新分母}}$$
所以，我们得到了完整的关系式：

$$L_{new}(x) = \beta \cdot L_{old}(x) \cdot (x - x_{r-1})$$
然后 我们又可以利用$$\begin{aligned} w^{(k)}(z) = \frac{u^{(k)}(z)(x_i - z) - k \cdot u^{(k-1)}(z)}{x_i - x_j} \end{aligned}$$
这个公式来求解了，c就是我们所要求的权重。