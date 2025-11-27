
# Levenberg方法

- 拟牛顿法 Broyden公式 前向差分 与Levenberg的关系
- 理论公式 怎么来的？
- 实现

## 1、拟牛顿法 Broyden公式 前向差分 与Levenberg的关系

- ### broyden公式满足拟牛顿条件
	 在FNC书本中我们可以看到broyden的公式如下
     $$ \mathbf{A}_{k+1} = \mathbf{A}_k + \frac{1}{\mathbf{s}_k^T \mathbf{s}_k} (\mathbf{y}_{k+1} - \mathbf{y}_k - \mathbf{A}_k \mathbf{s}_k) \mathbf{s}_k^T $$
	 那它是怎么来的呢？
	 