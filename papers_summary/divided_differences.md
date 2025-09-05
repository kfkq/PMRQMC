# Calculating the Divided Differences of the Exponential Function by Addition and Removal of Inputs

## Paper Metadata
- **Title**: Calculating the divided differences of the exponential function by addition and removal of inputs
- **Authors**: Lalit Gupta, Lev Barash, Itay Hen
- **Publication Date**: arXiv:1912.12157v2 [physics.comp-ph] 23 May 2020
- **Key Focus**: Efficient computation of divided differences for the exponential function (DDEF) via sequential addition/removal of inputs, enabling handling of very long input lists. Primary application: weights in off-diagonal series expansion quantum Monte Carlo (QMC).

## Abstract and Main Contribution
The paper introduces an algorithm to compute the divided differences $\exp[z_0, \dots, z_n]$ of the exponential function $f(z) = e^z$ by treating the input list $[z_0, \dots, z_n]$ as a stack. It leverages a recent identity from Zivcovich (2019) to update the DDEF when adding or removing an input, requiring only $O(sn)$ floating-point operations and $O(sn)$ memory, where $s \propto \max_{i,j} |z_i - z_j|$ (specifically, $s = \lceil \frac{1}{3.5} \max_{i,j} |z_i - z_j| \rceil$ after shifting by the mean $\bar{z}$).

This overcomes limitations of prior methods, which scale as $O(sn^2)$ for full recomputation and fail for large $n$ due to precision issues. The method handles input lists orders of magnitude longer than state-of-the-art, without needing arbitrary-precision arithmetic for most cases. A custom "ExExFloat" data type (extended exponent float) is introduced for large $n$ and $s$, combining a double-precision mantissa with a 32-bit integer exponent for range $\sim [10^{-646456992}, 10^{646456992}]$.

## Background: Divided Differences of the Exponential Function (DDEF)
Divided differences for a function $f(\cdot)$ over inputs $[z_0, \dots, z_n]$ are defined as:
\[
f[z_0, \dots, z_n] \equiv \sum_{j=0}^n \frac{f(z_j)}{\prod_{k \neq j} (z_j - z_k)}.
\]
For repeated inputs $z_0 = \dots = z_n = x$:
\[
f[x, \dots, x] = \frac{f^{(n)}(x)}{n!}.
\]
For $f(z) = e^z$, the DDEF is $\exp[z_0, \dots, z_n]$.

Traditional computation uses recurrence:
\[
f[z_i, \dots, z_{i+j}] = \frac{f[z_{i+1}, \dots, z_{i+j}] - f[z_i, \dots, z_{i+j-1}]}{z_{i+j} - z_i},
\]
with $f[z_i] = f(z_i)$. However, this is numerically unstable in double precision.

Opitz's formula relates DDEFs to matrix exponentials:
\[
F(z_0, \dots, z_n) = \begin{pmatrix}
\exp[z_0] & \exp[z_0, z_1] & \cdots & \exp[z_0, \dots, z_n] \\
& \exp[z_1] & \cdots & \exp[z_1, \dots, z_n] \\
& & \ddots & \vdots \\
& & & \exp[z_n]
\end{pmatrix},
\]
but prior methods (e.g., McCurdy et al., 2004) are inefficient for large $n$.

## Zivcovich's Algorithm (Basis for the Proposed Method)
Zivcovich's method computes the vector $\mathbf{f} = (\exp[z_0], \exp[z_0, z_1], \dots, \exp[z_0, \dots, z_n])^T$ in $O(sn^2)$ operations.

- Scale inputs by $s \geq \lceil \frac{1}{3.5} \max_i |z_i - \mu| \rceil$, shifting by $\mu = \bar{z}$ (mean of $z_i$).
- Initialize $\mathbf{h}^{(-1)} = (1, (1! \cdot s)^{-1}, (2! \cdot s^2)^{-1}, \dots, (N! \cdot s^N)^{-1})^T$, where $N \geq n + 30$.
- Update vectors $\mathbf{h}^{(j)}$ sequentially:
  \[
  \mathbf{h}^{(j)} = H_1(z_j - z_{j-1}) H_2(z_j - z_{j-2}) \cdots H_N(z_j - z_{j-N}) \mathbf{h}^{(j-1)},
  \]
  where $H_i(z)$ is identity with $z$ added at position $(i, i+1)$.
- Builds matrix $G = R F(z_0/s, \dots, z_n/s) R^{-1}$, where $R = \operatorname{diag}(1, s, \dots, s^n)$.
- $\mathbf{f}^T = \mathbf{g}^{(s)}$, obtained by powering the first row of $G$ ($s-1$ times matrix multiplication).

Accuracy ensured by $N \geq n + 30$ and $s$ choice, based on Taylor truncation at order 30 for $|z| \leq 3.5$.

## Proposed Method: Addition and Removal of Inputs
Treat inputs as a stack; update DDEF on add/remove.

### Addition (Pushing $z_n$ to $[z_0, \dots, z_{n-1}]$)
- Assume $\mathbf{h}^{(n-1)}$ and rows $\mathbf{g}^{(i)}$ ($i=1,\dots,s$) stored.
- Compute $\mathbf{h}^{(n)} = H_1(z_n - z_{n-1}) \cdots H_N(z_n - z_{n-N}) \mathbf{h}^{(n-1)}$ in $O(n)$ ops.
- Append to $\mathbf{g}^{(1)}$: new element from $\mathbf{h}^{(n)}$.
- For $i=2,\dots,s$: append $\mathbf{g}^{(i-1)} \cdot \mathbf{w}^{(n)}$, where $\mathbf{w}^{(n)}$ is derived from $\mathbf{h}^{(n)}$.
- Total: $O(sn)$ ops, $O(sn)$ memory.
- If new $z_n$ increases $s$, restart with updated $\mu, s$. Double $N$ if $n > N-30$.

### Removal (Popping $z_n$)
- Inverse transform: $\mathbf{h}^{(n-1)} = H_N(z_{n-N} - z_n) \cdots H_1(z_{n-1} - z_n) \mathbf{h}^{(n)}$ in $O(n)$ ops.
- Remove last element from each $\mathbf{g}^{(i)}$.
- Total: $O(n)$ ops (marginal $s$ dependence).

## Modifications for Higher Precision
DDEFs decay as $\sim 1/n!$, causing underflow for large $n$.

- Compute modified $\tilde{\mathbf{f}} = (\exp[z_0], 1! \cdot \exp[z_0, z_1], \dots, n! \cdot \exp[z_0, \dots, z_n])^T$.
- Rescale vectors: $\tilde{h}^{(j)}_i = h^{(j)}_i \cdot \frac{j!}{(j-i)!}$ for $i \leq j$, else $h^{(j)}_i \cdot i!$.
- Modified matrices $\tilde{G}_{ij} = G_{ij} \cdot j!/i!$, rows $\tilde{\mathbf{g}}^{(i)}_j = \mathbf{g}^{(i)}_j \cdot j!$.
- Modified operators $\tilde{H}^{(j)}_i(z)$ adjust for factorials.
- For $s=1$, accurate up to $n \approx 10^5$ in double precision (vs. $n \approx 100$ without modification).

For large $n,s$: Use ExExFloat (double mantissa + 32-bit int exponent) for extended range without mantissa precision loss. Arithmetic ops cost similar to doubles; overall slowdown $\approx 2.7\times$ vs. double precision.

## Numerical Testing
- **Precision Comparison**: For $s=1$, standard double fails at $n \approx 100$; improved lasts to $n \approx 10^5$. For $s=2$, double fails at $n \approx 1700$; ExExFloat unlimited.
- **Runtimes**: Addition scales $O(sn)$; removal $O(n)$. Benchmarked on Intel i5-8257U; linear in $n$ and $s$.
- Error analysis: Relative error after $N_f$ ops $\approx 10^{-17} \sqrt{N_f} \ll 10^{-6}$ for $N_f=10^{22}$.

## Application to Off-Diagonal Expansion Quantum Monte Carlo (ODE-QMC)
In ODE-QMC, partition function $Z = \operatorname{Tr}(e^{-\beta H})$ expands as sum over configurations with weights involving DDEFs: $w = \exp[z_0, \dots, z_n]$ where $z_i$ relate to off-diagonal matrix elements.

- Configurations updated by adding/removing terms, requiring fast DDEF updates.
- Method enables efficient weight recalculation, avoiding $O(sn^2)$ recomputes.
- Demonstrates applicability to quantum many-body simulations (e.g., thermal properties).

## Significance and Outlook
Enables DDEF for unlimited $n$ without arbitrary precision slowdown. Open-source C++ code on GitHub. Future: Broader QMC applications, complex inputs.