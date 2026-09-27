# AccelNet CPU/GPU speedup methods

**Document version 1.13 — 2026-09-27 (JST).**

This document records the mathematics, implementation decisions, and measurements
behind the CPU/GPU optimizations in this working tree. All numerical kernels use
FP64. These changes reorganize exact formulas, remove redundant work, and improve
data placement and execution; they do not approximate the potential. Floating-point
operation order changes, so agreement means agreement within tested tolerances,
not bitwise equality.

## 1. Versions and scope

### 1.1 Software and hardware

| Item | Version or identity |
|---|---|
| AccelNet / AccelNetPredictor | **1.0.1**, as declared by CMake |
| Base Git commit | `c6631460a1bbb990c82e3e0ff5e73c36c52f6f9b` |
| Base `git describe --tags --always` | `1.0.0-6-gc663146` |
| Optimization source revision | `gpu` checkpoint **`1d6985d`** (revision 1.6), followed by the common G5 moments in Section 18 and n2p2 extensions in Section 19 and grouped/LAMMPS evaluation in Section 20, exact high-order moments/threading in Section 21, atomic removal in Section 22, and unified inference APIs in Section 23, and energy-only recovery in Section 24; each validation archive identifies its measured sources |
| LAMMPS | **29 Aug 2024 Update 4**, with this repository's ACCELNET/GPU adapter and triclinic patch |
| GNU Fortran | **11.4.0**, Ubuntu `11.4.0-1ubuntu1~22.04`; CPU `-O3` |
| NVIDIA HPC SDK / nvfortran | **25.3 / 25.3-0**; CPU `-fast -O3`; GPU `-mp=gpu -gpu=cc90,cc120` |
| CUDA toolkit | **12.8**, `nvcc V12.8.93` |
| LAMMPS GPU configuration | `GPU_API=cuda`, `GPU_PREC=double` |
| CPU | Intel Xeon Gold **6526Y**; single-core timing pinned to CPU 6 |
| GPU A | NVIDIA **H100 NVL**, UUID `GPU-2644154d-7268-af42-6631-59e1f3c6e7f3` |
| GPU B | NVIDIA **RTX PRO 6000 Blackwell Max-Q**, UUID `GPU-c8a7b723-d4eb-a67d-5ebd-139cedd981c9` |

The document version, library version, and Chebyshev descriptor `version` are
separate identifiers. Exact compiler versions were checked when writing this
record. Each validation report and its JSON files specify the commands and inputs
used in that measurement campaign.

The validation archives retain the original hashes and descriptions of the
pre-commit working trees used for their measurements. The initial `gpu` checkpoint `cd00106`
collects the implementation, tests, and records; commit preparation only removed
trailing whitespace from one source line and added Python-cache ignore rules.

The **moment-Horner baseline** LAMMPS executable, deployed before revision 1.1,
has SHA256:

```text
fa20e636e98e4f7229233f3e25f3eba24d7d24808a38069214113421bcddd61e
```

The [Horner archive](docs/validation/moment-horner-2026-09-26/README.md) includes
[binary hashes](docs/validation/moment-horner-2026-09-26/binaries.json), its
[implementation diff](docs/validation/moment-horner-2026-09-26/implementation.diff),
and baseline/prototype sources. That diff covers the Horner change, not all earlier
uncommitted changes. It does not identify later generic-descriptor experiments.

The revision-1.2 validated/deployed LAMMPS executable has SHA256
`9a7200fba689020f9f53381de9ed6dd22e973d3e9c24bb440ed9af4700f8b5b5`. Its
[source/build identities and validation](docs/validation/behler-contraction-2026-09-26/README.md)
include the final reuse of NN-gradient slots, not the intermediate coefficient-array trial.

The revision-1.3 validated/deployed executable has SHA256
`6aea0e455b53f6033a55ab4f450db27f2231c9cc6453c30a37b6f524eb27ce8c`. Its source, timings, and validation
are recorded in the [G4 report](docs/validation/g4-common-2026-09-27/README.md).

### 1.2 Adoption status

| Method | Status and scope |
|---|---|
| CSR batching and persistent GPU model/work buffers | Implemented; descriptors, NN, forces, and virial execute on the GPU |
| Chebyshev direct and moment | Implemented in a common numerical source for supported CPU CSR batches and GPU execution |
| Reused unit directions and differentiated Clenshaw | Adopted in common code |
| Precomputed moment force coefficients and contracted polynomial gradient | Adopted in common code |
| Differentiated multivariate Horner and lexicographic coefficient packing | Adopted and validated |
| LJ and Behler G1–G5 | GPU support implemented; G4 uses direct pairs; G5 supports shared direct/moment evaluation (Section 18) |
| Generic radial caching and LJ component fusion | Implemented in the common serial CPU/GPU source; measurements in Section 12 |
| Behler angular coefficient contraction and differentiated Horner | Retained for G5; the former G4 path in Sections 13–14 is superseded by Section 15 |
| G4 value/Jacobian evaluation | Shared value/derivative loop with scalar caches and disjoint descriptor owners; CPU uses one owner, GPU uses a flat center/owner launch; force contraction reads the saved Jacobian; Sections 15–16 |
| LAMMPS GPU package integration | CUDA adapter plus Fortran OpenMP target; AMD/HIP/OpenCL interoperability is not implemented |

The default CPU shared path applies to the **CSR batch API** for supported
Chebyshev, LJ and Behler G1–G5 models (revision 1.6, Section 17).
`evaluate_batch_reference` retains the independent old CPU path and is the
fallback for mixed/multiple Chebyshev components within one element. G5 moments
now use the common source. Auto retains the original per-component threshold of
16 angular neighbors and the order-10 bound (revision 1.7, Section 18). Object APIs, atomic Fortran/C APIs, CLI, and
ordinary LAMMPS `pair_style accelnet` still use their established CPU implementation.

## 2. Notation and force contraction

For center atom $i$, let $j\in\mathcal N_i$ enumerate neighbors, including periodic
images. We suppress $i$ when one center is fixed:

$$
\mathbf d_j=\mathbf R_j+\mathbf H\mathbf n_j-\mathbf R_i,\qquad
r_j=\lVert\mathbf d_j\rVert,\qquad
\mathbf u_j=\frac{\mathbf d_j}{r_j},\qquad
c_{jk}=\mathbf u_j\cdot\mathbf u_k.
$$

Let $N$ be the number of centers, $z$ the average angular-neighbor count, $p$ the
maximum angular degree, and $B$ the number of NN inputs. Write $f_j=f_c(r_j)$,
$f'_j=df_c(r_j)/dr_j$, and $s_j$ for the species weight. Chebyshev has channels
$w_j^{(0)}=1$ and $w_j^{(1)}=s_j$; omit channel 1 when it is not used.
Repeated atom indices can denote distinct periodic images. Do not deduplicate CSR
edges by atom index alone.

### 2.1 Chain rule, including normalization

With normalized inputs $x_b=(G_b-\mu_b)\sigma_b$ and NN output $y(\mathbf x)$,

$$
E_i=\frac{y(\mathbf x_i)}{s_E}+e_{\rm shift}+e_{{\rm ref},i},\qquad
 g_{ib}=-\frac{\partial E_i}{\partial G_{ib}}
 =-\frac{\sigma_b}{s_E}\frac{\partial y}{\partial x_{ib}}.
$$

Here $g$ is the **negative** energy gradient. The neighbor force contribution is

$$
\mathbf F_{i\to j}=\sum_b g_{ib}\frac{\partial G_{ib}}{\partial\mathbf d_j}.
$$

Materializing the full descriptor Jacobian requires approximately $3NzB$ scalars,
plus their writes and reads. Instead, compute $g$ by NN backpropagation and contract
descriptor derivatives directly into three force components per edge. This removes
the $O(NzB)$ Jacobian storage, though moment and other work arrays still remain.

### 2.2 Force and virial accumulation

Compute edge contributions independently, then scatter by CSR row:

$$
\mathbf F_j\mathrel{+}=\mathbf F_{i\to j},\qquad
\mathbf F_i\mathrel{-}=\sum_j\mathbf F_{i\to j},\qquad
W_{ab}\mathrel{+}=\sum_j d_{j,a}F_{i\to j,b}.
$$

Shared GPU outputs use FP64 atomic updates. Form virial contributions using each
image displacement before combining repeated atom indices: periodic self images
can produce a virial even when their atom forces cancel. No extra factor $1/2$
belongs here; these are derivatives of the individual contributions in $E=\sum_iE_i$.
The public batch API adds forces and virial to caller accumulators and overwrites
row energies.

## 3. Direct Chebyshev evaluation and descriptor versions

### 3.1 Descriptor definitions

For radial degree $p_r$ and radial cutoff $R_c$,

$$
G_{n,\mathrm{rad}}^{(h)}
=\sum_jw_j^{(h)}f_c^{\rm rad}(r_j)
T_n\!\left(\frac{2r_j}{R_c}-1\right),\qquad 0\le n\le p_r.
$$

The angular descriptors sum **unordered** neighbor pairs:

$$
G_{n,\mathrm{ang}}^{(h)}
=\sum_{j<k}w_j^{(h)}w_k^{(h)}f_jf_kT_n(a c_{jk}+b),\qquad 0\le n\le p.
$$

| Descriptor `version` | Angular argument $a c+b$ | Radial behavior |
|---|---|---|
| `0` | $c$; $a=1,b=0$ | Neighbor sum above |
| `1` | $2c/\pi-1$; $a=2/\pi,b=-1$ | Neighbor sum above |
| `10` | $c$ | Also adds the central contribution $T_n(-1)=(-1)^n$ |

Version 1 uses the legacy **affine function of cosine**, not $\arccos(c)$.
For multielement version 10, the historical center-weight lookup through the
neighbor array and its minimum-neighbor-count validation are preserved. The
central term has zero coordinate derivative but affects NN inputs and hence $g$.
Fortran target supports 0/1/10. The current LAMMPS GPU interface uses version 0,
matching the existing LAMMPS CPU interface.

### 3.2 Reuse geometry

$$
\frac{\partial c_{jk}}{\partial\mathbf d_j}
=\frac{\mathbf u_k-c_{jk}\mathbf u_j}{r_j},\qquad
\frac{\partial\mathbf u_j}{\partial\mathbf d_j}
=\frac{I-\mathbf u_j\mathbf u_j^{\mathsf T}}{r_j}.
$$

Precompute $r_j,\mathbf u_j,f_j,f'_j,s_j$. Reuse unit directions instead of
renormalizing displacements inside pair loops. Radial forces also multiply saved
$\mathbf u_j$ rather than computing $\mathbf d_j/r_j$ again. Small savings matter
when repeated for $O(Nz^2)$ pairs.

### 3.3 Differentiated Clenshaw for direct angular forces

For one pair, define

$$
A_n=g_{n,\mathrm{ang}}^{(0)}+s_js_kg_{n,\mathrm{ang}}^{(1)},\qquad
S(x)=\sum_{n=0}^p A_nT_n(x),\qquad x=a c_{jk}+b.
$$

Instead of forming every $T_n,T'_n$ and then contracting, use the backward recurrence

$$
B_{p+1}=B_{p+2}=D_{p+1}=D_{p+2}=0,
$$
$$
B_k=A_k+2xB_{k+1}-B_{k+2},\qquad
D_k=2B_{k+1}+2xD_{k+1}-D_{k+2}\quad(k=p,\ldots,1),
$$
$$
S=A_0+xB_1-B_2,\qquad
S_x=B_1+xD_1-D_2,\qquad S_c=aS_x.
$$

Here $D_k=dB_k/dx$. The force contribution is

$$
\mathbf F_{i\to j}^{(jk)}
=f'_jf_kS\,\mathbf u_j
+\frac{f_jf_k}{r_j}S_c(\mathbf u_k-c_{jk}\mathbf u_j).
$$

This evaluates value and derivative with a constant number of temporaries and no
additional array or kernel launch. Degree zero gives $S=A_0,S_c=0$. Complexity
remains $O(p)$, with a smaller arithmetic constant. The adopted change applies to
**direct angular force contraction**, not all radial and descriptor-value loops.
Values still visit $j<k$ once, while edge forces visit both directed pairs.

## 4. Moment method: replace pair sums with moment products

### 4.1 Polynomial coefficients and moments

Expand the Chebyshev polynomial in powers of cosine:

$$
T_n(a c+b)=\sum_{q=0}^n P_{nq}c^q.
$$

Compute coefficients once during model initialization. With out-of-range entries zero,

$$
P_{00}=1,\quad P_{10}=b,\quad P_{11}=a,\qquad
P_{nq}=2bP_{n-1,q}+2aP_{n-1,q-1}-P_{n-2,q}.
$$

For multi-index $\boldsymbol\alpha=(\alpha_x,\alpha_y,\alpha_z)$, write
$|\boldsymbol\alpha|=q$ and
$\mathbf u^{\boldsymbol\alpha}=u_x^{\alpha_x}u_y^{\alpha_y}u_z^{\alpha_z}$. Define

$$
m_{\boldsymbol\alpha}=\frac{q!}{\alpha_x!\alpha_y!\alpha_z!},\qquad
M_{\boldsymbol\alpha}^{(h)}=\sum_jw_j^{(h)}f_j\mathbf u_j^{\boldsymbol\alpha}.
$$

The number of components through degree $p$ is

$$
K(p)=\sum_{q=0}^p\binom{q+2}{2}=\binom{p+3}{3}
=\frac{(p+1)(p+2)(p+3)}6.
$$

For example, $K(5)=56$ and $K(10)=286$: doubling the degree more than doubles the work.

### 4.2 Pair identity and diagonal correction

The multinomial theorem gives

$$
c_{jk}^q=\sum_{|\boldsymbol\alpha|=q}
 m_{\boldsymbol\alpha}\mathbf u_j^{\boldsymbol\alpha}\mathbf u_k^{\boldsymbol\alpha}.
$$

Therefore,

$$
D_q^{(h)}=\sum_{j<k}w_j^{(h)}w_k^{(h)}f_jf_kc_{jk}^q
=\frac12\left[
\sum_{|\boldsymbol\alpha|=q}m_{\boldsymbol\alpha}(M_{\boldsymbol\alpha}^{(h)})^2
-\sum_j(w_j^{(h)}f_j)^2\right],
$$
$$
G_{n,\mathrm{ang}}^{(h)}=\sum_{q=0}^nP_{nq}D_q^{(h)}.
$$

The subtracted term removes $j=k$. It is independent of degree because
$\mathbf u_j\cdot\mathbf u_j=1$, and is needed even at degree zero or with only one
neighbor. This is an exact identity, not an approximation to the angular distribution.

### 4.3 Complexity and automatic selection

Ignoring radial and NN work,

$$
T_{\rm direct}=O(Nz^2(p+1)),\qquad
T_{\rm moment}=O(NzK(p)+NK(p)+N(p+1)^2).
$$

Both moment construction and force evaluation contain $NzK$ work, but scale linearly
with neighbor count. High density favors moments; high degree, few neighbors, or
large work arrays can favor direct evaluation. The shared target code chooses
moments for row $i$ when

$$
z_i(p+1)\ge K(p).
$$

This operation-count heuristic uses the actual angular-cutoff neighbor count.
It is not device-specific autotuning and does not model launch costs, bandwidth,
parallelism, or NN fraction.

## 5. Moment forces as one polynomial and its gradient

### 5.1 Prepare coefficients once per center

Transform the negative NN gradient into the power basis:

$$
\delta_q^{(h)}=\sum_{n=q}^p g_n^{(h)}P_{nq},\qquad
C_{\boldsymbol\alpha}^{(h)}
=m_{\boldsymbol\alpha}\delta_{|\boldsymbol\alpha|}^{(h)}M_{\boldsymbol\alpha}^{(h)}.
$$

These quantities are independent of neighbor $j$: compute them **once per center**,
not once per edge and monomial. Raw moments are dead after descriptor/NN evaluation,
so overwrite their storage with force coefficients and rebuild moments on the next
call. No additional large coefficient array is required.

For a neighbor with species weight $s$, define

$$
\mathcal P_s(\mathbf u)=\sum_{|\boldsymbol\alpha|\le p}
(C_{\boldsymbol\alpha}^{(0)}+sC_{\boldsymbol\alpha}^{(1)})\mathbf u^{\boldsymbol\alpha},
\qquad
Q_s=\sum_{q=0}^p(\delta_q^{(0)}+s^2\delta_q^{(1)}).
$$

Pre-sum the two channels needed for $Q_s$ once per center as well.

### 5.2 Contract the gradient before projecting

The angular force becomes

$$
\boxed{
\mathbf F_{i\to j}^{\rm ang}
=f'_j[\mathcal P_{s_j}(\mathbf u_j)-f_jQ_{s_j}]\mathbf u_j
+\frac{f_j}{r_j}\left[
\nabla_{\mathbf u}\mathcal P_{s_j}(\mathbf u_j)
-\mathbf u_j(\mathbf u_j\cdot\nabla_{\mathbf u}\mathcal P_{s_j}(\mathbf u_j))
\right].}
$$

Differentiate $\sum_{h,q}\delta_q^{(h)}D_q^{(h)}$ to obtain this formula. The square's
factor 2 cancels $1/2$, the diagonal term gives
$-f_jf'_jQ_{s_j}\mathbf u_j$, and differentiating the direction supplies
$(I-\mathbf u\mathbf u^{\mathsf T})/r$.

Sum the Cartesian polynomial gradient first and project onto the unit sphere's
tangent plane **once per edge**, rather than once per monomial. Linearity preserves
the result while reducing geometric arithmetic and inner products.

### 5.3 Reuse powers and order construction work by center

Precompute $u_\gamma^0=1$ and $u_\gamma^{k+1}=u_\gamma^ku_\gamma$. Each monomial then
uses three saved powers. GPU construction collapses a loop with `row` outside and
`entry` inside, placing components of one center near one another in execution
order. This aims to reuse that center's neighbor geometry/powers before eviction.

Measurements support the change, but hardware cache-miss counters were not collected.
This is not a transpose of every array or a claim that every memory access becomes
contiguous.

## 6. Differentiated multivariate Horner for moment forces

### 6.1 Why Horner rather than direct reuse of Clenshaw

Direct forces retain $\sum_nA_nT_n(c)$, so Clenshaw applies naturally. Moment forces
already require a power-basis polynomial

$$
\mathcal P(x,y,z)=\sum_{a+b+c\le p}C_{abc}x^ay^bz^c
$$

and its three derivatives. Here $x,y,z$ are unit-direction components and
$C_{abc}=C_{abc}^{(0)}+s_jC_{abc}^{(1)}$. Use multivariate Horner on this expression.

### 6.2 Simultaneous value/derivative recurrence

Differentiating the one-variable update gives

$$
d_{\rm new}=td_{\rm old}+v_{\rm old},\qquad
v_{\rm new}=tv_{\rm old}+c.
$$

The derivative must use the old value. Nest these recurrences in $z$, then $y$, then
$x$. Execute assignments in the following order:

```text
v = dx = dy = dz = 0
for a = p, ..., 0:
    hy = dhy = hyz = 0
    for b = p-a, ..., 0:
        hz = dhz = 0
        for c = p-a-b, ..., 0:
            dhz = z*dhz + hz
            hz  = z*hz  + C[a,b,c]
        dhy = y*dhy + hy
        hyz = y*hyz + dhz
        hy  = y*hy  + hz
    dx = x*dx + v
    dy = x*dy + dhy
    dz = x*dz + hyz
    v  = x*v  + hy
# v = P; (dx,dy,dz) = grad(P)
```

Complexity remains $O(K(p))$, but fewer products and power-array reads are needed.
Unlike $au^a/u$, the recurrence has no division by a direction component, so zeros
need no special treatment. Degree zero gives a zero gradient. Saved powers remain
necessary for moment construction, but are no longer read by moment forces.

### 6.3 Pack coefficients in evaluation order

The first prototype preserved degree-major storage and looked up each coefficient's
index. The adopted version packs moments lexicographically in ascending $a,b,c$.
With $R=p+1-a$, the one-based index is

$$
I(a,b,c)=K(p)-\frac{R(R+1)(R+2)}6+bR-\frac{b(b-1)}2+c+1.
$$

Reverse traversal then reads coefficient numbers $K,K-1,\ldots,1$, without an index
lookup per term. The relative $a,b$ order within each total degree remains unchanged,
preserving descriptor degree-reduction order. Element-specific degrees use their
own $p$ and $K(p)$.

The Fortran layout remains `moments(row,entry,channel)`. Coefficient **numbers** are
consecutive, but one row reads with a stride equal to the row capacity; this does
not imply physically contiguous loads. No additional GPU kernel or large scratch
array is introduced.

### 6.4 What the slower prototype taught us

For H100, 512 atoms, degree 10, force-stage time relative to the baseline was
**1.120** with index lookup and **0.822** with lexicographic packing. The prototype
added dependent index-to-coefficient loads and widened metadata from four to five
integers. The adopted version removes the lookup and restores four integers.
Coefficient layout mattered even with the same Horner recurrence.

Static force-kernel register counts changed from 86 to 80 on H100 and 96 to 92 on
Blackwell. Both Horner variants had the same register counts; stack sizes were
unchanged. Increased register use therefore does not explain the prototype's
slowdown. The separate contributions of lookup and metadata width were not isolated;
cache misses, stalls, and achieved occupancy were not measured.

## 7. Other descriptors: evidence for removing repeated work

This section records the original CPU diagnostic experiments. Subsequent production
CPU/GPU changes must be distinguished from these diagnostic measurements.

### 7.1 Cache shared radial factors for G4/G5

For matching radial parameters,

$$
q(r)=e^{-\eta(r-R_s)^2}f_c(r),\qquad
q'(r)=e^{-\eta(r-R_s)^2}[f'_c(r)-2\eta(r-R_s)f_c(r)]
$$

can be computed once per edge. Define

$$
A(c)=\left(\frac{1+\lambda c}{2}\right)^\zeta,\qquad
A'(c)=\frac{\zeta\lambda}{2}\left(\frac{1+\lambda c}{2}\right)^{\zeta-1}.
$$

This implementation's unordered-pair values are

$$
G^{\rm G5}_{jk}=2A(c_{jk})q(r_j)q(r_k),\qquad
G^{\rm G4}_{jk}=2A(c_{jk})q(r_j)q(r_k)q(r_{jk}).
$$

For example,

$$
\frac{\partial G^{\rm G5}_{jk}}{\partial\mathbf d_j}
=2q(r_k)\left[Aq'(r_j)\mathbf u_j
+\frac{A'q(r_j)}{r_j}(\mathbf u_k-c_{jk}\mathbf u_j)\right].
$$

G4 also needs the neighbor–neighbor factor and its derivative; an edge cache alone
does not remove this pair-dependent work. Grouping identical
$(\eta,R_s,R_c,\mathrm{cutoff},\alpha_{\rm cutoff})$ avoids repeated exponentials,
cutoffs, and derivatives and limits storage to parameter-group count rather than
feature count. Existing shared scalar formulas already specialize small integer
angular powers.

The original generic implementation visits each unordered pair once for values and
twice for directed edge forces, and its value stage calls a value/derivative helper.
The old CPU G4 implementation obtains values and both Jacobians in one visit; CPU
G5 direct uses one visit for values and one for forces. Independent edge parallelism
avoids a full Jacobian and force conflicts but performs extra arithmetic.

OpenMP-disabled diagnostic time divided by old CPU time, 512 atoms:

| Family | GNU original common | GNU radial cache | NVHPC original common | NVHPC radial cache |
|---|---:|---:|---:|---:|
| G4 | 4.801 | 2.558 | 3.590 | 2.336 |
| G5 | 3.281 | 1.509 | 2.092 | 1.285 |
| Mixed Behler | 2.561 | 1.299 | 1.726 | 1.112 |

These diagnostics include cache construction/allocation and demonstrate redundant
arithmetic as a cause of slowdown. They do not establish GPU performance or parity
with the optimized legacy CPU loops.

### 7.2 Fuse LJ components

These are NN descriptors based on $r^{-6}$ and $r^{-12}$, not a standalone classical
LJ pair potential:

$$
t=r^{-2},\quad h_6=t^3,\quad h_{12}=h_6^2,\qquad
G_m=\sum_jf_c(r_j)r_j^{-m}\quad(m=6,12),
$$
$$
\frac{d}{dr}[f_c(r)r^{-m}]=r^{-m}\left[f'_c(r)-\frac mr f_c(r)\right].
$$

Evaluate the two components together to share $f_c,f'_c,1/r,h_6$.
For 512 atoms, diagnostic ratios to old CPU were GNU
1.618 → 1.639 → 1.056 and NVHPC 1.499 → 1.424 → 0.999 for original, powers-only,
and fused implementations. Rewriting powers alone was insufficient; sharing work
between components mattered. See the [diagnostic report and sources](docs/validation/gpu-descriptors-cpu-comparison-2026-09-26/README.md).

## 8. GPU residency, data layout, and LAMMPS

### 8.1 Avoid repeated allocation and model transfers

Pack model weights, metadata, normalization, and species maps into plain numerical
arrays rather than deeply mapping nested derived types. `target_workspace` owns
persistent device mappings and grows buffers only when capacity is insufficient.
Upload the model only when it changes. Compare packed values, not merely shapes,
to detect same-size weight reloads and species-map changes. Tests verify stable
allocation/upload counters in steady state.

Copies own model values independently. Assignment, release, and reinitialization
must preserve workspace ownership and avoid stale cache reuse. The ordinary CPU
shared batch repacks the publicly mutable model each call while reusing scratch;
benchmark this separately from an already prepared model.

NN arrays keep rows in the first dimension. Parallel work is organized separately
for per-center NN, center/component moments, per-edge forces, and per-center scatter.
The complete computation is not one large serial neighbor loop per GPU thread.

### 8.2 LAMMPS GPU package responsibilities

`pair_style accelnet/gpu` uses the GPU package for device assignment and neighbor
arrays, while Fortran OpenMP target computes descriptors, NN, and forces.

- `neigh no`: pack the CPU full neighbor list as CSR and upload it.
- `neigh yes/hybrid`: borrow lib/gpu device coordinates/neighbors and convert them
  to CSR/displacements on device. No neighbor-list download/reupload is needed;
  one integer edge count is returned to the CPU.
- In this LAMMPS release's standard CUDA configuration, binning is still on the CPU
  even with `neigh yes`; candidate-neighbor searching runs on the GPU.
- Synchronize at the CUDA/OpenMP boundary and download local/ghost force contributions
  within `compute()`. LAMMPS then performs ghost reverse communication.

The tested interface requires FP64, `newton on`, and `split 1` for this GPU package.
LAMMPS GPU uses n2p2 2G models after conversion to embedded AccelNet models.
See [build instructions and restrictions](docs/lammps-gpu.md). Portable Fortran
arithmetic does not imply that CUDA interoperability has been ported to AMD.

### 8.3 Why speedup changes with atom count

A simplified timing model is

$$
T_{\rm GPU}=T_{\rm check/pack}+T_{\rm H2D}+L\tau_{\rm launch}
+T_{\rm descriptor}+T_{\rm NN}+T_{\rm force/scatter}+T_{\rm D2H/sync}.
$$

With arithmetic work $F$, memory traffic $D$, and effective throughput/bandwidth,

$$
T_{\rm kernel}\gtrsim
\max\left(\frac F{P_{\rm effective}},\frac D{B_{\rm effective}}\right).
$$

Small systems expose fixed costs and insufficient parallelism. Larger systems
amortize those costs, but also enlarge moment buffers and change locality and
achievable throughput. **Speedup need not increase monotonically with atom count.**
Even at fixed $N$, neighbor count and degree change the direct/moment balance.
The initial moment code also had avoidable coefficient recomputation and work-order
problems; these were improved rather than explained away as inherent GPU behavior.

For a fraction $f$ of total time accelerated by $s$,

$$
S_{\rm total}=\frac1{(1-f)+f/s}.
$$

For example, the early Blackwell residency/moment change reduced a 512-atom,
degree-8 GPU batch from 12.631 to 1.393 ms. Yet the route including CPU neighbor
construction took 123.161 ms versus 137.666 ms for CPU, only 1.12× overall.
Kernel time, synchronous API time, and full MD time are different measurements.
Source: [residency report](docs/validation/gpu-residency-moments-2026-09-25/README.md).

## 9. Fair CPU comparison and correctness checks

### 9.1 OpenMP off is not OpenMP with one thread

`OMP_NUM_THREADS=1` can retain runtime entry, synchronization, atomics, and code
selection overhead. Compare three paths:

1. Independent established CPU code.
2. Common numerical code compiled **without OpenMP compiler flags**.
3. Common code with OpenMP enabled, explicitly executing on one host thread.

For the default shared Chebyshev CPU batch, CMake generates renamed modules with
OpenMP directives removed. CPU and GPU instances can coexist in one program.
The diagnostic serial target uses `ACCELNET_TARGET_SERIAL=ON`. GNU has no libgomp
dependency; NVHPC may link libnvomp through its Fortran runtime but the serial
kernels contain no OpenMP launch calls.

The no-neighbor one-atom experiment showed roughly 2–4.5 microseconds of additional
OpenMP time, including NN and other fixed costs; it was not a pure fork measurement.
Millisecond-scale LJ/G4/G5 gaps remained with OpenMP removed and required arithmetic
analysis.

In NVHPC 25.3, `device(omp_get_initial_device())` alone unexpectedly executed a
probe on the GPU. Those initial host timings were discarded. The final path uses
`target if(...)`, checks actual execution location at initialization, and was
verified to issue no CUDA launches on the host route. GPU tests reject accidental
CPU fallback.

### 9.2 Paired timing and measurement boundaries

Pin the CPU, warm up, use identical structures/models/modes/outputs, and alternate
reference and candidate timings. Report

$$
R=\operatorname{median}_k\left(
\frac{t_{{\rm candidate},k}}{t_{{\rm reference},k}}\right).
$$

For successive common-code versions, also compare $R_{\rm after}/R_{\rm before}$.
Paired ratios reduce sensitivity to shared-machine clock changes; a ratio of
medians obtained in separate runs can be misleading. Reverse binary order between
rounds. Do not claim significance for very small changes.

Steady-state batch timings exclude neighbor construction, first allocation, and
initial model packing, but include ordinary validation, copies, transfers,
synchronization, and all E/F/virial outputs. Ordinary CPU API timings additionally
include per-call packing. LAMMPS timings use warmed-up Loop time per step and exclude
dump writing from the measured interval.

### 9.3 Numerical validation

- Compare all energies, forces, and virial components with an independent CPU
  reference. Batch comparisons use absolute and relative tolerances of `2e-10`.
- Check coordinate finite differences for forces and strain finite differences for
  virial; agreement between two instances of the same new source is insufficient.
- Include degree zero, higher degrees, zero direction components, cosine ±1,
  isolated atoms, periodic images, mixed element degrees, and different NN shapes.
- Exercise cutoffs, descriptor versions, mixed direct/moment rows, model reloads,
  buffer growth, and multiple handles.
- Use GNU bounds checking, GPU Compute Sanitizer, and MPI/ghost, triclinic, and short
  NVE trajectory comparisons where applicable.

$$
F_a\simeq-\frac{E(\mathbf R+\epsilon\mathbf e_a)
-E(\mathbf R-\epsilon\mathbf e_a)}{2\epsilon}.
$$

For deformation $\mathbf d'=(I+\varepsilon)\mathbf d$, this virial convention gives
$\partial E/\partial\varepsilon_{ba}=-W_{ab}$. Finite-difference steps and tolerances
are test-specific. Sanitizer zero errors applies to executed cases; it does not
establish that runtime memory pools are leak-free.

## 10. Measured effects

These campaigns have different baselines and conditions. **Do not multiply their
speedups to construct a cumulative speedup.**

### 10.1 Major changes

| Change | Conditions | Before → after | Evidence |
|---|---|---|---|
| Residency, moments, parallelism | Blackwell, synthetic 512 atoms, degree 8, GPU batch | 12.631 → 1.393 ms, 9.07× | [Report](docs/validation/gpu-residency-moments-2026-09-25/README.md) |
| Moment work order, coefficients, gradient contraction | H100, Ti/O 24,000 atoms, LAMMPS moment | 25.069 → 10.813 ms/step, 2.32× | [Report](docs/validation/gpu-moment-force-2026-09-26/README.md) |
| Same | Blackwell, same real model | 54.249 → 17.189 ms/step, 3.16× | Same report |
| Unit-direction reuse and direct Clenshaw | H100, Ti/O 24,000 atoms, LAMMPS direct | 27.700 → 24.331 ms/step, 12.2% reduction | [Report](docs/validation/chebyshev-common-2026-09-26/README.md) |
| Moment Horner and coefficient packing | H100, same real model, LAMMPS moment | 9.754 → 8.925 ms/step, 8.5% reduction | [Report](docs/validation/moment-horner-2026-09-26/README.md) |

During the Clenshaw campaign, moment also changed from 10.546 to 9.745 ms/step.
Moment does not execute the Clenshaw branch: that change also reused directions
and recompiled a shared force kernel. It is not evidence of Clenshaw alone speeding
up moments.

### 10.2 Horner: synchronous GPU batch

Two species, input–8–4–1 NN, spacing 1.7, forced moments. The baseline already
includes the direct/Clenshaw optimization. Source: [Horner report](docs/validation/moment-horner-2026-09-26/README.md).

| GPU | Atoms | Degree | Before [ms] | After [ms] | Total reduction | Force/scatter reduction |
|---|---:|---:|---:|---:|---:|---:|
| H100 | 512 | 5 | 0.6293 | 0.6278 | 0.2% | 9.6% |
| H100 | 512 | 10 | 1.0713 | 1.0316 | 3.7% | 17.8% |
| H100 | 4096 | 5 | 1.4496 | 1.4170 | 2.2% | 16.2% |
| H100 | 4096 | 10 | 2.3171 | 2.1466 | 7.4% | 38.0% |
| Blackwell | 512 | 5 | 0.7973 | 0.7634 | 4.3% | 31.7% |
| Blackwell | 512 | 10 | 1.5239 | 1.3614 | 10.7% | 55.3% |
| Blackwell | 4096 | 5 | 1.7549 | 1.6223 | 7.6% | 48.3% |
| Blackwell | 4096 | 10 | 3.1987 | 2.4726 | 22.7% | 72.4% |

H100's 512-atom degree-5 total is effectively unchanged. Do not substitute force
stage improvements for end-to-end improvements.

### 10.3 Horner: CPU with OpenMP disabled

| Compiler | Atoms | Degree | After/before, reference-normalized | Time reduction | Final / independent old CPU |
|---|---:|---:|---:|---:|---:|
| GNU | 512 | 5 | 0.779 | 22.1% | 0.640 |
| GNU | 512 | 10 | 0.831 | 16.9% | 0.874 |
| GNU | 4096 | 5 | 0.762 | 23.8% | 0.619 |
| GNU | 4096 | 10 | 0.724 | 27.6% | 0.876 |
| NVHPC | 512 | 5 | 0.809 | 19.1% | 0.616 |
| NVHPC | 512 | 10 | 0.915 | 8.5% | 0.979 |
| NVHPC | 4096 | 5 | 0.796 | 20.4% | 0.612 |
| NVHPC | 4096 | 10 | 0.844 | 15.6% | 1.027 |

Horner improved all measured CPU cases relative to the preceding common code.
NVHPC at 4096 atoms, degree 10 still trails the independent old CPU by about 3%.
This is not a guarantee of beating old CPU code for every model.
Ordinary CPU batches including per-call packing also improved at degree 5:
reference-normalized after/before ratios were GNU 0.772/0.778 and NVHPC 0.798/0.794
for 512/4096 atoms.

### 10.4 LAMMPS direct versus moment after Horner

H100, Ti/O 24,000 atoms, radial degree 20, angular degree 6, GPU neighbors, one rank;
20 warmup steps plus 100 measured steps, three samples with rotating mode order.

| Mode | Before Horner [ms/step] | After Horner [ms/step] |
|---|---:|---:|
| direct | 24.366 | 23.990 |
| moment | 9.754 | 8.925 |
| auto | 9.759 | 8.926 |

Moment is approximately **2.69×** faster than direct in this case. The small direct
change is not a Horner algebra benefit. All nine before/after final snapshots were
compared; the maximum force difference was `1.770e-13`.
This campaign did not remeasure the latest single-core CPU full-MD baseline, so
older CPU timings should not be presented as a current matched CPU/GPU speedup.

## 11. Implementation map, evidence, and remaining work

### 11.1 Source map

| Responsibility | Source |
|---|---|
| Geometry, moments, NN, Clenshaw/Horner, force scatter | [accelnet_target_kernels.f90](AccelNetPredictor/src/accelnet_target_kernels.f90) |
| Model packing, coefficient order, residency/lifetime | [accelnet_batch_target.f90](AccelNetPredictor/src/accelnet_batch_target.f90) |
| Shared CPU selection and independent reference | [accelnet_batch.f90](AccelNetPredictor/src/accelnet_batch.f90) |
| Generic descriptors and device scalar formulas | [accelnet_target_math.f90](AccelNetPredictor/src/accelnet_target_math.f90) |
| Shared scalar expressions | [shared](AccelNetDescriptors/src/shared/) |
| Chebyshev definitions and moment coefficients | [accelnet_descriptors.f90](AccelNetDescriptors/src/accelnet_descriptors.f90) |
| Serial module generation | [CMakeLists.txt](AccelNetPredictor/CMakeLists.txt) |
| LAMMPS CUDA neighbor adapter | [lib-gpu](interfaces/lammps/29Aug2024/lib-gpu/) |

### 11.2 Reproduction and validation records

- [OpenMP target build/comparison instructions](docs/openmp-target.md) and [CSR batch API](docs/batch-api.md).
- [Initial CPU batches](docs/validation/cpu-batch-2026-09-25/README.md).
- [Initial GPU](docs/validation/openmp-target-2026-09-25/README.md) and [residency/moments](docs/validation/gpu-residency-moments-2026-09-25/README.md).
- [LAMMPS/n2p2/CabanaMD design research](docs/lammps-n2p2-cabanamd-gpu-research.md) and [GPU package integration plan](docs/lammps-gpu-package-plan.md).
- [LAMMPS GPU validation](docs/validation/lammps-gpu-2026-09-26/README.md) and [initial direct/moment profiling](docs/validation/lammps-gpu-modes-2026-09-26/README.md).
- [Moment coefficient/gradient contraction](docs/validation/gpu-moment-force-2026-09-26/README.md).
- [G1–G5/LJ CPU diagnostics](docs/validation/gpu-descriptors-cpu-comparison-2026-09-26/README.md).
- [Chebyshev sharing/Clenshaw](docs/validation/chebyshev-common-2026-09-26/README.md) and [Horner validation](docs/validation/moment-horner-2026-09-26/README.md).

Use raw logs/JSON with their inputs, baselines, and timing boundaries. Horner
validation passed H100 GPU 29 + host 2, Blackwell GPU 29, GNU/NVHPC serial 10 each,
and existing CPU numerical 42 and performance 2 tests, among other checks.
These are recorded results, not tests rerun merely to write this document.

### 11.3 Further work

Grouped G4/G5 caching, LJ fusion, and value-only helpers have now been implemented
(Section 12). Remaining candidates include pair-level direct parallelism, removal
of directed-pair duplication, and reduced CSR/packing transfers. Computing both pair derivatives together can add
atomics, reductions, or scratch; fewer arithmetic operations alone do not prove a
speedup. High-degree moment work must also consider power-basis cancellation and
buffer growth; do not extrapolate tested accuracy to arbitrary degrees.

For subsequent changes, compare the independent CPU reference, common CPU with
OpenMP disabled, and actual GPUs under matched conditions. Update this document's
revision identity, adoption status, and measured results when adopting a change.

## 12. Revision 1.1: common LJ/Behler optimization

The CPU diagnostic ideas in Section 7 were implemented in the **same numerical
source used by serial CPU and OpenMP target GPU execution**. The default LJ/Behler
CPU API remains on the independent old CPU implementation because the optimized
common code still trails it for some families.

### 12.1 Adopted changes

G4/G5 features now share a radial cache when all of
$(R_c,\eta,R_s,\mathrm{cutoff},\alpha_{\rm cutoff})$ match exactly. Neither species
pair nor angular parameters change $q(r)$, and G4/G5 can share the same group.
With $E$ CSR edges and $G$ radial groups, values and derivatives require

$$
M_{\rm cache}=2EG\times8\ \mathrm{bytes}
$$

in FP64, rather than allocating per angular feature. The allocator retains at least
one group slot. Model packing stores group/representative indices, and the existing
geometry stage fills the device-resident cache once per call. No additional GPU
launch is needed; the value and force stages reuse it. Capacity growth, model reload,
and release follow the existing workspace ownership rules.

Value-only Behler helpers omit discarded derivatives. For G4, let
$\mathbf v_{jk}=(r_k\mathbf u_k-r_j\mathbf u_j)/r_{jk}$; its remaining pair derivative is

$$
\frac{\partial G^{\rm G4}_{jk}}{\partial\mathbf d_j}
=2\left[
\frac{A'q_jq_kq_{jk}}{r_j}(\mathbf u_k-c_{jk}\mathbf u_j)
+Aq_k(q'_jq_{jk}\mathbf u_j-q_jq'_{jk}\mathbf v_{jk})
\right].
$$

The $q_{jk}$ terms still require pair-dependent evaluation. Three directed/value
visits are not eliminated by this change.

LJ6/LJ12 values and forces are evaluated together using the formulas in Section 7.
Radial coefficients are contracted to a scalar before multiplying $\mathbf u_j$.
G1/G2/G3 also use the value-only path and scalar radial force contraction; the
mixed-Behler benchmark includes them, but does not isolate their individual speedups.

### 12.2 Measured common-code improvements

These are synchronous batches on the same two-species synthetic inputs used in the
[generic optimization report](docs/validation/generic-common-2026-09-26/README.md).
CPU ratios compare **new common / previous common**, using same-run reference ratios.
GPU ratios compare previous/new absolute batch times. They do not compare with
full CPU MD or imply universal model-independent speedups.

| Atoms | Family | GNU CPU new/old | NVHPC CPU new/old | H100 speedup | Blackwell speedup |
|---|---|---:|---:|---:|---:|
| 512 | lj | 0.609 | 0.538 | 1.05× | 1.18× |
| 512 | g4 | 0.471 | 0.593 | 1.77× | 2.11× |
| 512 | g5 | 0.362 | 0.550 | 1.91× | 2.76× |
| 512 | behler | 0.425 | 0.579 | 1.79× | 2.11× |
| 4096 | lj | 0.600 | 0.554 | 1.00× | 1.11× |
| 4096 | g4 | 0.466 | 0.592 | 1.52× | 1.89× |
| 4096 | g5 | 0.359 | 0.537 | 1.27× | 2.46× |
| 4096 | behler | 0.428 | 0.564 | 1.90× | 2.08× |

H100's 4096-atom LJ total is effectively unchanged despite smaller descriptor and
force phases; upload is about 0.66–0.68 ms of the roughly 1.20 ms total. G4/G5 and
mixed Behler improve substantially on both GPUs and both CPU compilers. The common
G4 CPU path still takes roughly 2.3–2.5 times the established CPU time in these
cases. This supports the optimization, but not replacing every CPU path yet.

The new group-sharing regression cases split and merge radial groups without
changing model dimensions, verify buffer growth/reuse, and compare energies,
forces, and virial including finite differences. See the linked report for the
complete validation log, binary/source identities, and reproduction commands.

### 12.3 Document revision history

| Document revision | Scope |
|---|---|
| 1.0 | CSR/residency, moment identities and force contraction, direct Clenshaw, moment Horner, CPU/GPU comparison methodology, and historical generic CPU diagnostics |
| 1.1 | Common G4/G5 radial caching, value-only Behler helpers, fused LJ, and matched CPU/H100/Blackwell measurements |
| 1.2 | CPU-style Jacobian experiments; Behler angular coefficient contraction, differentiated Horner, grouped values, and measured force-accumulation choices |

These are documentation revisions, not new AccelNet release numbers; the CMake
project version remains 1.0.1.

## 13. Reusing the Chebyshev lessons for Behler descriptors

The relevant Chebyshev optimizations are algebraic contraction, reuse of common
radial/geometry factors, evaluation of a polynomial and its derivative together,
and avoiding a large descriptor Jacobian. Clenshaw itself is specific to a
recurrence basis: Behler's angular power is more naturally evaluated by Horner in
its original nonnegative variable.

### 13.1 Exact angular groups and contracted coefficients

For a group sharing descriptor kind (G4 or G5), the unordered species pair,
radial parameters, cutoff settings, and lambda, define

$$
t=\frac{1+\lambda c}{2},\qquad
R_{jk}=q_jq_kq_{jk},\qquad
G_b=2\sum_{j<k}R_{jk}t_{jk}^{\zeta_b}.
$$

Here $q_{jk}=1$ for G5. The existing radial cache supplies $q_j,q'_j,q_k,q'_k$.
Integer powers up to the implementation limit (16) form one polynomial group;
other powers keep their original formula and may only combine identical zeta.
Polynomial grouping requires exactly integer zeta; near-integer values retain the
original derivative prefactor even if the scalar helper recognizes an integer power.
The limit bounds scratch/register usage; it does not truncate a series or
approximate higher powers.

After NN backpropagation, with $g_b=-\partial E/\partial G_b$, collect

$$
a_m=\sum_{b:\zeta_b=m}g_b,\qquad
H(t)=\sum_{m=1}^{p}a_mt^m.
$$

Initialize $h=a_p$, $d=0$ and, descending from $m=p-1$ to zero, use

$$
d\leftarrow td+h_{\rm old},\qquad
h\leftarrow th_{\rm old}+a_m.
$$

Then $h=H(t)$ and $H_c=(\lambda/2)d$. No conversion to powers of $c$ is needed;
$t$ stays in $[0,1]$ for the supported $|\lambda|\leq1$ and clipped cosine. The
force contribution for the whole group is

$$
\mathbf F_j^{\rm group}=2\left[
\frac{H_cR_{jk}}{r_j}(\mathbf u_k-c\mathbf u_j)
+Hq_k(q'_jq_{jk}\mathbf u_j-q_jq'_{jk}\mathbf v_{jk})
\right].
$$

For an exact noninteger-power group, $H=(\sum_bg_b)t^\zeta$ and
$H_c=(\sum_bg_b)(\lambda\zeta/2)t^{\zeta-1}$. Use the original power helper,
including its integer specializations for degrees above the polynomial limit.
This also keeps $t=0$, zeta=1 well defined without division by $t$.

Descriptor values can share the same pair geometry and radial factor. Accumulate
$S_m=2\sum_{j<k}R_{jk}t_{jk}^m$ locally using $t^m=t\,t^{m-1}$, then write
$G_b=S_{\zeta_b}$ once. This avoids repeated global output updates inside the pair
loop. Fractional-power groups use one scalar sum. Species maps and output indices
are preserved; grouping does not remove NN inputs or alter normalization.

### 13.2 Why copying the old CPU strategy is a useful experiment

The old CPU descriptor implementation evaluates an unordered pair once, stores
both derivatives, then contracts them with the NN gradient. The original common
G4/G5 implementation evaluated values once and revisited the pair for each directed
force. Radial caching alone did not remove these three visits.

A Jacobian trial changed the common implementation to the old CPU strategy. It
improved CPU performance substantially with **OpenMP compilation disabled**, and
showed that the CPU algorithm was a valid optimization candidate. Its global
scratch, however, scales as $3EB$ FP64 values, where $E$ is the CSR edge count.
The trial with GPU threads sharing a center/feature also needed atomic Jacobian
updates. The measured result depended on GPU and problem size: a small-system win
did not establish a large-system win. The archived phase timings locate the
large H100 regression in descriptor/Jacobian construction; they are not hardware
counter evidence proving a particular memory-bandwidth limit.

A second experiment forms both group forces after NN contraction:

$$
\mathbf F_k^{\rm group}=2\left[
\frac{H_cR_{jk}}{r_k}(\mathbf u_j-c\mathbf u_k)
+Hq_j(q'_kq_{jk}\mathbf u_k+q_kq'_{jk}\mathbf v_{jk})
\right].
$$

It visits $j<k$ once in the force stage and accumulates into the two edge-force
slots. CPU builds strip the OpenMP directives, while GPU edge threads require
atomic additions. This reduces arithmetic but introduces a different execution
cost. It must be measured against the independently owned directed-edge version,
which evaluates each direction separately without these atomics.

All CPU comparisons in this investigation use no `-fopenmp`/`-mp` compilation flags,
the same compiler and optimization options for both candidates, and affinity to
one core. `OMP_NUM_THREADS=1` alone is not the serial-build criterion. NVHPC may
link `libnvomp` indirectly through its normal Fortran runtime even without `-mp`;
that dependency alone does not mean the CPU kernels contain parallel regions.

The [Behler contraction report](docs/validation/behler-contraction-2026-09-26/README.md)
records the implementations, results, and adoption decision. Performance claims
must distinguish a faster common implementation from a win over the established
CPU reference, and distinguish an isolated synchronous batch from full LAMMPS MD.

G5's separable $q_jq_k$ factor permits a further multinomial-moment construction
for integer angular powers, analogous to Chebyshev. That is a distinct algorithm
and is not implemented by the contraction above. G4 contains the additional
$q_{jk}$ depending on both neighbor vectors, so the same independent-neighbor
moment identity cannot be applied unchanged. Neither fact prevents G4/G5 from
sharing radial caches, angular coefficient contraction, or differentiated Horner.

### 13.3 Adopted execution choices

Angular groups containing one distinct power retain the original value-only
pair kernel; a polynomial recurrence is unnecessary for such groups. Models with
multiple distinct integer powers in a group use local power sums. GPU teams share
neighbors within a center/group and reduce these sums, then write the NN inputs
once. With OpenMP disabled, the same source executes a serial loop. Equal-power coefficients
are summed in their first NN-gradient slot at the end of the existing NN kernel.
Distinct powers already have the required coefficient and stay in place. A bounded
local vector gathers these slots for Horner; there is no additional coefficient
array, workspace allocation, coefficient-formation launch, or host transfer.

The force kernel uses one common scalar formula and Horner implementation.
The host path visits unordered pairs once and updates both edge-force slots.
The GPU path gives each directed edge its own output slot, avoiding the atomic
force updates that slowed the tested two-sided GPU variant. This is an explicit
execution choice in the same source, not a claim that identical GPU and CPU
thread ownership always wins. The global descriptor-Jacobian trial is not adopted. Purely radial/LJ host batches
also keep directed-edge assignment, avoiding unnecessary zeroing and additions.

The default CPU LJ/Behler dispatch remains the independent established path.
Improving common code alone does not justify changing that default when some
models still run faster on the established CPU algorithm. The report separates
common-before/common-after comparisons from comparisons with that CPU reference.

### 13.4 Final measurements

The table uses 4096 atoms and the 48-input integer-power models defined in the
[report](docs/validation/behler-contraction-2026-09-26/README.md). CPU entries are
**final common / established CPU** time with OpenMP compilation disabled (less
than one is faster). GPU entries are **previous common / final common** speedups;
they are not GPU-versus-CPU speedups or full MD measurements.

| Model | GNU CPU ratio | NVHPC CPU ratio | H100 speedup | Blackwell speedup |
|---|---:|---:|---:|---:|
| g4-series | 1.048 | 0.957 | 2.09x | 4.19x |
| g5-series | 0.364 | 0.334 | 2.13x | 2.48x |

The smaller four-input G4 models still favor the established CPU path; the report
also includes those cases, G5, mixed Behler, and LJ at 512 and 4096 atoms. Hence
these numbers support algebraic reuse for compatible bases, not a universal
performance guarantee or an unconditional CPU-dispatch replacement. The final
99-case descriptor checks, device memcheck, and LAMMPS comparisons passed.


## 14. Reducing G4 pair work in the common kernel (revision 1.3)

These changes modify the same numerical source compiled for CPU and GPU. They
add neither a CPU-specific G4 implementation nor an additional device kernel,
workspace array, approximation, or CPU OpenMP region. The software/toolchain
versions in Section 1 are unchanged. Exact source identities, commands, and
measurements are recorded in the
[G4 validation report](docs/validation/g4-common-2026-09-27/README.md).

### 14.1 Equal-power value reuse and early pair rejection

Revision 1.2 combined identical angular NN-gradient coefficients but still
recomputed their descriptor values in the single-power value kernel. For a group
of identical unnormalized descriptor functions,

$$
G_{b_1}=G_{b_2}=\cdots=G_{b_m},
$$

one representative now evaluates the neighbor-pair sum and writes that value to
all group members. Each member keeps its own subsequent normalization and NN
weights. The existing group key includes descriptor kind, unordered species pair,
radial parameters, cutoff configuration, lambda, and the relevant power; merely
having the same angular exponent does not make two descriptors interchangeable.
The multi-power value kernel already performed the corresponding grouped work.

For G4, define $\mathbf d_{jk}=r_k\mathbf u_k-r_j\mathbf u_j$ and
$s_{jk}=\mathbf d_{jk}\cdot\mathbf d_{jk}$. Test

$$
s_{jk}\le\epsilon^2\quad\text{or}\quad s_{jk}\ge R_c^2
$$

before calculating $\sqrt{s_{jk}}$, a cutoff/exponential, or an angular power.
Such pairs make no G4 contribution under the existing cutoff convention. This
retains the distance test used by the established CPU implementation. It does
not substitute the less numerically stable difference
$r_j^2+r_k^2-2r_jr_k\cos\theta$ for the squared Cartesian displacement.

### 14.2 Scalar contraction before Cartesian components

Let $a(c)$ be the NN-gradient-contracted angular polynomial (or the single-power
fallback), $a'(c)=da/dc$, $c=\mathbf u_j\cdot\mathbf u_k$, and
$P=q_jq_kq_{jk}$. With $\mathbf u_{jk}$ directed from neighbor $j$ to $k$, define

$$
A_j=\frac{2a'P}{r_j},\qquad A_k=\frac{2a'P}{r_k},\qquad
B=2a q_jq_kq'_{jk},
$$

$$
R_j=2a q_kq'_jq_{jk}-c A_j,\qquad
R_k=2a q_jq'_kq_{jk}-c A_k.
$$

The same pair force contribution is then

$$
\mathbf f_j=A_j\mathbf u_k+R_j\mathbf u_j-B\mathbf u_{jk},\qquad
\mathbf f_k=A_k\mathbf u_j+R_k\mathbf u_k+B\mathbf u_{jk}.
$$

This explicitly reuses scalar products and performs scalar division before
forming the three Cartesian components. The coefficient $a$ already includes
the negative energy-gradient sign; no additional force-sign change is made.
For G5, $q_{jk}=1$ and $q'_{jk}=0$, so the same formula applies. CPU unordered
pair updates and GPU directed-edge ownership remain the revision 1.2 execution
policies around this shared expression.

The packed integer feature argument also has an explicit extent of 14, matching
the packer's layout, instead of an assumed-shape argument. This avoids passing
an array descriptor for this small fixed-layout argument on every pair call.
It is an internal interface change, not a change to the public C/Fortran batch API.

These are complementary changes. The incremental scalar-force rewrite gave
only about one percent in the initial alternating GNU comparison; it must not
be credited with the full improvement from value reuse, rejection, and explicit
argument shape. No default CPU LJ/Behler dispatch switch is implied by this
optimization: compare the common kernel against the established single-core CPU
path separately from comparing common-kernel revisions.

### 14.3 Measurements and remaining scope

At 4096 atoms, speedups relative to the **previous common implementation** were:

| Model | GNU CPU | NVHPC CPU | H100 | Blackwell |
|---|---:|---:|---:|---:|
| g4 | 1.29x | 1.30x | 1.09x | 1.11x |
| g4-series | 1.08x | 1.14x | 1.18x | 1.06x |

CPU measurements disable OpenMP compilation and pin execution to one core. The
four-input `g4` fixture includes two identical mixed-species inputs; `g4-series`
has 48 distinct inputs. Reuse of duplicate values explains part of the small
fixture's improvement, not the distinct-basis result.

Relative to the established CPU algorithm, the new common four-input `g4` kernel
still takes about 1.08x (GNU) or 1.03x (NVHPC) at 4096 atoms. For the new
`g4-distinct` four-input control without duplicates, the ratios are about 1.44x
and 1.28x in that campaign. A subsequent two-round GNU comparison using the
saved revision-1.2 numerical libraries found a 1.10x common-kernel speedup
even for this distinct basis, while its time relative to the legacy CPU varied
to about 1.58x; both campaigns are retained in the report. The larger `g4-series` common kernel takes about 0.71x and 0.70x.
Thus this revision improves common G4 execution but does not justify an
unconditional default CPU switch. The established CPU G4 algorithm's fused
value/Jacobian traversal remains advantageous for a small distinct basis;
the common value/contracted-force pipeline traverses pairs twice.

The final validation includes 130 descriptor cases, H100 and Blackwell tests,
zero H100 memcheck errors, 63 LAMMPS comparisons, and the ordinary CPU correctness
and performance gates. Exact timings, control experiments, hashes, and commands
are in the linked report. These batch speedups are not full LAMMPS MD speedups.


## 15. G4 value/Jacobian traversal shared by CPU and GPU (revision 1.4)

This revision implements the requested comparison against the established CPU
G4 algorithm. The immediate before baseline is the `gpu` checkpoint `cd00106`
(revision 1.3). Both CPU and GPU now use the same center-to-unordered-pair-to-feature
loop for G4, including independent cutoff, exponential, and angular reuse. This
is not the earlier feature-to-pair Jacobian experiment.

For center $i$ and descriptor $b$, accumulate its value and edge derivatives
in the same pair traversal:

$$
G_{ib}=\sum_{j<k}g_b(\mathbf r_{ij},\mathbf r_{ik}),\qquad
J_{\mu b e}=\frac{\partial G_{ib}}{\partial r_{ij,\mu}},\quad e=(i,j).
$$

Each pair adds both $\partial g_b/\partial\mathbf r_{ij}$ and
$\partial g_b/\partial\mathbf r_{ik}$ to their respective edge slots. The NN
then produces the scaled negative descriptor gradients $w_{ib}$, and the force
stage only contracts

$$
f_{\mu e}=\sum_b w_{ib}J_{\mu b e}.
$$

The existing scatter step adds neighbor forces, the opposite center force, and
the image-displacement virial. **No G4 pair geometry, radial factor, or angular
power is recalculated in the force stage.** Unlike revision 1.3, G4 keeps each
NN gradient separately; aggregating duplicate gradients and also reading every
Jacobian column would count duplicate contributions twice. G5 retains its
previous coefficient-contracted force algorithm.

Packed integer fields 15–19 identify G4 cutoff, exponential, and angular
representatives, the next descriptor for an unordered species pair, and the G4
radial representative. Cutoffs include cutoff type and alpha in their key;
exponentials include eta and shift; angular sharing is within a species pair
and includes lambda and the original real zeta. Integer/near-integer/fractional
power semantics remain unchanged. A named packed-field extent replaces literal
array extents in the internal helper interfaces.

The Jacobian uses the CPU layout `(Cartesian, descriptor, CSR edge)`. Each center
owns its edges and writes both sides of each pair without atomics. The GPU target
loop distributes centers; the OpenMP-disabled CPU build executes the same loop
serially. Pair caches and species-pair heads occupy persistent per-center
workspace arrays. These explicit mapped buffers avoid the invalid-address issue
observed with variable-sized private arrays in the first offload prototype.
There are no allocations inside the device kernel and no Jacobian/cache transfers
between host and GPU during evaluation.

For $E$ edges, $D$ input columns in species containing G4, $R$ centers, and $S$
local species, additional resident payload is approximately

$$
M_J=8\cdot3DE,\qquad M_{\rm scratch}=8\cdot12DR+4S^2R\quad\text{bytes}.
$$

The host also reserves corresponding arrays. Mixed models reserve all input
columns of G4-containing species, including unused non-G4 columns, to retain a
simple direct index. Capacities grow and persist; models without G4 need only
minimal placeholders. This memory cost is part of the one-traversal design and
must be reported alongside steady-state timings.

The [comparison report](docs/validation/g4-fused-2026-09-27/README.md) records
OpenMP-disabled single-core CPU measurements, H100/Blackwell timings, memory,
numerical checks, and source/binary identities. Earlier revision tables remain
historical measurements of the former two-traversal common implementation.


### 15.1 Measured outcome

At 4096 atoms, the timings below compare revision 1.3 against the shared
one-traversal implementation. Values are milliseconds per batch evaluation,
including transfers for GPU, but excluding neighbor construction. CPU builds
disable OpenMP compilation and use one core. These are two-round medians with
five samples per round, not full LAMMPS MD timings.

| Backend | 4 distinct inputs: before → fused ms | 48 distinct inputs: before → fused ms |
|---|---:|---:|
| gnu | 197.500 → 104.196 | 247.388 → 296.201 |
| nvhpc | 139.868 → 80.665 | 197.842 → 260.596 |
| h100 | 2.919 → 6.839 | 5.586 → 30.406 |
| blackwell | 9.943 → 14.309 | 11.330 → 34.829 |

The four-input common CPU implementation now beats the established CPU evaluator
as well as the previous common code. The 48-input result is less favorable, and
both GPUs regress in every measured G4 case. The force stage is cheaper, but the
descriptor/Jacobian stage is more expensive; the report includes phase timings.
The GPU implementation currently distributes only centers, with batch-resident
derivatives and per-center pair caches. These observations do not prove a
one-traversal algorithm inherently unsuitable for GPU; separating parallelism
from memory/cache costs requires further measurement.

The pair loop and derivative calculation now follow the established CPU strategy,
but the complete pipelines still differ: the CPU evaluator consumes each center's
Jacobian immediately, whereas the common backend keeps all centers' derivatives
across a batch-wide NN stage. This storage lifetime is a remaining optimization
opportunity. At 4096 atoms and 48 inputs, the Jacobian and pair-cache/head buffers
add about 148.5 MiB of resident payload, with corresponding host allocations.

Validation passed 135 descriptor cases, the GNU/NVHPC and H100/Blackwell suites,
63 LAMMPS comparisons, and ordinary CPU correctness/performance gates. H100
memcheck passed three selected benchmark smoke cases. The common one-traversal
source is retained, and a separate `lmp-g4-fused` executable preserves this
measured candidate alongside the previous installed `lmp`.

## 16. G4 scalar caches and flat descriptor ownership (revision 1.5)

Revision 1.4 proved that a shared value/Jacobian traversal could improve the small
CPU G4 models, but its GPU implementation regressed. Revision 1.5 retains the
same value/derivative mathematics and saved-Jacobian force contraction, and
changes the reuse and execution layout. The reference checkpoint is `2bcc603`;
the older two-traversal reference remains `cd00106` (revision 1.3).

### 16.1 Scalar pair reuse

For each matching species-pair list, retain the last cutoff, exponential, radial
group, and angular group identities and their computed values in scalar locals.
Recompute a factor when its identity changes. The radial group caches

$$
P=q_jq_kq_{jk},\qquad
\mathbf R_j=q'_j\mathbf u_jq_kq_{jk}-q_jq_kq'_{jk}\mathbf u_{jk},\qquad
\mathbf R_k=q_jq'_k\mathbf u_kq_{jk}+q_jq_kq'_{jk}\mathbf u_{jk}.
$$

With $A=[(1+\lambda c)/2]^\zeta$ and $A'=dA/dc$, each descriptor contributes

$$
g_b=2AP,\qquad
\nabla_jg_b=2A'P\,\frac{\mathbf u_k-c\mathbf u_j}{r_j}+2A\mathbf R_j,\qquad
\nabla_kg_b=2A'P\,\frac{\mathbf u_j-c\mathbf u_k}{r_k}+2A\mathbf R_k.
$$

The pair loop visits only its matching descriptor list; it no longer scans all
G4 descriptors to fill global cutoff/exponential/radial scratch for unrelated
species pairs. Reuse is consecutive within an owner. Interleaved groups are
recomputed, so arbitrary model ordering is supported without a stale cache.
Cutoff type/alpha, eta/shift, and the original real angular power remain separate
keys. A cutoff-excluded representative need not run before a later descriptor:
the current descriptor initializes its own required factor.

Removing `g4_scratch` saves $8\cdot12DN$ bytes of resident storage and the same
host reservation: **18 MiB at 4096 atoms and 48 inputs**. The saved Jacobian
$J_{\mu b e}$ remains; its storage is $8\cdot3DE$ bytes. This is not a
Jacobian-free implementation.

### 16.2 One owner per descriptor column

Flatten each species-pair descriptor list into packed positions $p=1,\ldots,D_4$.
Fields 19–21 store the descriptor index at a packed position and each list's
start/count; the packed integer extent becomes 21. For owner $\ell$ among $L$
owners, assign

$$
\ell=(p-1)\bmod L.
$$

Initialization and every pair visit use this same assignment. Each $(i,b)$
descriptor value and every edge derivative in its Jacobian column have one
writer. Consequently the G4 value/Jacobian kernel needs **no atomic additions**
and no synchronization inside the pair loop. The subsequent force/virial scatter
still uses its existing atomics. Models wider than $L$ assign multiple columns
to each owner; there is no descriptor-count limit of 32.

The flat work index is

$$
t=(i-1)L+\ell,\qquad t=0,\ldots,NL-1.
$$

One `target teams distribute parallel do` distributes these work items, with
`thread_limit(32)`. A small preceding kernel builds immutable species-pair heads.
The CPU compiled without OpenMP uses $L=1$. The GPU chooses a power of two up to
32, increasing it until it covers the largest matching G4 list and exposes at
least 16384 center/owner work items, or reaches the cap. This portable launch
heuristic was measured on H100 and Blackwell; it is not a promise of optimality
on every GPU or model. An OpenMP conditional source line/block changes only the
launch count. There is one numerical loop, with no separate CPU/GPU G4 formula.

GPU owners independently compute pair geometry and maintain their own scalar
caches. This intentionally trades some repeated geometry/radial work for more
parallelism and nearby descriptor-column writes. Each descriptor's value and
both derivatives are still computed together; the force stage never revisits
G4 pairs. The CPU retains the one-owner, pair-first sharing of geometry across
all matching descriptors.

### 16.3 Experiments that informed the layout

The [validation report](docs/validation/g4-tuning-2026-09-27/README.md) preserves
successful and rejected trials, commands, raw timings, compiler resource reports,
and Nsight Systems launch traces.

- A 32-thread limit improved the original center-only kernel over 64/128.
- Scalar caches removed global intermediate traffic and unrelated group work,
  improving both GPU and OpenMP-disabled CPU timings.
- A bounded private Jacobian increased GPU stack use to about 99 KB per thread
  and slowed the 48-input case; it was rejected. Declaring an array private does
  not guarantee register or shared-memory placement.
- Distributing neighbor pairs with atomic accumulation was slower in the tested
  nested implementation. The experiment does not isolate atomic cost from the
  nested execution overhead.
- An eight-center interleaved Jacobian layout did not improve GPU timing. It was
  rejected; the final Jacobian remains `(Cartesian, descriptor, CSR edge)`.
- Nested descriptor worksharing avoided atomics but generated 64-thread CUDA
  blocks despite a 32-thread OpenMP limit, with about 1.5 KB of stack per thread.
  Explicit team counts helped but were insufficient. Flat ownership removed the
  nested parallel region. Grid/block sizes were measured with Nsight Systems;
  stack/register counts come from `cuobjdump`, not inferred hardware counters.
- Nsight Compute counter collection was denied with `ERR_NVGPUCTRPERM`. No
  bandwidth, occupancy, or cache-hit measurements are claimed from that attempt.

### 16.4 Final measurements

The report compares revision 1.5 against **both** the one-traversal revision 1.4
and the older two-traversal revision 1.3. The established CPU evaluator is a third
reference, timed in the same single-core executable with OpenMP compilation off.
These baselines answer different questions and must not be conflated.

At 4096 atoms, the following table compares the starting one-traversal v1.4
against the final v1.5. Times are milliseconds per complete batch evaluation.

| Backend | 4 distinct inputs: v1.4 → v1.5 ms | 48 inputs: v1.4 → v1.5 ms |
|---|---:|---:|
| gnu | 104.175 → 77.746 | 296.399 → 213.038 |
| nvhpc | 80.511 → 59.906 | 259.779 → 204.053 |
| h100 | 6.840 → 3.420 | 30.262 → 4.620 |
| blackwell | 14.224 → 7.680 | 34.911 → 10.837 |

Every measured CPU/GPU case improved relative to v1.4. Final common CPU time
relative to the established CPU evaluator ranges from **0.592 to
0.853** across the tested compilers, sizes and models.
GPU speedups over v1.4 range from **1.85x to 17.62x**.


Some comparisons against v1.3 remain slower, notably the tested H100 four-input
4096-atom case; consult the full table rather than interpreting the v1.4 speedup
as a universal win over the older two-traversal algorithm. Numerical validation
passed 137 descriptor cases, all stated CPU/GPU suites, three selected memcheck
cases, 63 LAMMPS comparisons, and the ordinary CPU correctness/performance gates.
The validated `lmp-g4-optimized` candidate is kept alongside the previous binaries.

## 17. Default common CPU batch dispatch (revision 1.6)

The G5 dispatch limitation described in this historical section is superseded
by revision 1.7 below.

Revision 1.6 adopts the shared numerical implementation for every CSR batch
model accepted by the existing target-model packer. In addition to the already
shared Chebyshev direct/moment path, this includes LJ, Behler G1–G5, multiple
LJ/Behler components, and supported per-element combinations. The CPU uses the
serial module instance generated from the GPU source: OpenMP directives are
removed at configure time and no OpenMP runtime or GPU compiler is required.
This dispatch change introduces no new descriptor formula.

The public `evaluate_batch` repacks the model on each call to honor edits to its
public fields. It retains the workspace, including the G4 Jacobian. Therefore
its measured cost is

$$
T_{\mathrm{batch}} = T_{\mathrm{pack}} + T_{\mathrm{prepare}}
 + T_{\mathrm{descriptor}} + T_{\mathrm{NN}} + T_{\mathrm{force/virial}}.
$$

Unlike prepared-`target_model` timings, these measurements include
$T_{\mathrm{pack}}$. Fixed-neighbor timings exclude neighbor construction;
the new end-to-end regression gate includes it for both implementations.
Initial allocations are warmed up in either case. G4 retains the full
$(3,\mathrm{descriptor},\mathrm{edge})$ Jacobian, so its storage grows with the
batch's edge count, unlike the reference's per-center temporary Jacobian.

Unsupported packing falls back to `evaluate_batch_reference`. In particular,
mixed/multiple Chebyshev components within one element and explicitly requested
G5 moments keep their established CPU behavior. G5 auto on the common path
uses direct pairs; the reference's auto policy may select moments. The
benchmark driver now allows generic mode 0 as well as mode 1 so that this
policy difference is included in the comparison. The measured cases do not
establish the optimal choice for arbitrary models or densities.

The adoption concerns the CSR batch API. Object/atomic APIs, CLI and ordinary
LAMMPS `pair_style accelnet` keep their existing evaluators. The old CPU code
also remains an independent numerical/performance reference, rather than being
deleted and then compared with itself. Earlier sections describe the adoption
status at their respective revisions.

The [revision-1.6 validation report](docs/validation/common-cpu-default-2026-09-27/README.md)
records GNU/NVHPC single-core timings with OpenMP compilation disabled,
8/512/4096-atom fixtures, two reversed rounds of five paired samples, and denser
G5 cases. It also separates the pre-existing high-order Chebyshev direct
slowdown from this dispatch change. Default CPU energy/force/virial checks were
added to the target descriptor suite; batch tests cover shared/fallback mode
transitions, model reloads and partitioned/additive output. The opt-in
`ACCELNET_TEST_COMMON_CPU_PERFORMANCE` CTest gate compares default batches with
the independent structure evaluator, including real n2p2 G4/G5 models.


## 18. Common G5 moments without losing the original advantage (revision 1.7)

The earlier degree-8, short-cutoff synthetic tests do not answer whether the
original G5 moment advantage survives. The original scaling benchmark uses
orders $\zeta\in\{1,2,4\}$, three radial groups and controlled neighbor counts.
Revision 1.7 reproduces those descriptor parameters and local environments, then
compares **old direct, old moment, common direct and common moment** with the
same network. CPU comparisons compile OpenMP out entirely. The
[validation archive](docs/validation/g5-moments-2026-09-27/README.md) records all
four methods, both CPU compilers, H100 and Blackwell.

### 18.1 Exact species-resolved moment contraction

For radial group $g=(\eta,R_s,R_c,\text{cutoff type},\alpha)$, define

$$
h_j^{g}=e^{-\eta(r_j-R_s)^2}f_c(r_j),\qquad
M_{\boldsymbol\beta}^{g,t}=\sum_{j:s_j=t}h_j^{g}\,\mathbf u_j^{\boldsymbol\beta},
\qquad S^{g,t}=\sum_{j:s_j=t}(h_j^{g})^2.
$$

Here $\boldsymbol\beta=(\beta_x,\beta_y,\beta_z)$ is a nonnegative multi-index,
$|\boldsymbol\beta|=q$ and
$\mathbf u^{\boldsymbol\beta}=u_x^{\beta_x}u_y^{\beta_y}u_z^{\beta_z}$.
For descriptor $b$ with integer $\zeta_b$,

$$
a_{bq}=2^{1-\zeta_b}{\zeta_b\choose q}\lambda_b^q,\qquad
G_b=\sum_{q=0}^{\zeta_b}a_{bq}B_q^{g,t_1,t_2},
$$

$$
B_q^{g,a,b}=
\begin{cases}
\displaystyle\sum_{|\boldsymbol\beta|=q}{q!\over\boldsymbol\beta!}
M_{\boldsymbol\beta}^{g,a}M_{\boldsymbol\beta}^{g,b},&a\ne b,\\[4pt]
\displaystyle\frac12\left(\sum_{|\boldsymbol\beta|=q}{q!\over\boldsymbol\beta!}
(M_{\boldsymbol\beta}^{g,a})^2-S^{g,a}\right),&a=b.
\end{cases}
$$

The same-species term removes self pairs and counts each unordered pair once.
It follows from $\mathbf u_j\cdot\mathbf u_j=1$; a self correction is required at
every degree, including zero. Distinct species need neither the factor $1/2$
nor a self correction. No neighbor-pair approximation is introduced.

Raw moments are built once and kept through NN evaluation. The bilinear
$B_q$ terms are computed **once per radial group, species pair and degree**,
then reused by every descriptor with that radial group. Repeating that reduction
for every descriptor would waste much of the moment method's benefit.

### 18.2 NN contraction before monomial and neighbor loops

Let $w_b=-\partial E/\partial G_b$ include the input-normalization derivative.
Accumulate the symmetric species-pair coefficients

$$
C_q^{g,a,b}=\sum_{d\text{ matching }g,\{a,b\}}w_d a_{dq}.
$$

The raw-moment adjoint and self coefficient are

$$
A_{\boldsymbol\beta}^{g,a}={|\boldsymbol\beta|!\over\boldsymbol\beta!}
\sum_b C_{|\boldsymbol\beta|}^{g,a,b}M_{\boldsymbol\beta}^{g,b},\qquad
d^{g,a}=\sum_q C_q^{g,a,a}.
$$

Thus each group's contribution for neighbor $j$ of species $a$ is

$$
\mathbf F_j^g=
(h_j^g)'\,[P^{g,a}(\mathbf u_j)-d^{g,a}h_j^g]\,\mathbf u_j
+{h_j^g\over r_j}
\left[\nabla P^{g,a}(\mathbf u_j)
-\mathbf u_j\big(\mathbf u_j\cdot\nabla P^{g,a}(\mathbf u_j)\big)\right],
\qquad
P^{g,a}(\mathbf u)=\sum_{\boldsymbol\beta}A_{\boldsymbol\beta}^{g,a}
\mathbf u^{\boldsymbol\beta}.
$$

These are edge-force contributions; the existing scatter accounts for the
central atom and virial. $P$ and its three derivatives are evaluated together by
nested differentiated Horner, reusing the Chebyshev moment optimization. A step
$H\leftarrow xH+c$ carries $D\leftarrow xD+H_{\rm old}$. This does not divide by
coordinates, so zero components are safe. Descriptor work ends at the degree
contraction; neither the monomial adjoint nor force loop scans all descriptors.

With $K={p+3\choose3}$ monomials through degree $p$, $R$ radial groups, $S$
species and $D$ descriptors, the main costs are moment construction
$O(RSKN_n)$ in the current species-channel loops, bilinear contractions
$O(RS^2K)$, descriptor/NN coefficient work $O(RS^2Dp)$ in the current matching
loops, and forces $O(RKN_n)$. There is no $N_n^2$ pair loop for active moments.
The padded maximum degree across groups is retained; group-specific degree
bounds and skipping unused species-pair channels remain possible optimizations.

### 18.3 Shared implementation, selection and storage

CPU and GPU compile the same moment construction, degree contractions, adjoints
and differentiated Horner source. GPU directives distribute disjoint
center/channel/monomial work; CPU compilation removes those directives. Raw
moments and adjoints occupy separate halves of the persistent moment buffer;
self corrections use an extra entry. Degree coefficients reuse the NN delta
buffer. Radial values/derivatives and unit-coordinate powers are cached. Rows
containing only active moments skip the irrelevant direct feature scans.

G5 selection is independent of Chebyshev selection. Auto (0) and thresholded
moment (2) use moments at 16 valid neighbors inside **each component's maximum
angular cutoff**. Direct (1) never uses moments. Forced moment (3) bypasses the
neighbor threshold. All moment modes require exact integer $1\le\zeta\le10$;
fractional, near-integer and higher powers retain the shared direct formula.
Mixed modes, radial groups, cutoffs and element models can coexist. The original
CPU evaluator remains an independent reference; G5 moment requests no longer
need its fallback.

`target_model%initialize` adds optional `g5_mode`. The C API adds
`accelnet_target_create_modes` without changing the old creation ABI. LAMMPS uses
`pair_style accelnet/gpu auto ... g5 moment`; the first mode controls Chebyshev,
and the trailing mode controls G5. Auto retains the original neighbor threshold;
that threshold is a compatibility policy, not a promise of optimal performance
for every degree, density, compiler or GPU batch size.

### 18.4 Validation and performance regression protection

The new 185-case suite compares energy, every force component and virial against
independent original CPU direct evaluation. It covers ten cutoffs, all four
modes, finite differences, integer/noninteger/high orders, collinear and zero
components, isolated atoms, periodic images, mixed components and central-element
models, per-component cutoff thresholds, and workspace reuse. GNU, NVHPC,
checked GNU, H100 and Blackwell validation passed. H100 memcheck found no errors;
63 LAMMPS comparisons include forced G5 moments against CPU direct, multiple
neighbor modes, MPI/empty ranks, triclinic cells and short trajectories.

The opt-in `ACCELNET_TEST_G5_MOMENT_PERFORMANCE` CTest requires the serial target
benchmark. Its controlled 64-center/64-neighbor degree-4 fixture measures all
four paths with alternating method order. It requires common moment time to be
at most 1.10 times original moment time and no slower than common direct. It
checks numerical agreement and rejects a binary importing known OpenMP runtime
symbols. The general default-CPU performance gate remains separate.


On the original degree-4 scaling fixture at 512 centers and 64 neighbors, the
final times (ms per complete fixed-neighbor evaluation) are:

| Backend | Old direct | Old moment | Common direct | Common moment |
|---|---:|---:|---:|---:|
| GNU single core, OpenMP off | 189.614 | 145.568 | 171.331 | 22.418 |
| NVHPC single core, OpenMP off | 165.276 | 149.391 | 151.821 | 18.684 |
| H100 | — | — | 4.412 | 1.937 |
| Blackwell | — | — | 6.667 | 2.384 |

The shared CPU timing includes per-call model packing; GPU timing uses a prepared
model and includes transfers. Neighbor construction is excluded. Original
moment beats original direct here, and common moment improves further. The
separate degree-8/small-cutoff H100 control gives 1.345 ms direct versus 1.361 ms
moment: almost equal, with moment slightly slower. Neither the original threshold
nor this optimization guarantees a moment win for every workload. Full tables,
measurement scope, ordinary direct controls, and regression gates are in the
revision-1.7 validation archive.


## 19. n2p2 weighted/compact descriptors and multi-element validation (revision 1.8)

This revision adds n2p2 types 12, 13, and 20--25 to the Fortran importer,
shared CPU/OpenMP-target evaluation, and Fortran conversion in both directions.
It builds on `gpu` checkpoint `f8c19a0` (revision 1.7). The independent oracle is
n2p2 **v2.3.0**, built with GNU C++ **11.4.0**, `-O3`, Eigen **3.4.0**, and no
OpenMP. The CPU candidate is GNU Fortran **11.4.0**, `-O3`, with OpenMP compiled
out. GPU builds use NVHPC **25.3** and the flags in Section 1.1.

### 19.1 Exact formulas and moment applicability

Let $Z_j$ denote atomic number, $r_j=|\mathbf r_{ij}|$,
$\mathbf u_j=\mathbf r_{ij}/r_j$, $c=\mathbf u_j\cdot\mathbf u_k$, and
$\theta=\arccos c$. Species-resolved descriptors restrict the neighbor sum;
weighted descriptors sum all species with $Z_j$ or $Z_jZ_k$.

For type 12, $q(r)=\exp[-\eta(r-R_s)^2]f_c(r)$ and
$G_{12}=\sum_j Z_j q(r_j)$. Type 13 is

$$G_{13}=\sum_{j<k} Z_j Z_k\,2^{1-\zeta}(1+\lambda c)^\zeta
q(r_j)q(r_k)q(r_{jk}).$$

For compact types, define

$$t=\left|\frac{x-(L+R)/2}{(R-L)/2}\right|,\qquad
W(x;L,R)=\begin{cases}P(t),&L<x<R,\\0,&\text{otherwise}.\end{cases}$$

The supported cores are

$$\begin{aligned}
P_1(t)&=1-3t^2+2t^3,\\
P_2(t)&=1-10t^3+15t^4-6t^5,\\
P_3(t)&=1-35t^4+84t^5-70t^6+20t^7,\\
P_4(t)&=1-126t^5+420t^6-540t^7+315t^8-70t^9,\\
P_e(t)&=\exp\!\left(1+\frac{1}{t^2-1}\right).
\end{aligned}$$

Asymmetric radial subtypes replace $t$ by $t(2-t)$ and apply its derivative
$2(1-t)$. Angular windows retain the symmetric core. Put
$q(r)=W(r;R_{\rm low},R_c)$ and $A(c)=W(\arccos c;\theta_L,\theta_R)$.
Types 20/23 are radial sums of $q$. Types 21/24 sum
$A(c)q(r_j)q(r_k)q(r_{jk})$, and types 22/25 omit $q(r_{jk})$.
Types 23--25 apply the atomic-number weights above.

An exact finite G5 moment contraction requires a finite polynomial in $c$
and separable radial factors. The G5/type-9 integer-$\zeta$ direct and moment
paths remain implemented and tested. Types 12/20/23 already cost $O(N_n)$;
there is no pair sum for moments to remove. Types 13/21/24 depend on $r_{jk}$.
Types 22/25 contain a window polynomial in **angle**, $\arccos c$, with compact
support; it is not a finite polynomial in $c$. Thus these general families
have no exact finite moment implementation analogous to G5. A truncated
expansion would change the potential and is not used.

### 19.2 Value and derivative in one pair traversal

For $V=A(c)q_jq_kq_{jk}$, put $P=q_jq_kq_{jk}$ and
$\mathbf u_{jk}=(\mathbf r_{ik}-\mathbf r_{ij})/r_{jk}$. Then

$$\nabla_{\mathbf r_{ij}}V=
 A'(c)P\frac{\mathbf u_k-c\mathbf u_j}{r_j}
 +A(c)q_k\left(q'_jq_{jk}\mathbf u_j-q_jq'_{jk}\mathbf u_{jk}\right),$$

$$\nabla_{\mathbf r_{ik}}V=
 A'(c)P\frac{\mathbf u_j-c\mathbf u_k}{r_k}
 +A(c)q_j\left(q'_kq_{jk}\mathbf u_k+q_kq'_{jk}\mathbf u_{jk}\right).$$

For wide types, set $q_{jk}=1$, $q'_{jk}=0$. Compact angular differentiation
uses $A'(c)=-W'(\theta)/\sqrt{1-c^2}$. The implementation caches $q_j,q'_j$
by radial-parameter group, evaluates value and both derivatives together, and
contracts the retained Jacobian after NN backpropagation. A center/descriptor
owns its Jacobian column, so pair accumulation needs no atomics. CPU and GPU
share the same code; OpenMP changes execution placement only.

This layout repeats some pair geometry across descriptor channels. n2p2's
grouped evaluator can reuse that work more widely, so shared angular CPU code
is not uniformly faster than n2p2. Reusing pair geometry while retaining enough
GPU parallelism is a concrete further optimization; no speedup is inferred
merely from removing transcendental radial functions.

### 19.3 Multi-element bugs caught by real models

Official Anisole has 354/351/331 descriptors for C/H/O, while DMABN has
333/334/219 for C/H/N. The tests caught two integration issues:

* The converter treated differing **input dimensions** as differing hidden
  topology, emitting redundant per-element overrides that upstream rejected.
  Hidden topology comparison now excludes the descriptor dimension.
* A minimum-distance packing bound could request tens of millions of neighbor
  slots per atom. A second bound uses actual atom count and cell geometry:
  $N\prod_a[\lceil2(R_c+\delta)\|L^{-1}_{a,:}\|\rceil+1]$, with the existing
  neighbor skin $\delta$. Dynamic CSR storage checks actual integer capacity.
  This changes allocation bounds, not selected neighbors.

The importer preserves n2p2 descriptor order and the associated scaling and
first-layer weights. New native metadata separates type-specific parameters
from global cutoff metadata. Force and strain finite differences, converted
models, two/four-element mixtures, and real trained models exercise the full
path. The known n2p2 collinear-window discontinuity is documented explicitly in
[model compatibility](docs/model-compatibility.md#10-weightedcompact-validation-methods-revision-18).

Measurements, source hashes, commands, oracle model hashes, and numerical
errors are in [the revision-1.8 archive](docs/validation/n2p2-extensions-2026-09-27/README.md).

## 20. Grouped extended angular evaluation and LAMMPS CPU batches (revision 1.9)

The baseline for this revision is `gpu` commit **18d1ac3** (library version
**1.0.1**). The measurement protocol and before/after results are archived in
[the CPU/LAMMPS validation report](docs/validation/n2p2-cpu-2026-09-27/README.md).
Both CPU inference binaries are compiled without OpenMP. The GPU implementation
uses the same Fortran source with OpenMP target directives enabled.

### 20.1 Reuse geometry before evaluating descriptor members

For each center, the previous extension kernel traversed its unordered neighbor
pairs once **per descriptor**. It already evaluated each descriptor's value and
Jacobian together, but repeated the cosine, Cartesian angle derivatives,
`acos`, and neighbor-neighbor distance for different descriptors. The G4 scalar
cache/ownership strategy in Section 16 applies to these extensions as well.

Let $g$ index identical radial parameters, $a$ identical angular parameters,
and $h\in\{0,1\}$ select a wide or narrow descriptor. For each pair, cache

$$P_{g,h}=q_g(r_j)q_g(r_k)\,[q_g(r_{jk})]^h,$$

$$\mathbf R_{j,g,h}=q'_g(r_j)q_g(r_k)[q_g(r_{jk})]^h\mathbf u_j
 -h q_g(r_j)q_g(r_k)q'_g(r_{jk})\mathbf u_{jk},$$

$$\mathbf R_{k,g,h}=q_g(r_j)q'_g(r_k)[q_g(r_{jk})]^h\mathbf u_k
 +h q_g(r_j)q_g(r_k)q'_g(r_{jk})\mathbf u_{jk}.$$

A descriptor member then adds $w A_a(c)P_{g,h}$ to its value, and

$$w\left[A'_a(c)P_{g,h}\frac{\mathbf u_k-c\mathbf u_j}{r_j}
       +A_a(c)\mathbf R_{j,g,h}\right]$$

to its $j$ derivative (analogously for $k$). Here $w=1$ for species-selected
compact functions and $w=Z_jZ_k$ for weighted functions. These scalar products
are evaluated before expanding Cartesian components, following the G4 lesson.

The packer builds an auxiliary member order by unordered species pair, angular
representative, narrow/wide flag, and radial group. It **does not reorder**
descriptor outputs, scaling arrays, or NN inputs. A weighted list is separate
from the species-pair list. Each constant-angular run evaluates its angular window once. If both its value
and derivative vanish, the whole run is skipped. The narrow path also caches
the last radial product and cutoff in scalars; no $O(N_n^2)$ geometry buffer is
introduced. Compact radial
asymmetry and angular symmetry use distinct keys: the angular key contains the
base polynomial subtype, while the radial key retains the full subtype.

A serial CPU owner visits a pair once and processes its matching members. GPU
owners partition the ordered descriptor indices into disjoint lanes. Every
Jacobian column has a single writer, avoiding atomics. Geometry is reused among
the members owned by each lane; it can be repeated across lanes in exchange for
more parallel work. The value and derivative formulas remain shared. GPU owners also retain the
union of their radial supports and a species-pair mask, allowing irrelevant
pairs to be rejected before angle evaluation. The mask accelerates the first
ten local species; higher indices retain the unrestricted, correct path.
The owner count is

$$L=\min\left(256,N_{\mathrm{angular,max}},
              \max(1,\lfloor65536/N_{\mathrm{rows}}\rfloor)\right).$$

The serial build always uses $L=1$. Limiting GPU execution to 32 owners left
small, descriptor-rich models underoccupied; additional owners plus early
rejection restored and improved their GPU throughput.

For wide types 22/25, let $\mathbf d_j=r_j\mathbf u_j$ and
$\mathbf d_k=r_k\mathbf u_k$. Precompute the angular coefficients once per run:

$$a_j=\frac{wA}{r_j},\quad a_k=\frac{wA}{r_k},\quad
 b_j=\frac{wA'c}{r_j^2},\quad b_k=\frac{wA'c}{r_k^2},\quad
 b_{jk}=\frac{wA'}{r_jr_k}.$$

For each radial member, $P=q_jq_k$ and

$$p_j=a_jq'_jq_k-b_jP,\quad
 p_k=a_kq_jq'_k-b_kP,\quad p_c=b_{jk}P,$$

$$\nabla_{\mathbf d_j}V=p_j\mathbf d_j+p_c\mathbf d_k,\qquad
  \nabla_{\mathbf d_k}V=p_k\mathbf d_k+p_c\mathbf d_j.$$

Compact angular windows additionally pack their center $m=(\theta_l+\theta_r)/2$
and inverse half-width $h=2/(\theta_r-\theta_l)$ at model load. Their normalized
argument is $y=(\theta-m)h$, avoiding repeated divisions in the pair loop.

This removes repeated angular scaling and narrow-only distance/cache branches
from the wide inner loop. It retains one pair pass for values and derivatives,
as in the optimized CPU G4 implementation. A trial array of radial products
increased CPU time and was removed: saving arithmetic does not guarantee a win
when it adds indexing, loads, and cache traffic.

### 20.2 Exact powers for weighted angular type 13

With $t=(1+\lambda c)/2$, the weighted angular factor is
$A=2t^\zeta$ and $A'=\lambda\zeta t^{\zeta-1}$.
Exact integer exponents 1--16 reuse the G4/G5 integer-power helper. For other
valid exponents, compute $p=t^{\zeta-1}$ once, then $A=2tp$ and
$A'=\lambda\zeta p$. Near-integer exponents are never rounded. This reduces
power evaluations without a polynomial approximation. It does not create a
finite exact moment representation for the nonseparable narrow radial term.

The same scalar Cartesian contraction now applies to G4. For a narrow member,
define

$$C=\frac{q_jq_kq'_{jk}}{r_{jk}},\quad
 R_j=\frac{q'_jq_kq_{jk}}{r_j}+C,\quad
 R_k=\frac{q_jq'_kq_{jk}}{r_k}+C.$$

With $P=q_jq_kq_{jk}$, the coefficients are

$$p_j=wA R_j-wA'P\frac{c}{r_j^2},\quad
 p_k=wA R_k-wA'P\frac{c}{r_k^2},\quad
 p_c=wA'P\frac{1}{r_jr_k}-wA C.$$

The Cartesian derivatives again equal $p_j\mathbf d_j+p_c\mathbf d_k$ and
$p_k\mathbf d_k+p_c\mathbf d_j$. Cache $R_j,R_k,C$ by radial group before
processing angular members. This removes repeated vector normalizations and
vector-valued radial intermediates from the pair loop.

Value/derivative evaluation also shares the hyperbolic tangent cutoff:

$$t=\tanh(1-r/R_c),\qquad f_c=t^3,\qquad
 f'_c=3t^2(t^2-1)/R_c.$$

For normalized tanh, divide both outputs by $\tanh^3(1)$. A single `tanh`
evaluation supplies both outputs; the cutoff boundary and all other cutoff
formulas are unchanged. This matters for the published water model, whose G4
functions use the unnormalized tanh cutoff.

G2 uses exactly the same $q(r)=\exp[-\eta(r-R_s)^2]f_c(r)$ as the radial
factors of G4/G5. The packed radial cache therefore includes G2, sharing
identical $(\eta,R_s,R_c,\mathrm{cutoff},\alpha)$ groups across chemical
species and angular functions. Both its value pass and its force pass read
$q,q'$ instead of recomputing the exponential and cutoff independently.
Angular polynomial grouping still applies only to G4/G5; it must not merge
G2 output slots into an angular contraction group. If a radial group contains
only G2 functions for one neighbor species, cache preparation skips all other
species. Unrestricted caching would perform extra work for these small models.

When all standard descriptors for a central species have the same cutoff,
geometry preparation also stores $f_c(r),f'_c(r)$ once per edge. Different
$\eta,R_s$ groups then apply their Gaussian factor to these cached cutoffs.
Compact or mixed-cutoff models retain their general path; exact metadata
comparison selects the optimization without changing any cutoff support.

The G4 angular derivative prefactor $\lambda\zeta/2$ is packed once. Its
factor of two is applied when the angular cache changes. Hot angular-power and
paired-cutoff helpers, and compact angular/window helpers, include the same
scalar source in the kernel translation
unit, allowing GNU Fortran to inline them without whole-program LTO. The
Jacobian's Cartesian extent is explicitly **3**, matching its allocation;
this exposes fixed-size stores and strides to both compilers.

### 20.3 Remove repeated neighbor candidate searches

The standalone orthogonal linked-cell builder previously checked candidate
membership with `any(nblist(1:count)==candidate)`. For $K$ visited candidates,
this can require $O(K^2)$ comparisons per center. Track visited wrapped cell
indices instead: each cell list is traversed once. With $C$ candidate cells,
the work is $O(K+C^2)$ per center, without resetting an $N$-atom membership
array for every center. A repeated central cell inserts the central atom once
at its original first-revisit position. Candidates retain their original
first-visit order; periodic
image expansion and distance filtering are unchanged. The independent neighbor
test enumerates all relevant images, including small cells, self images,
unwrapped coordinates, and multiple linked-cell grid shapes.

This optimization affects AccelNet's standalone neighbor builder. **LAMMPS
constructs its own neighbor lists**, so it cannot account for a LAMMPS speedup.

### 20.4 Connect the LAMMPS CPU interface to the shared batch implementation

The previous `pair_style accelnet` called the retained atomic C API separately
for every center. Consequently, changing the shared batch kernel alone did not
change that LAMMPS execution path. The new CSR C entry point uses the loaded
model and reuses a serial batch workspace. It preserves additive local/ghost
forces and one energy per center; finalization releases the workspace.

The LAMMPS adapter packs full, image-aware neighbor rows in chunks of at most
4 centers, bounding Jacobian scratch space while retaining the common CPU/GPU
numerical implementation. Edges in the neighbor-list skin beyond the model's
maximum cutoff are excluded while packing; their descriptor contribution is
zero. Each chunk remaps only its touched center/neighbor targets, then scatters
its forces back to LAMMPS local/ghost slots. The target map is initialized once
per step and only touched entries are reset after each chunk. This avoids
scanning all $N$ force targets once per chunk, an $O(N^2/B)$ overhead for fixed
chunk size $B$, while preserving ordinary Newton reverse communication.

The C API retains packed metadata for its private
loaded model. Every successful model load and evaluation-mode setter invalidates
this cache; finalization releases both packed metadata and scratch storage.
The public Fortran batch API still repacks metadata because callers can edit
its model directly. Smaller chunks improve locality without repeating the
relatively expensive metadata preparation. Normal LAMMPS Newton
reverse communication and energy/virial tallying remain in place. GPU inference
continues through the existing independent target-model interface.

### 20.5 Compare the actual LAMMPS force path

Standalone n2p2 and LAMMPS do not have identical force-assembly costs.
`Mode::calculateForces` in the standalone oracle searches neighboring centers'
neighbor lists, whereas `InterfaceLammps::getForces` scatters their stored
contributions directly to local/ghost force slots. Matching standalone timing
therefore does not establish LAMMPS parity. This revision measures both paths
and uses the actual LAMMPS integration loop for the MD comparison.

LAMMPS measurements use one CPU core with OpenMP compiled out, unmodified
published models, matching neighbor settings, and the same short NVE
trajectory. Initial/final forces, positions, energy, and six virial-pressure
components must agree before timings are accepted. GPU timings include the
LAMMPS GPU neighbor path and host/device transfers. They use the same model
converted to native format without retraining or descriptor changes. Detailed
timing samples and build versions are retained in the validation archive.


## 21. Exact high-order G5 moments (revision 1.10)

Baseline: `gpu` commit **2a96dcc**, AccelNet library **1.0.1**. This revision
extends the shared G5/type-9 moment evaluator from exact integer orders 1--10
to **1--16 in explicit moment modes**. It does not approximate any descriptor.
[Validation and CPU/H100 timings](docs/validation/high-g5-moments-2026-09-27/README.md)
record the supported range, fallbacks, and workloads where moments are slower.

### 21.1 Which existing descriptors admit this reduction?

| Descriptor | Exact fixed-order moment reduction |
|---|---|
| Chebyshev angular | Already implemented: finite polynomial in the cosine, separable radial weights |
| G5 / n2p2 type 9, integer zeta | Implemented; this revision adds explicit orders 11--16 |
| G1/G2/G3, LJ, n2p2 types 2/12/20/23 | Already single-neighbor sums, O(Nn); no pair sum to eliminate |
| G4 / n2p2 type 3, types 13/21/24 | General neighbor-neighbor radial dependence prevents the separable G5 reduction |
| n2p2 types 22/25 | Compact window in the angle, not a finite polynomial in its cosine |
| Fractional G5 zeta | Not a finite polynomial; direct is retained |

Here "fixed-order" means that the number of moments is independent of the
number of neighbors. For types 22/25, the angular factor is
$A(c)=W(\arccos c)$. A nonzero finite polynomial in $c$ cannot vanish on an
open interval of compact angular support. Even a window covering the full
physical angular range generally retains nonpolynomial $\arccos c$ dependence.
A truncated Chebyshev/Legendre expansion would change the potential and its
forces; it is deliberately excluded from this exact-equivalence revision.

### 21.2 Extend the finite polynomial, not the descriptor definition

For integer $p=\zeta$, the same identity used in Section 18 holds:

$$A_p(c)=2^{1-p}(1+\lambda c)^p
 =\sum_{q=0}^{p} a_q c^q,\qquad
 a_q=2^{1-p}\binom pq\lambda^q.$$

With $M_{\boldsymbol\alpha}^{s}=\sum_{j\in s}q(r_j)
\mathbf u_j^{\boldsymbol\alpha}$, the degree contraction is

$$S_q^{st}=\sum_{|\boldsymbol\alpha|=q}
 \frac{q!}{\alpha_x!\alpha_y!\alpha_z!}
 M_{\boldsymbol\alpha}^{s}M_{\boldsymbol\alpha}^{t}.$$

For equal species, replace this with
$\tfrac12(S_q^{ss}-\sum_{j\in s}q(r_j)^2)$ to remove self pairs and double
counting. Sum $a_q S_q^{st}$ for the descriptor. The existing NN-adjoint
contraction and differentiated three-variable Horner evaluation give forces
without a second neighbor-pair loop. Integer eligibility uses **exact equality**;
a near-integer parameter is not rounded into this polynomial path.

The packed coefficient region now has room for $a_0,\ldots,a_{16}$. The
component cutoff, previously immediately after $a_{10}$, moves after $a_{16}$.
A named field constant is used by both packing and runtime selection. Both
local degree-accumulator arrays use the same maximum-order constant, avoiding
an order-dependent out-of-bounds write. Packed data are internal, so this does
not change native model files or the public C ABI.

### 21.3 Shared execution and conservative automatic selection

CPU batch evaluation, the LAMMPS CPU batch adapter, and OpenMP target GPU
execution compile the same moment construction, coefficient contraction and
force routines. CPU builds remove the OpenMP directives. The retained atomic
CPU evaluator is still an independent reference and keeps its old order-10
moment limit; it evaluates higher orders directly. Its existing diagnostic now
explicitly identifies that legacy limit, rather than implying that the shared
batch path cannot evaluate higher-order moments.

* Auto (0) retains exact integer orders 1--10 and its existing 16-neighbor rule.
* Direct (1) always keeps the direct evaluator.
* Thresholded moment (2) supports orders 1--16 above the existing neighbor rule.
* Forced moment (3), including LAMMPS `g5 moment`, supports orders 1--16 without
  the neighbor threshold.
* Fractional, near-integer, and orders above 16 use direct evaluation.

High-order eligibility is excluded during auto-mode packing, so an order-16
function cannot enlarge the automatic moment basis for unrelated low-order
functions. This preserves the existing selection policy rather than assuming
that every exact moment transform is faster.

The number of Cartesian monomials through degree $p$ is

$$K(p)=\binom{p+3}{3},\qquad K(10)=286,\quad K(16)=969.$$

Moment arithmetic and storage therefore increase by about **3.39 times**
between those limits. Reuse across angular powers and species pairs matters:
a single high-order channel can lose to direct pairs, while a group sharing
moments across many powers can amortize this work. GPU transfer and launch
costs also matter for small center counts. The benchmark reports both isolated
high-order functions and a shared degree-1--16 family, and does not assert an
unconditional moment speedup.


### 21.4 Measure actual CPU parallel execution

On a combined `target teams distribute parallel do`, an unqualified `if(...)`
applies to both the target and parallel constituents. Thus the original
`if(device /= omp_get_initial_device())` also serialized the parallel region
on the host. Use `if(target:...)` to select the target alone, leaving host
parallelism controlled by the OpenMP runtime. This follows the
[OpenMP specification for the if clause](https://www.openmp.org/spec-html/5.2/openmpse17.html).
The data mapping conditions and scalar numerical routines are unchanged.

The thread benchmark compares prepared common kernels compiled with OpenMP
OFF against ON at 1, 2, 4 and 8 physical cores on one socket. It verifies actual
worker affinities, fixes dynamic teams off, reverses process order, and checks
energies, every force component and virial against the independent retained
reference. Report both $T_{\mathrm{OFF},1}/T_{\mathrm{ON},p}$ and the ON/1 cost;
using ON/1 alone as the baseline would hide directive/runtime overhead.

This is an explicit host-target path. The normal CPU batch and LAMMPS CPU
adapter still use the serial compilation. GPU target execution keeps the same
math, while the retained atomic CPU evaluator remains the compatibility and
independent reference implementation.


## 22. Remove fine-grained force atomics (revision 1.11)

Baseline: `gpu` **0e79af3**, AccelNet **1.0.1**. Revision 1.10 exposed the cost
of enabling OpenMP on the host, including expensive locked updates even with
one worker. This revision changes ownership and reduction of already-computed
forces; descriptor and neural-network formulas stay shared and unchanged.
[Validation and paired before/after timings](docs/validation/atomic-reduction-2026-09-27/README.md)
include OpenMP OFF, ON/1/2/4/8 and H100.

### 22.1 Give a center exclusive ownership of its CSR edge forces

For center row $i$, define its edge interval

$$E_i=\{\mathrm{offsets}(i),\ldots,\mathrm{offsets}(i+1)-1\}.$$

These intervals are disjoint even when two rows refer to the same physical
center or periodic images refer to the same target atom. The host G5 direct
method evaluates an unordered pair $(j,k)$ once and writes both derivatives
into `edge_force(:,j)` and `edge_force(:,k)`. Previously separate edge workers
could update the same $k$, requiring atomic additions. Assigning the whole
$E_i$ to one worker makes both writes exclusive and removes those atomics.
Each worker zeros only its own edge interval before processing it.

Two short scheduling wrappers use one Fortran arithmetic body,
`AccelNetPredictor/src/shared/generic_force_edge_body.inc`. Host pair-once work
items span a center's edge interval; GPU work items span one edge. Both
schedules use `target teams distribute parallel do`, including the host
target context selected by `if(target:...)`. Scalar pair geometry,
angular/radial evaluation and contraction are not duplicated. Non-G5 families
keep the one-edge work mapping.

CMake generates the host center-owned include by removing only OpenMP atomic
directives from that single body. The normal serial module also uses this
stripped include even if global compiler flags enable OpenMP. Generated copies
are build artifacts, not independently maintained numerical implementations.
The GPU edge schedule preserves the original guarded pair-update directives;
`pair_once` is false in this schedule, so they are not executed. Only the final
physical-atom scatter executes explicit atomic additions on GPU.

This distinction was required by an NVHPC 25.3 performance result. Removing
atomic directives even from the GPU's inactive pair-once branch doubled the
G5 direct force kernel's measured duration, about 0.95 to 1.96 ms. Whole
low-order G5 evaluation regressed about 22%. Removing a one-iteration inner
loop, changing the host schedule, and specializing the traversal flag did not
cure it. Nsight Systems showed identical 1024-block/128-thread launches and
156 registers per thread on H100; the final scatter stayed about 0.27 ms.
Restoring the guarded GPU directives while stripping them only for host center
ownership restored GPU performance. This isolates a compiler code-generation
sensitivity to those directives; the exact compiler optimization responsible
has not been established. It is not evidence that executing extra GPU atomics
would be beneficial.

### 22.2 Separate edge ownership from physical-atom scatter

Exclusive edge ownership does not imply exclusive physical-atom ownership.
With edge force $\mathbf f_e$ and edge displacement $\mathbf r_e$, assembly is

$$\mathbf F_a=\sum_{e:\,\mathrm{target}(e)=a}\mathbf f_e
 -\sum_{i:\,\mathrm{center}(i)=a}\sum_{e\in E_i}\mathbf f_e,$$

$$W_{\alpha\beta}=\sum_i\sum_{e\in E_i}
 r_{e,\alpha} f_{e,\beta}.$$

On the host this final assembly is a serial streaming pass, with **no per-edge,
per-center or per-pair atomic force updates**. Descriptor, NN, direct and moment
force computation remain parallel. The serial pass avoids both contention and
an additional reverse-neighbor structure or thread-private force arrays.
Its cost must still be measured as atom and thread counts grow.

On GPU, physical-atom scatter remains parallel and retains atomic updates for
both incoming edge forces and center forces. Removing these safely would need
a reverse adjacency/gather or another conflict-free ownership scheme, including
periodic images and MPI ghost indices. This revision does not substitute a
quadratic scan or introduce unmeasured preprocessing into that path.

The nine virial components use an OpenMP reduction instead of nine explicit
global atomics per center. The runtime may implement its final reduction with
atomics internally; the claim is fewer explicit contended updates, not zero
synchronization inside the OpenMP runtime. Host scatter intentionally uses one
team and one worker; its unqualified `if(parallel_scatter)` deliberately applies
to both the target and parallel constructs. The other compute loops retain
`if(target:...)` so host computations remain parallel.

### 22.3 Keep host worker management consistent and measure clocks

Center ownership alone is insufficient to predict OpenMP performance. A native
host `parallel do` mixed with the other stages' host target teams left 15 live
threads for an eight-thread request on GNU 11.4. The unchanged descriptor stage
slowed as well. Using host target teams for both ownership schedules restored
eight live workers. This scheduling-only fix leaves the serial kernel's compiled
`.text` unchanged. It is a runtime interaction observed with this compiler and
configuration, not a claim about all OpenMP runtimes.

The CPU governor also moved between 3.9 and 2.8 GHz. At matched frequencies the
G5 direct OFF/1 before/after difference is below 1%; mixed-frequency wall times
must not be mistaken for an arithmetic regression. The diagnostic script
`probe_openmp_ownership.py` records live worker counts and frequency samples
without changing machine-wide clock settings.

In the final paired sweep, eight-thread before/after improvements range from
1.03 to 1.48 times, while GPU times change by less than 0.5%. Thread scaling
remains nonmonotonic for cheap workloads. The validation report retains the
full OFF/1, ON/1/2/4/8 and H100 table, raw samples, and limitations. Numerical
checks passed; a separate strict n2p2 type-21 direct parity gate narrowly missed
its unchanged limit (CPU/n2p2 1.101346 versus 1.10). This is recorded as a failure,
not hidden by the successful correctness and other performance tests.

## 23. Unify potential-inference entry points (revision 1.12)

Baseline: `gpu` **543b180**, AccelNet **1.0.1**. Structure/file evaluation,
Fortran/C atomic energy and force/virial calls, ordinary CSR batches, the
ænet-compatible structural-fingerprint API, and LAMMPS inference now reach the
common numerical backend. The ordinary CPU instance still strips OpenMP at
build time. The public atomic ABI and additive force/virial semantics remain.
Former evaluators and exact old API snapshots live in `legacy/cpu-reference/`;
`evaluate_batch_reference` is an explicit reference call, never a fallback.
Standalone low-level descriptor/NN compatibility utilities retain those
reference routines. Model types, loading and shared scalar formulas are not
forked into a second data model.

### 23.1 Compose descriptor blocks before the neural network

For element $s$, let $D_{sc}$ denote component $c$ and $I_{sc}$ its original
input-coordinate range. Preserve the original network $N_s$ and assemble

$$G_{I_{sc}}=D_{sc}(R),\qquad E_i=N_s(S_s(G-b_s))/a_s+e_s.$$

Evaluate the NN once, then scatter its input adjoint back into each block:

$$q_{I_{sc}}=-\frac{\partial E_i}{\partial G_{I_{sc}}},\qquad
 F=\sum_c J_{sc}^{T}q_{I_{sc}}.$$

Each component uses the existing common descriptor and contraction kernels.
GPU descriptor values, adjoints and component forces stay on device between
phases. Component offsets preserve input ordering even when descriptor families
are interleaved or differ by element. Missing components receive zero adjoints.
The same staging supports several Chebyshev components and Chebyshev/LJ/Behler
mixtures. It adds no approximation and no per-component independent NN.
Existing single-component/grouped models keep their normal fused execution.

Energy-only calls skip the NN reverse pass and force stage. SFB calls stop after
descriptor assembly. Descriptor implementations which fuse values and saved
Jacobians may still perform derivative preparation in their value stage; this
revision does not claim every energy-only instruction is eliminated.

### 23.2 Adapt per-atom calls with one CSR row

Represent a central environment by local slots $0,1,\ldots,n$, with one distinct
slot for each periodic image. Evaluate the common kernel once and fold its
forces onto the caller's indices only after forming the image-displacement
virial. Duplicate physical target indices therefore preserve both additive
forces and $W_{ab}=\sum_j r_{j,a}f_{j,b}$. The private atomic API caches packed
model metadata and scratch, invalidated by loading and mode setters. Public
mutable object models are repacked so direct edits cannot leave stale weights.

### 23.3 Workspace lifetime and independent validation

A flat component array avoids a GNU 11 recursive-finalization compiler failure.
Each component owns its own model mapping and scratch. Release component device
mappings explicitly before deallocating/resizing component arrays: relying only
on inherited/array finalization left stale NVHPC 25.3 present-table entries in
the first GPU trial. Scalar/rank finalizers also cover scope exit and C handles.

Tests compare the new entry points against the retained reference, upstream
ænet/n2p2, and coordinate/strain finite differences. Separate before/after
executables use the identical public-API benchmark source, with OpenMP compiled
out, so moving the old structure API does not turn the baseline into another
call to the new kernel. See [the migration report](docs/validation/unified-api-2026-09-27/README.md)
for the measured public-API and GPU timings and their limits.


### 23.4 Tile independent moment accumulators without splitting CPU/GPU source

For monomial index $a=(a_x,a_y,a_z)$ and species weight $s_j$, construct

$$M_a=\sum_j f_j u_{jx}^{a_x}u_{jy}^{a_y}u_{jz}^{a_z},\qquad
  \widetilde M_a=\sum_j s_jf_j u_{jx}^{a_x}u_{jy}^{a_y}u_{jz}^{a_z}.$$

One work item now accumulates four consecutive indices $a$ while traversing
neighbors. Each lane preserves its original neighbor summation order; lanes
are independent. Geometry/cutoff loads are shared, and fixed-size accumulator
arrays expose independent instructions to the CPU compiler. The identical
loop runs on GPU. The last tile uses padded lanes, discarding their outputs.
Eight lanes did not improve the paired CPU or GPU probe and were rejected.

The initial unified per-atom energy-only path took about 8.05 ms per 192-atom
Ti/O structure against 4.24 ms for the former evaluator. Profiling attributed
most time to moment construction, not the network. Four-lane tiling, skipping
unused Chebyshev/activation derivatives, and avoiding redundant packed-model
comparisons in the explicitly invalidated private C cache reduced this to
about 5.19 ms. This diagnostic still leaves a roughly 22% energy-only regression;
force-inclusive and bulk API timings must be reported separately. Final paired
measurements and raw samples are in the migration report. This is not a claim
that all migrated entry points are faster.


## 24. Restore energy-only CPU performance in the common kernels (revision 1.13)

### 24.1 One accumulation per species, then form weighted channels

The four-moment tile in Section 23.4 repeats a neighbor traversal for each tile
and loads cached powers with a stride equal to the edge capacity. Both effects
matter for a one-center CPU call. A diagnostic split on the real Ti/O model
also found substantial geometry/radial time: about 2.50 ms there versus 2.08 ms
in moment construction/contraction per 192 individual calls. The optimization
therefore addresses both stages rather than assuming all time is in moments.

For species $t$ and monomial $a=(a_x,a_y,a_z)$, accumulate

$$M_a^{(t)}=\sum_{j:t_j=t} f_j
u_{jx}^{a_x}u_{jy}^{a_y}u_{jz}^{a_z}.$$

Here $u_j=r_j/|r_j|$, and the angular cutoff factor is $f_j$. For a central
species $s$, the two required channels are exactly

$$M_a^{(0)}=\sum_t M_a^{(t)},\qquad
M_a^{(w)}=\sum_t w_{ts}M_a^{(t)}.$$

These identities hold for zero, negative and fractional weights. Pair moments
are still formed with the same diagonal subtraction,

$$P_a^{(0)}=\tfrac12\left[(M_a^{(0)})^2-
\sum_j f_j^2u_j^{2a}\right],\qquad
P_a^{(w)}=\tfrac12\left[(M_a^{(w)})^2-
\sum_j w_{t_js}^2f_j^2u_j^{2a}\right].$$

The existing polynomial coefficients and multinomial factors contract these
moments into Chebyshev values. The existing differentiated contraction supplies
forces; only the moment construction has changed. Accumulation order changes
between species, so the guarantee is FP64 tolerance agreement, not bitwise
identity. No angular series is truncated.

### 24.2 Contiguous powers and independent center/species work

For each center, a reverse CSR traversal prepends active neighbors to a list
for their species. Each resulting list retains its forward CSR order. Building
these lists costs $O(N)$; traversing all species lists also visits each active
neighbor once. A work item owns one center/species pair and its contiguous
moment block. There are no shared writes or atomic additions.

Generate the three one-dimensional powers by recurrence. At fixed $a_x,a_y$,
hoist $c=f_j u_{jx}^{a_x}u_{jy}^{a_y}$ and update the consecutive $a_z$ entries:

$$M_{(a_x,a_y,a_z)}^{(t)}\mathrel{+}=c\,u_{jz}^{a_z}.$$

The contiguous inner update exposes vector operations to the serial CPU
compiler. Per-edge global Chebyshev power construction is no longer needed;
G5 retains its own power cache. Device scratch is allocated/mapped with the
workspace, and the device helper receives explicit dimensions. An early
assumed-shape helper produced an NVHPC 25.3 invalid device read; memcheck
identified a host-like address, and the explicit-dimension version removed
that failure.

CPU and GPU execute the same species-list and power/moment loops. Only the
OpenMP team-size limit differs (32 on host, 4 on device). Serial compilation
removes the directives entirely. An earlier whole-center work item recovered
CPU speed but lost GPU parallelism; splitting by species reduced that GPU
penalty in paired probes. This is a scheduling choice, not a separate formula.

### 24.3 Remove work that energy-only callers do not consume

Energy-only evaluation does not clear, download or add force/virial arrays.
Sentinel tests call energy-only immediately after a force calculation and then
perform further force calls, checking both untouched caller accumulators and
correct scratch reuse. The cutoff value/derivative helpers include the existing
shared scalar bodies in the kernel translation unit to allow inlining.

Workspace capacities are cached only after validation against the cached
model. Repeated private atomic calls skip shape scans when that model and the
requested dimensions still fit. Model changes bypass the fast path, and buffer
release clears its capacities. Public mutable model handling and cache
invalidation rules are unchanged.

The leading row capacity is padded to $N\equiv1\pmod 8$ for all descriptor
families. A fixed-row traversal of many columns at $N=512$ otherwise has a
4096-byte stride, which can repeatedly map to the same CPU cache sets. This
padding does not change active rows or arithmetic, and is used on GPU too.
A paired 512-atom probe reduced structure energy from 162.06 to 151.11 ms
(legacy 146.70 ms), and force-inclusive batch from 45.87 to 33.43 ms. This is
evidence for a cache-layout effect, not a hardware-counter attribution.

### 24.4 Exactness boundary for n2p2 types 13, 21 and 24

All three descriptors already have exact direct values, analytic forces and
virial in the common CPU/GPU implementation. They are not unsupported model
types. The absence is a general finite single-neighbor moment factorization:
type 13 has a third-distance Gaussian and cutoff; types 21/24 have a
third-distance compact window and an angular window. General parameters do not
reduce to the finite separable polynomial used for integer-power G5. A finite
series approximation would violate the requested exactness. Special parameter
cases need a separate algebraic justification before adding a moment dispatch.
See the [upstream definitions](https://compphysvienna.github.io/n2p2/api/symmetry_function_types.html)
and Sections 19--21 for the existing direct implementations and exact G5 scope.

Final numerical checks, paired timings, version identities and limitations are
recorded in the [energy-only validation report](docs/validation/energy-common-2026-09-27/README.md).
