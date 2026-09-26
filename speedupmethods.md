# AccelNet CPU/GPU speedup methods

**Document version 1.5 — 2026-09-27 (JST).**

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
| Optimization source revision | `gpu` checkpoint **`2bcc603`** (revision 1.4), followed by the G4 scheduling and scalar-cache changes in Section 16; source hashes and the implementation diff identify the measured revision 1.5 |
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
| LJ and Behler G1–G5 | GPU support implemented; G4/G5 use direct pairs; forced GPU G5 moments are rejected |
| Generic radial caching and LJ component fusion | Implemented in the common serial CPU/GPU source; measurements in Section 12 |
| Behler angular coefficient contraction and differentiated Horner | Retained for G5; the former G4 path in Sections 13–14 is superseded by Section 15 |
| G4 value/Jacobian evaluation | Shared value/derivative loop with scalar caches and disjoint descriptor owners; CPU uses one owner, GPU uses a flat center/owner launch; force contraction reads the saved Jacobian; Sections 15–16 |
| LAMMPS GPU package integration | CUDA adapter plus Fortran OpenMP target; AMD/HIP/OpenCL interoperability is not implemented |

The default CPU shared path applies to the **CSR batch API when every element has
one Chebyshev component**. `evaluate_batch_reference` retains the independent old
CPU path. Object APIs, atomic Fortran/C APIs, CLI, and ordinary LAMMPS
`pair_style accelnet` still use their established CPU implementation. CPU LJ/Behler
and composite Chebyshev batches also retain that implementation. Multiple or mixed
Chebyshev components within one element are not supported by the GPU backend.

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
