//! Banded storage for assembled stiffness matrices, with the graph helpers
//! that decide how wide the band has to be.
//!
//! The lower triangle is stored row-wise: `A(i, j)` with `i >= j` sits at
//! `data[i * bw + (i - j)]`, so walking along a row is a contiguous scan. The
//! upper triangle is not stored; `get` reads it from the same slot as its
//! lower twin.

use super::error::SolveError;
use super::index::{NewIndex, OldIndex, Permutation};
use super::tolerance::{PIVOT_RELATIVE, RESIDUAL};
use super::Real;
use std::collections::VecDeque;

/// Symmetric matrix stored by lower-triangle rows in a band of half-width `bw`.
#[derive(Clone, Debug)]
pub struct BandedMatrix {
    n: usize,
    bw: usize,
    data: Vec<Real>,
}

impl BandedMatrix {
    /// Zero matrix of dimension `n` and half-bandwidth `bw` (at least 1).
    ///
    /// # Panics
    ///
    /// Panics when `n * bw` overflows. A bandwidth that large is almost always
    /// a numbering mistake; run [`renumber`] on the element connectivities
    /// before allocating.
    pub fn new(n: usize, bw: usize) -> Self {
        let bw = bw.max(1);
        let slots = n
            .checked_mul(bw)
            .expect("banded matrix size overflow; check the equation numbering");
        Self {
            n,
            bw,
            data: vec![0.0; slots],
        }
    }

    /// Identity matrix of dimension `n`, at the narrowest band.
    pub fn identity(n: usize) -> Self {
        let mut m = Self::new(n, 1);
        for i in 0..n {
            m.add_symmetric(i, i, 1.0);
        }
        m
    }

    /// Matrix dimension.
    #[inline]
    pub fn dim(&self) -> usize {
        self.n
    }

    /// Half-bandwidth: entries with `|i - j| < bandwidth` are stored.
    #[inline]
    pub fn bandwidth(&self) -> usize {
        self.bw
    }

    /// Number of coefficient slots, `n * bandwidth`.
    #[inline]
    pub fn stored(&self) -> usize {
        self.data.len()
    }

    /// Row-major lower-triangle slot for `i >= j`.
    #[inline]
    fn idx(i: usize, j: usize, bw: usize) -> usize {
        i * bw + (i - j)
    }

    /// First column in row `i`'s stored band, inclusive.
    #[inline]
    fn band_start(i: usize, bw: usize) -> usize {
        (i + 1).saturating_sub(bw)
    }

    /// One past the last column in row `i`'s stored band, exclusive.
    #[inline]
    fn band_end(i: usize, bw: usize, n: usize) -> usize {
        (i + bw).min(n)
    }

    /// Adds `v` to both `(i, j)` and `(j, i)`. Call this **once per symmetric
    /// pair**: an assembly loop over all ordered pairs doubles every
    /// off-diagonal, which no solve will detect.
    ///
    /// # Panics
    /// Panics when `i` or `j` is outside the matrix dimensions, when the entry
    /// lies outside the stored half-band, or when `v` is not finite. The band
    /// check runs in every build profile: a dropped stiffness contribution
    /// would yield plausible but wrong displacements, so an out-of-band write
    /// is refused rather than silently ignored.
    #[inline]
    pub fn add_symmetric(&mut self, i: usize, j: usize, v: Real) {
        assert!(i < self.n && j < self.n, "banded index out of range");
        assert!(v.is_finite(), "stiffness contribution must be finite");
        let (hi, lo) = if i >= j { (i, j) } else { (j, i) };
        assert!(
            hi - lo < self.bw,
            "entry ({i},{j}) outside bandwidth {}",
            self.bw
        );
        self.data[Self::idx(hi, lo, self.bw)] += v;
    }

    /// [`BandedMatrix::add_symmetric`] with indices already in the solver's
    /// numbering.
    #[inline]
    pub fn add_symmetric_at(&mut self, i: NewIndex, j: NewIndex, v: Real) {
        self.add_symmetric(i.get(), j.get(), v);
    }

    /// Reads `A(i, j)`, returning zero outside the band.
    #[inline]
    pub fn get(&self, i: usize, j: usize) -> Real {
        assert!(i < self.n && j < self.n, "banded index out of range");
        let (hi, lo) = if i >= j { (i, j) } else { (j, i) };
        if hi - lo >= self.bw {
            0.0
        } else {
            self.data[Self::idx(hi, lo, self.bw)]
        }
    }

    /// [`BandedMatrix::get`] with indices already in the solver's numbering.
    #[inline]
    pub fn get_at(&self, i: NewIndex, j: NewIndex) -> Real {
        self.get(i.get(), j.get())
    }

    /// Matrix-vector product `self * rhs`.
    ///
    /// # Panics
    /// Panics when `rhs` has a different length from the matrix.
    pub fn mul_vec(&self, rhs: &[Real]) -> Vec<Real> {
        assert_eq!(rhs.len(), self.n, "rhs length must match matrix size");
        let mut out = vec![0.0; self.n];
        for i in 0..self.n {
            let mut sum = self.data[i * self.bw] * rhs[i];
            let first = Self::band_start(i, self.bw);
            for j in first..i {
                sum += self.data[i * self.bw + (i - j)] * rhs[j];
            }
            // Upper entries live in the later row that owns the symmetric pair.
            let last = Self::band_end(i, self.bw, self.n);
            for j in i + 1..last {
                sum += self.data[j * self.bw + (j - i)] * rhs[j];
            }
            out[i] = sum;
        }
        out
    }

    /// Factors the matrix as `A = L·D·Lᵀ` and returns the owning [`Factored`].
    ///
    /// The matrix is consumed: on success the result holds both the factor and
    /// an unfactored copy for the residual check, so a stale factor cannot be
    /// reused after the stiffness changes. On error everything is dropped.
    ///
    /// The factorisation is unpivoted, which needs a positive definite matrix.
    /// The pivot test is relative to each equation's original diagonal, so the
    /// verdict does not depend on the units of the stiffness matrix. A pivot
    /// that is negative, non-finite or too small is refused rather than
    /// propagated, so a caller never receives NaNs.
    pub fn factor(mut self) -> Result<Factored, SolveError> {
        let n = self.n;
        let bw = self.bw;
        let a = self.clone();
        let original_diagonal: Vec<Real> = (0..n).map(|i| self.data[i * bw]).collect();
        let mut d = vec![0.0; n];

        for i in 0..n {
            let k0 = Self::band_start(i, bw);
            let row = i * bw;

            // Off-diagonals first: row i only needs D and the rows already
            // finished, so the diagonal can come last.
            for j in k0..i {
                let mut t = self.data[row + (i - j)];
                // k < j and i - k < bw; the second implies the first.
                for (k, dk) in d[..j].iter().enumerate().skip(k0) {
                    t -= self.data[row + (i - k)] * self.data[j * bw + (j - k)] * dk;
                }
                self.data[row + (i - j)] = t / d[j];
            }

            let mut s = self.data[row];
            for (k, dk) in d[..i].iter().enumerate().skip(k0) {
                let l = self.data[row + (i - k)];
                s -= l * l * dk;
            }
            if !s.is_finite() {
                return Err(SolveError::NotFinite { equation: i });
            }
            if s < 0.0 {
                return Err(SolveError::NotPositiveDefinite { equation: i });
            }
            if s <= PIVOT_RELATIVE * original_diagonal[i] {
                return Err(SolveError::Singular { equation: i });
            }
            d[i] = s;
        }
        Ok(Factored { a, l: self, d })
    }

    /// Infinity norm of `A·x - b`, relative to `‖b‖∞` when `b` is non-zero and
    /// absolute otherwise. Compare with [`RESIDUAL`](super::tolerance::RESIDUAL).
    pub fn residual(&self, x: &[Real], b: &[Real]) -> Real {
        self.residual_report(x, b).1
    }

    /// The relative residual and the equation that produced it.
    fn residual_report(&self, x: &[Real], b: &[Real]) -> (usize, Real) {
        assert_eq!(x.len(), self.n, "solution length must match matrix size");
        assert_eq!(b.len(), self.n, "rhs length must match matrix size");
        let y = self.mul_vec(x);
        let (equation, worst) = y
            .iter()
            .zip(b)
            .enumerate()
            .map(|(i, (&yi, &bi))| (i, (yi - bi).abs()))
            .fold(
                (0, 0.0),
                |best, current| {
                    if current.1 > best.1 {
                        current
                    } else {
                        best
                    }
                },
            );
        let scale = b
            .iter()
            .fold(0.0 as Real, |largest, &value| largest.max(value.abs()));
        let value = if scale > 0.0 { worst / scale } else { worst };
        (equation, value)
    }

    /// Expands the band to a dense row-major matrix. Test-only: the kernel
    /// never needs the dense form, and expanding it would defeat the point of
    /// banded storage.
    #[cfg(test)]
    fn to_dense(&self) -> Vec<Vec<Real>> {
        let mut dense = vec![vec![0.0; self.n]; self.n];
        for i in 0..self.n {
            let first = Self::band_start(i, self.bw);
            for j in first..=i {
                dense[i][j] = self.data[Self::idx(i, j, self.bw)];
                dense[j][i] = dense[i][j];
            }
        }
        dense
    }
}

/// A factorised [`BandedMatrix`]. Only this type can solve, which makes
/// "solve before factor", "factor twice" and "add a stiffness then reuse the
/// old factor" unrepresentable.
#[derive(Clone, Debug)]
pub struct Factored {
    a: BandedMatrix,
    l: BandedMatrix,
    d: Vec<Real>,
}

impl Factored {
    /// The diagonal `D` of `L·D·Lᵀ`, all positive for a factored SPD matrix.
    pub fn diagonal(&self) -> &[Real] {
        &self.d
    }

    /// Solves `A·x = b`, overwriting `b` with the solution.
    ///
    /// The right-hand side is checked for finite values before substitution,
    /// and the relative residual is checked afterwards. Both failures are
    /// typed [`SolveError`]s, so a bad load cannot put NaN into a result.
    ///
    /// # Panics
    /// Panics when `b` has a different length from the matrix.
    pub fn solve(&self, b: &mut [Real]) -> Result<(), SolveError> {
        assert_eq!(b.len(), self.l.n, "rhs length must match matrix size");
        for (equation, &value) in b.iter().enumerate() {
            if !value.is_finite() {
                return Err(SolveError::NonFiniteRightHandSide { equation });
            }
        }
        let rhs = b.to_vec();
        let n = self.l.n;
        let bw = self.l.bw;

        // Forward: L·y = b.
        for i in 0..n {
            let k0 = BandedMatrix::band_start(i, bw);
            let row = i * bw;
            let mut s = b[i];
            for (k, bk) in b[..i].iter().enumerate().skip(k0) {
                s -= self.l.data[row + (i - k)] * bk;
            }
            b[i] = s;
        }
        // Diagonal: z = D⁻¹·y.
        for i in 0..n {
            b[i] /= self.d[i];
        }
        // Backward: Lᵀ·x = z.
        for i in (0..n).rev() {
            let jmax = BandedMatrix::band_end(i, bw, n);
            let mut s = b[i];
            for (j, bj) in b.iter().enumerate().take(jmax).skip(i + 1) {
                s -= self.l.data[j * bw + (j - i)] * bj;
            }
            b[i] = s;
        }

        let (equation, residual) = self.a.residual_report(b, &rhs);
        if !residual.is_finite() || residual > RESIDUAL {
            return Err(SolveError::ResidualTooLarge { equation, residual });
        }
        Ok(())
    }
}

/// Reverse Cuthill-McKee ordering, as `p[new_index] = old_index`.
///
/// Disconnected components are ordered one after another. This is the
/// low-level routine; [`renumber`] wraps its result in a typed [`Permutation`].
fn rcm(n: usize, adj: &[Vec<usize>]) -> Vec<usize> {
    assert_eq!(adj.len(), n, "adjacency must have one entry per equation");
    for (u, neighbours) in adj.iter().enumerate() {
        for &v in neighbours {
            assert!(
                v < n,
                "adjacency of equation {u} references out-of-range equation {v}"
            );
        }
    }

    let mut visited = vec![false; n];
    let mut order: Vec<usize> = Vec::with_capacity(n);

    // Start with the lowest-degree node; disconnected components are handled
    // by continuing with the next unvisited node in degree order.
    let mut seeds: Vec<usize> = (0..n).collect();
    seeds.sort_by_key(|&i| adj[i].len());

    for &seed in &seeds {
        if visited[seed] {
            continue;
        }
        let mut queue = VecDeque::new();
        queue.push_back(seed);
        visited[seed] = true;

        while let Some(u) = queue.pop_front() {
            order.push(u);
            let mut neighbours: Vec<usize> =
                adj[u].iter().copied().filter(|&v| !visited[v]).collect();
            neighbours.sort_by_key(|&v| adj[v].len());
            for v in neighbours {
                if !visited[v] {
                    visited[v] = true;
                    queue.push_back(v);
                }
            }
        }
    }

    order.reverse();
    order
}

/// Build a bandwidth-reducing numbering from element connectivities expressed
/// in the model's original numbering.
///
/// Elements couple every pair of their DOFs, so each element becomes a clique
/// in the adjacency graph. Returns the permutation, the connectivities in the
/// solver's numbering, and the half-bandwidth those connectivities require.
///
/// # Panics
/// Panics when an element references an equation outside `0..n`.
pub fn renumber(n: usize, elements: &[Vec<OldIndex>]) -> (Permutation, Vec<Vec<NewIndex>>, usize) {
    let mut adj = vec![Vec::new(); n];
    for element in elements {
        for &dof in element {
            assert!(
                dof.get() < n,
                "element references out-of-range equation {}",
                dof.get()
            );
        }
        for (k, &a) in element.iter().enumerate() {
            for &b in element.iter().skip(k + 1) {
                if a != b {
                    adj[a.get()].push(b.get());
                    adj[b.get()].push(a.get());
                }
            }
        }
    }
    let old_order: Vec<OldIndex> = rcm(n, &adj).into_iter().map(OldIndex::new).collect();
    let permutation = Permutation::from_old_order(&old_order);
    let reordered: Vec<Vec<NewIndex>> = elements
        .iter()
        .map(|element| permutation.to_new_all(element))
        .collect();
    let raw: Vec<Vec<usize>> = reordered
        .iter()
        .map(|element| element.iter().map(|&dof| dof.get()).collect())
        .collect();
    let bw = bandwidth(&raw);
    (permutation, reordered, bw)
}

/// Half-bandwidth implied by a set of element-to-equation connectivities.
///
/// The result is the smallest `bw` such that every element spans fewer than
/// `bw` equations apart, matching [`BandedMatrix::new`]'s `|i - j| < bw`
/// contract. At least 1 is returned.
pub fn bandwidth(elements: &[Vec<usize>]) -> usize {
    let mut bw = 1usize;
    for el in elements {
        let (mut lo, mut hi) = (usize::MAX, 0usize);
        for &d in el.iter() {
            lo = lo.min(d);
            hi = hi.max(d);
        }
        if hi >= lo {
            bw = bw.max(hi - lo + 1);
        }
    }
    bw
}

#[cfg(test)]
mod tests {
    use super::super::error::SolveError;
    use super::super::index::OldIndex;
    use super::super::tolerance::RESIDUAL;
    use super::{bandwidth, rcm, renumber, BandedMatrix, Real};

    fn next(seed: &mut u64) -> Real {
        *seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        // Generated in f64 and cast, so the oracle's precision does not follow
        // the `Real` alias.
        (((*seed >> 11) as f64) / ((1u64 << 53) as f64)) as Real
    }

    fn seeded_spd(n: usize, bw: usize) -> BandedMatrix {
        let mut seed = 0x9e37_79b9_7f4a_7c15;
        let mut m = BandedMatrix::new(n, bw);
        for i in 0..n {
            m.add_symmetric(i, i, 2.0 + next(&mut seed));
            for j in i + 1..BandedMatrix::band_end(i, bw, n) {
                m.add_symmetric(i, j, 0.2 * next(&mut seed));
            }
        }
        m
    }

    fn dense_mul(dense: &[Vec<Real>], rhs: &[Real]) -> Vec<Real> {
        dense
            .iter()
            .map(|row| row.iter().zip(rhs).map(|(a, b)| a * b).sum())
            .collect()
    }

    /// Gaussian elimination with partial pivoting, on the dense oracle only.
    /// The banded solver under test never uses it.
    fn dense_solve(dense: &[Vec<Real>], b: &[Real]) -> Vec<Real> {
        let n = b.len();
        let mut a: Vec<Vec<Real>> = dense.to_vec();
        let mut x = b.to_vec();
        for col in 0..n {
            let pivot = (col..n)
                .max_by(|&i, &j| a[i][col].abs().partial_cmp(&a[j][col].abs()).unwrap())
                .unwrap();
            a.swap(col, pivot);
            x.swap(col, pivot);
            for row in col + 1..n {
                let f = a[row][col] / a[col][col];
                for k in col..n {
                    a[row][k] -= f * a[col][k];
                }
                x[row] -= f * x[col];
            }
        }
        for row in (0..n).rev() {
            let s: Real = (row + 1..n).map(|k| a[row][k] * x[k]).sum();
            x[row] = (x[row] - s) / a[row][row];
        }
        x
    }

    /// Factor an SPD matrix and solve `A·x = b`, returning `x`.
    fn solve_spd(m: &BandedMatrix, b: &[Real]) -> Vec<Real> {
        let factored = m.clone().factor().expect("SPD system must factor");
        assert!(
            factored.diagonal().iter().all(|&di| di > 0.0),
            "SPD factor is positive"
        );
        let mut x = b.to_vec();
        factored.solve(&mut x).expect("solve must converge");
        x
    }

    #[test]
    fn add_reports_both_halves_of_a_symmetric_entry() {
        let mut m = BandedMatrix::new(6, 2);
        m.add_symmetric(2, 1, 3.0);
        assert_eq!(m.get(2, 1), 3.0);
        assert_eq!(m.get(1, 2), 3.0);
        m.add_symmetric(1, 1, 2.0);
        m.add_symmetric(1, 1, 3.0);
        assert_eq!(m.get(1, 1), 5.0);
        assert_eq!(m.get(0, 5), 0.0);
    }

    #[test]
    fn dense_expansion_keeps_the_band_and_symmetry() {
        let mut m = BandedMatrix::new(4, 2);
        m.add_symmetric(0, 0, 1.0);
        m.add_symmetric(1, 0, 2.0);
        m.add_symmetric(1, 1, 3.0);
        m.add_symmetric(3, 2, 4.0);
        assert_eq!(
            m.to_dense(),
            vec![
                vec![1.0, 2.0, 0.0, 0.0],
                vec![2.0, 3.0, 0.0, 0.0],
                vec![0.0, 0.0, 0.0, 4.0],
                vec![0.0, 0.0, 4.0, 0.0],
            ]
        );
    }

    #[test]
    fn identity_and_stored_count_match_the_band() {
        let m = BandedMatrix::identity(4);
        assert_eq!(m.dim(), 4);
        assert_eq!(m.bandwidth(), 1);
        assert_eq!(m.stored(), 4);
        let v = vec![1.0, 2.0, 3.0, 4.0];
        assert_eq!(m.mul_vec(&v), v);
        assert_eq!(BandedMatrix::new(3, 0).bandwidth(), 1);
        assert_eq!(BandedMatrix::new(3, 0).stored(), 3);
    }

    #[test]
    fn mul_vec_matches_dense_on_seeded_spd_matrices() {
        // Strict diagonal dominance with a positive diagonal makes these SPD.
        for bw in 1..=3 {
            let m = seeded_spd(12, bw);
            let dense = m.to_dense();
            let mut seed = 42;
            for _ in 0..4 {
                let v: Vec<Real> = (0..12).map(|_| 2.0 * next(&mut seed) - 1.0).collect();
                let expected = dense_mul(&dense, &v);
                for (got, want) in m.mul_vec(&v).iter().zip(&expected) {
                    assert!(
                        (got - want).abs() < 1.0e-12,
                        "bandwidth {bw}: got {got}, dense oracle {want}"
                    );
                }
            }
        }
    }

    // Refuses the write in every build profile, debug and release alike: a
    // dropped stiffness contribution is the one failure mode that must not be
    // silent. This replaces the old debug-only version that let release
    // builds ignore the write.
    #[test]
    #[should_panic(expected = "outside bandwidth")]
    fn out_of_band_write_is_refused_in_every_profile() {
        let mut m = BandedMatrix::new(5, 2);
        m.add_symmetric(0, 4, 1.0);
    }

    #[test]
    #[should_panic(expected = "stiffness contribution must be finite")]
    fn non_finite_stiffness_is_refused_at_assembly() {
        let mut m = BandedMatrix::new(2, 1);
        m.add_symmetric(0, 0, Real::NAN);
    }

    #[test]
    #[should_panic(expected = "size overflow")]
    fn new_rejects_a_size_overflow() {
        BandedMatrix::new(usize::MAX, 2);
    }

    #[test]
    #[should_panic(expected = "rhs length")]
    fn mul_vec_rejects_wrong_rhs_length() {
        BandedMatrix::new(2, 1).mul_vec(&[1.0]);
    }

    #[test]
    fn bandwidth_covers_the_widest_element() {
        assert_eq!(bandwidth(&[]), 1);
        assert_eq!(bandwidth(&[vec![0, 3], vec![1, 2]]), 4);
        // A single element spanning equations 1..=6 needs half-bandwidth 6.
        assert_eq!(bandwidth(&[vec![1, 2, 6]]), 6);
    }

    #[test]
    fn rcm_orders_a_chain_and_the_band_stays_narrow() {
        // Chain of 20 nodes, each connected to its neighbours.
        let n = 20;
        let mut adj = vec![Vec::new(); n];
        for i in 0..n - 1 {
            adj[i].push(i + 1);
            adj[i + 1].push(i);
        }
        let perm = rcm(n, &adj);
        assert_eq!(perm.len(), n);

        // Elements connect consecutive DOFs only -> bandwidth 2 in the
        // reordered numbering as well.
        let els: Vec<Vec<usize>> = (0..n - 1).map(|i| vec![perm[i], perm[i + 1]]).collect();
        assert!(bandwidth(&els) <= 2);

        // The permutation visits every node exactly once.
        let mut sorted = perm.clone();
        sorted.sort_unstable();
        assert_eq!(sorted, (0..n).collect::<Vec<_>>());
    }

    #[test]
    #[should_panic(expected = "one entry per equation")]
    fn rcm_rejects_short_adjacency() {
        rcm(3, &[vec![1], vec![0]]);
    }

    #[test]
    #[should_panic(expected = "out-of-range equation")]
    fn rcm_rejects_an_out_of_range_neighbour() {
        rcm(2, &[vec![1], vec![9]]);
    }

    #[test]
    fn renumber_returns_a_typed_permutation_and_a_narrow_band() {
        let elements: Vec<Vec<OldIndex>> = (0..19)
            .map(|i| vec![OldIndex::new(i), OldIndex::new(i + 1)])
            .collect();
        let (permutation, reordered, bw) = renumber(20, &elements);
        assert_eq!(permutation.len(), 20);
        assert_eq!(reordered.len(), elements.len());
        assert!(bw <= 2, "unexpected bandwidth {bw}");
        for (original, renamed) in elements.iter().zip(&reordered) {
            assert_eq!(
                renamed,
                &permutation.to_new_all(original),
                "renumbered element must use typed indices"
            );
        }
    }

    #[test]
    fn connectivity_decides_the_band_and_admits_every_element_entry() {
        // Two elements: one spanning equations 0..=3, one spanning 3..=5.
        let elements = vec![vec![0, 1, 2, 3], vec![3, 4, 5]];
        let n = 6;

        // The band comes from the connectivity, not from a caller-chosen
        // number: the widest element spans four equations.
        let bw = bandwidth(&elements);
        assert_eq!(bw, 4);

        let mut m = BandedMatrix::new(n, bw);
        for el in &elements {
            for i in 0..el.len() {
                for j in i..el.len() {
                    // One call per symmetric pair: adding (j, i) as well would
                    // double every off-diagonal.
                    m.add_symmetric(el[i], el[j], 1.0);
                }
            }
        }

        assert_eq!(m.dim(), n);
        assert_eq!(m.bandwidth(), bw);
        assert_eq!(m.get(0, 0), 1.0);
        assert_eq!(m.get(0, 3), 1.0);
        assert_eq!(m.get(3, 5), 1.0);
        // The mirror entry is reported as well.
        assert_eq!(m.get(5, 3), 1.0);
    }

    #[test]
    fn factor_and_solve_reproduce_the_closed_form_3x3_system() {
        // A = [[4, 1, 0], [1, 3, 1], [0, 1, 2]] with x = [1, 2, -1] gives
        // b = [6, 6, 0], solved in closed form by the tridiagonal algorithm.
        let mut m = BandedMatrix::new(3, 2);
        m.add_symmetric(0, 0, 4.0);
        m.add_symmetric(1, 0, 1.0);
        m.add_symmetric(1, 1, 3.0);
        m.add_symmetric(2, 1, 1.0);
        m.add_symmetric(2, 2, 2.0);

        let x = solve_spd(&m, &[6.0, 6.0, 0.0]);
        let expected = [1.0, 2.0, -1.0];
        for (got, want) in x.iter().zip(&expected) {
            assert!((got - want).abs() < 1.0e-12, "got {got}, want {want}");
        }
        // The residual is measured on the unfactored matrix.
        assert!(m.residual(&x, &[6.0, 6.0, 0.0]) < RESIDUAL);
    }

    #[test]
    fn solve_agrees_with_the_dense_oracle_on_seeded_spd_systems() {
        for bw in 1..=3 {
            let m = seeded_spd(20, bw);
            let dense = m.to_dense();
            let mut seed = 7;
            let x_true: Vec<Real> = (0..20).map(|_| 2.0 * next(&mut seed) - 1.0).collect();
            let b = dense_mul(&dense, &x_true);

            let x = solve_spd(&m, &b);

            let oracle = dense_solve(&dense, &b);
            for (got, want) in x.iter().zip(&oracle) {
                assert!(
                    (got - want).abs() < 1.0e-10,
                    "bandwidth {bw}: got {got}, oracle {want}"
                );
            }
            for (got, want) in x.iter().zip(&x_true) {
                assert!(
                    (got - want).abs() < 1.0e-10,
                    "bandwidth {bw}: {got} vs {want}"
                );
            }
            assert!(m.residual(&x, &b) < RESIDUAL, "bandwidth {bw}");
        }
    }

    #[test]
    fn residual_is_zero_for_an_exact_solution_and_large_for_a_wrong_one() {
        let m = seeded_spd(8, 3);
        let x: Vec<Real> = (0..8).map(|i| 0.5 + 0.1 * i as Real).collect();
        let b = m.mul_vec(&x);
        assert!(m.residual(&x, &b) < RESIDUAL);
        let mut wrong = x.clone();
        wrong[3] += 1.0;
        assert!(m.residual(&wrong, &b) > RESIDUAL);
    }

    #[test]
    fn relative_residual_does_not_depend_on_the_scale_of_the_load() {
        let m = seeded_spd(8, 2);
        let x: Vec<Real> = (0..8).map(|i| 0.5 + 0.1 * i as Real).collect();
        let b = m.mul_vec(&x);
        let reference = m.residual(&x, &b);
        for &factor in &[1.0e-8, 1.0e8] {
            // Scale the solution as well: scaling only the load would leave a
            // large genuine residual rather than testing the same system.
            let scaled_x: Vec<Real> = x.iter().map(|&v| v * factor).collect();
            let scaled_b = m.mul_vec(&scaled_x);
            assert!((m.residual(&scaled_x, &scaled_b) - reference).abs() < 1.0e-12);
        }
    }

    #[test]
    fn non_finite_right_hand_side_is_reported() {
        let factored = BandedMatrix::identity(2)
            .factor()
            .expect("identity factors");
        let mut b = vec![1.0, Real::NAN];
        assert_eq!(
            factored.solve(&mut b),
            Err(SolveError::NonFiniteRightHandSide { equation: 1 })
        );
    }

    #[test]
    fn excessive_relative_residual_is_reported() {
        let mut factored = BandedMatrix::identity(2)
            .factor()
            .expect("identity factors");
        // Perturb the retained unfactored copy so the solve is correct for the
        // factor but no longer satisfies the matrix used for the check.
        factored.a.add_symmetric(0, 0, 1.0);
        let mut b = vec![1.0, 2.0];
        assert_eq!(
            factored.solve(&mut b).err(),
            Some(SolveError::ResidualTooLarge {
                equation: 0,
                residual: 0.5
            })
        );
    }

    #[test]
    fn zero_diagonal_is_reported_as_singular() {
        let mut m = BandedMatrix::new(3, 2);
        m.add_symmetric(0, 0, 1.0);
        m.add_symmetric(1, 1, 0.0);
        m.add_symmetric(2, 2, 1.0);
        assert_eq!(m.factor().err(), Some(SolveError::Singular { equation: 1 }));
    }

    #[test]
    fn negative_pivot_is_reported_as_not_positive_definite() {
        let mut m = BandedMatrix::new(2, 1);
        m.add_symmetric(0, 0, 1.0);
        m.add_symmetric(1, 1, -2.0);
        assert_eq!(
            m.factor().err(),
            Some(SolveError::NotPositiveDefinite { equation: 1 })
        );
    }

    #[test]
    fn non_finite_pivot_is_reported_as_not_finite() {
        let mut m = BandedMatrix::new(2, 1);
        m.add_symmetric(0, 0, 1.0);
        m.add_symmetric(1, 1, 1.0);
        // Write past the assembly-time finite check so the factorisation's
        // own error path is the one under test.
        m.data[0] = Real::NAN;
        assert_eq!(
            m.factor().err(),
            Some(SolveError::NotFinite { equation: 0 })
        );
    }

    #[test]
    fn pivot_near_the_relative_threshold_is_rejected() {
        // After eliminating equation 0 the second pivot is about 1e-13, one
        // hundredth of the original diagonal of equation 1.
        let mut m = BandedMatrix::new(2, 2);
        m.add_symmetric(0, 0, 1.0);
        m.add_symmetric(0, 1, 1.0);
        m.add_symmetric(1, 1, 1.0 + 1.0e-13);
        assert_eq!(m.factor().err(), Some(SolveError::Singular { equation: 1 }));
    }

    #[test]
    fn pivot_rejection_does_not_depend_on_the_units() {
        for &scale in &[1.0e-6, 1.0, 1.0e6] {
            let mut m = BandedMatrix::new(2, 2);
            m.add_symmetric(0, 0, scale);
            m.add_symmetric(0, 1, scale);
            m.add_symmetric(1, 1, scale * (1.0 + 1.0e-13));
            assert_eq!(
                m.factor().err(),
                Some(SolveError::Singular { equation: 1 }),
                "scale {scale}"
            );
        }
    }

    #[test]
    #[should_panic(expected = "rhs length")]
    fn solve_rejects_a_mismatched_rhs() {
        let factored = BandedMatrix::identity(3)
            .factor()
            .expect("identity factors");
        factored.solve(&mut [1.0, 2.0]).expect("length check first");
    }

    #[test]
    #[should_panic(expected = "outside bandwidth")]
    fn connectivity_derived_band_refuses_an_entry_beyond_it() {
        // The band is derived from the elements, and an entry one equation
        // past the widest element is still refused rather than dropped.
        let elements = vec![vec![0, 1, 2, 3]];
        let bw = bandwidth(&elements);
        let mut m = BandedMatrix::new(5, bw);
        m.add_symmetric(0, 4, 1.0);
    }
}
