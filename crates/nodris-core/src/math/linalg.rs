//! Banded storage for assembled stiffness matrices, with the graph helpers
//! that decide how wide the band has to be.
//!
//! The lower triangle is stored row-wise: `A(i, j)` with `i >= j` sits at
//! `data[i * bw + (i - j)]`, so walking along a row is a contiguous scan.
//! The upper triangle sits in the later row that owns the symmetric pair.

use super::Real;
use std::collections::VecDeque;

/// Symmetric matrix stored by lower-triangle rows in a band of half-width `bw`.
#[derive(Clone, Debug, PartialEq)]
pub struct BandedMatrix {
    n: usize,
    bw: usize,
    data: Vec<Real>,
}

impl BandedMatrix {
    /// Zero matrix of dimension `n` and half-bandwidth `bw` (at least 1).
    pub fn new(n: usize, bw: usize) -> Self {
        let bw = bw.max(1);
        Self {
            n,
            bw,
            data: vec![0.0; n * bw],
        }
    }

    /// Identity matrix of dimension `n`.
    pub fn identity(n: usize) -> Self {
        let mut m = Self::new(n, 1);
        for i in 0..n {
            m.add(i, i, 1.0);
        }
        m
    }

    /// Matrix dimension.
    pub fn dim(&self) -> usize {
        self.n
    }

    /// Half-bandwidth: entries with `|i - j| < bandwidth` are stored.
    pub fn bandwidth(&self) -> usize {
        self.bw
    }

    /// Number of coefficient slots, `n * bandwidth`.
    pub fn stored(&self) -> usize {
        self.data.len()
    }

    /// Row-major lower-triangle slot for `i >= j`.
    #[inline]
    fn idx(i: usize, j: usize, bw: usize) -> usize {
        i * bw + (i - j)
    }

    /// Adds `v` to the symmetric entry `(i, j)`.
    ///
    /// # Panics
    /// Panics when `i` or `j` is outside the matrix dimensions, or when the
    /// entry lies outside the stored half-band. The band check runs in every
    /// build profile: a dropped stiffness contribution would yield plausible
    /// but wrong displacements, so an out-of-band write is refused rather
    /// than silently ignored.
    pub fn add(&mut self, i: usize, j: usize, v: Real) {
        assert!(i < self.n && j < self.n, "banded index out of range");
        let (hi, lo) = if i >= j { (i, j) } else { (j, i) };
        assert!(
            hi - lo < self.bw,
            "entry ({i},{j}) outside bandwidth {}",
            self.bw
        );
        self.data[Self::idx(hi, lo, self.bw)] += v;
    }

    /// Reads `A(i, j)`, returning zero outside the band.
    pub fn get(&self, i: usize, j: usize) -> Real {
        assert!(i < self.n && j < self.n, "banded index out of range");
        let (hi, lo) = if i >= j { (i, j) } else { (j, i) };
        if hi - lo >= self.bw {
            0.0
        } else {
            self.data[Self::idx(hi, lo, self.bw)]
        }
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
            let first = (i + 1).saturating_sub(self.bw);
            for j in first..i {
                sum += self.data[i * self.bw + (i - j)] * rhs[j];
            }
            // Upper entries live in the later row that owns the symmetric pair.
            let last = (i + self.bw).min(self.n);
            for j in i + 1..last {
                sum += self.data[j * self.bw + (j - i)] * rhs[j];
            }
            out[i] = sum;
        }
        out
    }

    /// Expands the band to a dense row-major matrix for inspection and tests.
    pub fn to_dense(&self) -> Vec<Vec<Real>> {
        let mut dense = vec![vec![0.0; self.n]; self.n];
        for i in 0..self.n {
            let first = i.saturating_sub(self.bw - 1);
            for j in first..=i {
                dense[i][j] = self.data[Self::idx(i, j, self.bw)];
                dense[j][i] = dense[i][j];
            }
        }
        dense
    }
}

/// Reverse Cuthill-McKee ordering.
///
/// Returns a permutation `p` where `p[new_index] = old_index`. The ordering
/// minimises the envelope of the matrix, which is exactly the bandwidth the
/// banded solver pays for.
///
/// `adj` must have one entry per node: `adj[u]` lists the nodes connected to
/// `u`. Disconnected components are ordered one after another.
pub fn rcm(n: usize, adj: &[Vec<usize>]) -> Vec<usize> {
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
    use super::{bandwidth, rcm, BandedMatrix, Real};

    fn next(seed: &mut u64) -> Real {
        *seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((*seed >> 11) as Real) / ((1u64 << 53) as Real)
    }

    fn seeded_spd(n: usize, bw: usize) -> BandedMatrix {
        let mut seed = 0x9e37_79b9_7f4a_7c15;
        let mut m = BandedMatrix::new(n, bw);
        for i in 0..n {
            m.add(i, i, 2.0 + next(&mut seed));
            for j in i + 1..(i + bw).min(n) {
                m.add(i, j, 0.2 * next(&mut seed));
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

    #[test]
    fn add_reports_both_halves_of_a_symmetric_entry() {
        let mut m = BandedMatrix::new(6, 2);
        m.add(2, 1, 3.0);
        assert_eq!(m.get(2, 1), 3.0);
        assert_eq!(m.get(1, 2), 3.0);
        m.add(1, 1, 2.0);
        m.add(1, 1, 3.0);
        assert_eq!(m.get(1, 1), 5.0);
        assert_eq!(m.get(0, 5), 0.0);
    }

    #[test]
    fn dense_expansion_keeps_the_band_and_symmetry() {
        let mut m = BandedMatrix::new(4, 2);
        m.add(0, 0, 1.0);
        m.add(1, 0, 2.0);
        m.add(1, 1, 3.0);
        m.add(3, 2, 4.0);
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
        m.add(0, 4, 1.0);
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
                    m.add(el[i], el[j], 1.0);
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
    #[should_panic(expected = "outside bandwidth")]
    fn connectivity_derived_band_refuses_an_entry_beyond_it() {
        // The band is derived from the elements, and an entry one equation
        // past the widest element is still refused rather than dropped.
        let elements = vec![vec![0, 1, 2, 3]];
        let bw = bandwidth(&elements);
        let mut m = BandedMatrix::new(5, bw);
        m.add(0, 4, 1.0);
    }
}
