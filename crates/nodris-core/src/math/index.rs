//! Typed equation indices and the permutation between numberings.
//!
//! Assembly works in the model's original DOF numbering; the banded solver
//! works in a reordered numbering chosen by [`renumber`](super::linalg::renumber). Both
//! are `usize` underneath, so an untyped `Vec<usize>` lets a caller pass one
//! numbering for the other and get a plausible but wrong answer. These
//! newtypes make the direction of every conversion explicit.

/// Equation index in the model's original numbering.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct OldIndex(usize);

/// Equation index in the numbering used by the banded solver.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct NewIndex(usize);

impl OldIndex {
    pub const fn new(value: usize) -> Self {
        Self(value)
    }

    /// The underlying number. Pass an index through
    /// [`Permutation::to_new`] before using it in the solver's numbering.
    pub const fn get(self) -> usize {
        self.0
    }
}

impl NewIndex {
    pub const fn new(value: usize) -> Self {
        Self(value)
    }

    pub const fn get(self) -> usize {
        self.0
    }
}

impl From<usize> for OldIndex {
    fn from(value: usize) -> Self {
        Self(value)
    }
}

impl From<usize> for NewIndex {
    fn from(value: usize) -> Self {
        Self(value)
    }
}

impl From<NewIndex> for usize {
    fn from(value: NewIndex) -> Self {
        value.0
    }
}

// There is deliberately no `From<OldIndex> for usize`. Turning an original
// index into a solver index must go through `Permutation::to_new`.

/// Bijection between the original and the solver's equation numbering.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Permutation {
    /// `to_old[new.get()] == old`
    to_old: Vec<OldIndex>,
    /// `to_new[old.get()] == new`
    to_new: Vec<NewIndex>,
}

impl Permutation {
    /// Build a permutation from `p[new] = old`.
    ///
    /// # Panics
    ///
    /// Panics when `p` is not a bijection of `0..p.len()`.
    pub fn from_old_order(p: &[OldIndex]) -> Self {
        let n = p.len();
        let mut to_new = vec![None; n];
        for (new, &old) in p.iter().enumerate() {
            assert!(
                old.0 < n,
                "permutation maps equation {new} to out-of-range original index {}",
                old.0
            );
            assert!(
                to_new[old.0].is_none(),
                "permutation maps two equations to original index {}",
                old.0
            );
            to_new[old.0] = Some(NewIndex(new));
        }
        let to_new = to_new
            .into_iter()
            .map(|slot| slot.expect("every original index is mapped"))
            .collect();
        Self {
            to_old: p.to_vec(),
            to_new,
        }
    }

    pub fn len(&self) -> usize {
        self.to_old.len()
    }

    pub fn is_empty(&self) -> bool {
        self.to_old.is_empty()
    }

    /// Original index of a solver equation.
    pub fn to_old(&self, new: NewIndex) -> OldIndex {
        self.to_old[new.0]
    }

    /// Solver index of an original equation.
    pub fn to_new(&self, old: OldIndex) -> NewIndex {
        self.to_new[old.0]
    }

    /// Convert a whole element's connectivities into the solver numbering.
    pub fn to_new_all(&self, old: &[OldIndex]) -> Vec<NewIndex> {
        old.iter().map(|&index| self.to_new(index)).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn old(values: &[usize]) -> Vec<OldIndex> {
        values.iter().copied().map(OldIndex::new).collect()
    }

    #[test]
    fn permutation_maps_both_ways() {
        let p = Permutation::from_old_order(&old(&[3, 0, 2, 1]));
        assert_eq!(p.len(), 4);
        for new in 0..4 {
            let o = p.to_old(NewIndex::new(new));
            assert_eq!(p.to_new(o), NewIndex::new(new));
        }
        assert_eq!(p.to_old(NewIndex::new(0)), OldIndex::new(3));
        assert_eq!(p.to_new(OldIndex::new(3)), NewIndex::new(0));
    }

    #[test]
    fn to_new_all_converts_a_whole_element() {
        let p = Permutation::from_old_order(&old(&[3, 0, 2, 1]));
        assert_eq!(
            p.to_new_all(&old(&[3, 1])),
            vec![NewIndex::new(0), NewIndex::new(3)]
        );
    }

    #[test]
    #[should_panic(expected = "out-of-range")]
    fn permutation_rejects_an_out_of_range_index() {
        Permutation::from_old_order(&old(&[0, 2]));
    }

    #[test]
    #[should_panic(expected = "two equations")]
    fn permutation_rejects_a_duplicate_index() {
        Permutation::from_old_order(&old(&[0, 0]));
    }
}
