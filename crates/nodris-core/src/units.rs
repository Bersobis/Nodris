//! Unit boundaries for the kernel.
//!
//! Analysis quantities live in kN, m and kN/m². Section tables quote
//! dimensions in mm, mm², mm⁴ and mm⁶, and material strengths in MPa. Those
//! values cross into analysis units **once**, through the helpers here, so a
//! table lookup cannot silently mix mm⁴ into an m⁴ stiffness.
//!
//! The constants are named after the unit they are written in, not after an
//! abstract "scale". [`MM`] is metres per millimetre; [`KN`] is newtons per
//! kilonewton; [`MPA`] is pascals per megapascal. Analysis stress is kN/m²,
//! so [`MPA_TO_KN_M2`] connects the material table to the solver.

use crate::math::Real;

/// Metres per millimetre.
pub const MM: Real = 1.0e-3;

/// Newtons per kilonewton.
pub const KN: Real = 1.0e3;

/// Pascals (N/m²) per megapascal.
pub const MPA: Real = 1.0e6;

/// Analysis stress units (kN/m²) per MPa.
///
/// One MPa is 10⁶ N/m², which is 10³ kN/m². The kernel stores stress as
/// force per analysis length squared, so this factor, not [`MPA`], is what a
/// material strength needs before it reaches a check.
pub const MPA_TO_KN_M2: Real = MPA / KN;

const MM2: Real = MM * MM;
const MM4: Real = MM2 * MM2;
const MM6: Real = MM4 * MM2;

/// Millimetres → metres.
pub const fn length_from_mm(mm: Real) -> Real {
    mm * MM
}

/// Metres → millimetres.
pub const fn length_to_mm(m: Real) -> Real {
    m / MM
}

/// Kilonewtons → newtons.
pub const fn force_from_kn(kn: Real) -> Real {
    kn * KN
}

/// Newtons → kilonewtons.
pub const fn force_to_kn(n: Real) -> Real {
    n / KN
}

/// Megapascals → analysis stress (kN/m²).
pub const fn stress_from_mpa(mpa: Real) -> Real {
    mpa * MPA_TO_KN_M2
}

/// Analysis stress (kN/m²) → megapascals.
pub const fn stress_to_mpa(kn_m2: Real) -> Real {
    kn_m2 / MPA_TO_KN_M2
}

/// Millimetres² → metres². Use once when a table area enters the model.
pub const fn area_from_mm2(mm2: Real) -> Real {
    mm2 * MM2
}

/// Metres² → millimetres², for reporting a model area in table units.
pub const fn area_to_mm2(m2: Real) -> Real {
    m2 / MM2
}

/// Millimetres⁴ → metres⁴, for second moments and torsional constants.
pub const fn inertia_from_mm4(mm4: Real) -> Real {
    mm4 * MM4
}

/// Metres⁴ → millimetres⁴, for reporting a model inertia in table units.
pub const fn inertia_to_mm4(m4: Real) -> Real {
    m4 / MM4
}

/// Millimetres⁶ → metres⁶, for warping constants.
pub const fn warping_from_mm6(mm6: Real) -> Real {
    mm6 * MM6
}

/// Metres⁶ → millimetres⁶, for reporting a warping constant in table units.
pub const fn warping_to_mm6(m6: Real) -> Real {
    m6 / MM6
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(got: Real, want: Real) -> bool {
        (got - want).abs() <= 1.0e-12 * want.abs().max(1.0)
    }

    #[test]
    fn boundary_constants_state_their_units() {
        assert!(close(MM, 0.001));
        assert!(close(KN, 1000.0));
        assert!(close(MPA, 1.0e6));
        assert!(close(MPA_TO_KN_M2, 1000.0));
    }

    #[test]
    fn length_force_and_stress_convert_to_analysis_units() {
        assert!(close(length_from_mm(2500.0), 2.5));
        assert!(close(length_to_mm(2.5), 2500.0));
        assert!(close(force_from_kn(10.0), 10_000.0));
        assert!(close(force_to_kn(10_000.0), 10.0));
        // 235 MPa is 235_000 kN/m² in the analysis stress unit.
        assert!(close(stress_from_mpa(235.0), 235_000.0));
        assert!(close(stress_to_mpa(235_000.0), 235.0));
    }

    #[test]
    fn area_round_trips_through_mm2_and_m2() {
        // 1 m² is exactly 10⁶ mm².
        assert!(close(area_from_mm2(1.0e6), 1.0));
        for mm2 in [0.0, 1.0, 1_000.0, 1.0e6, 2.5e8] {
            assert!(close(area_to_mm2(area_from_mm2(mm2)), mm2));
        }
    }

    #[test]
    fn inertia_round_trips_through_mm4_and_m4() {
        // 1 m⁴ is exactly 10¹² mm⁴.
        assert!(close(inertia_from_mm4(1.0e12), 1.0));
        for mm4 in [0.0, 1.0, 1.0e4, 1.0e12, 8.0e14] {
            assert!(close(inertia_to_mm4(inertia_from_mm4(mm4)), mm4));
        }
    }

    #[test]
    fn warping_round_trips_through_mm6_and_m6() {
        // 1 m⁶ is exactly 10¹⁸ mm⁶.
        assert!(close(warping_from_mm6(1.0e18), 1.0));
        for mm6 in [0.0, 1.0, 1.0e6, 1.0e18, 3.0e21] {
            assert!(close(warping_to_mm6(warping_from_mm6(mm6)), mm6));
        }
    }

    #[test]
    fn section_properties_do_not_share_a_power() {
        // The same table number means different areas, inertias and warping
        // constants; a shared scale would make one of these fail.
        let value = 1.0;
        assert!(area_from_mm2(value) > inertia_from_mm4(value));
        assert!(inertia_from_mm4(value) > warping_from_mm6(value));
    }
}
