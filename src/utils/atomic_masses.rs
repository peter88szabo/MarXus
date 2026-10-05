/// Atomic masses (u) of the most abundant isotope of each element, for converting XYZ symbols to masses
/// (moments of inertia, fragment and reduced masses).
///
/// Values: Atomic Mass Evaluation AME2020, M. Wang, W. J. Huang, F. G. Kondev, G. Audi, S. Naimi,
/// Chin. Phys. C 45, 030003 (2021): 1H, 4He, 12C, 14N, 16O, 19F, 20Ne, 31P, 32S, 35Cl, 40Ar, 79Br, 84Kr, 127I,
/// 132Xe. Isotopic (not standard atomic-weight) masses, as used for rotational constants of a specific
/// isotopologue.
///
/// If you need broader coverage, extend this table with isotopic masses (u).
pub fn atomic_mass_amu(symbol: &str) -> Option<f64> {
    // Accept common capitalization variants ("c" -> "C", etc.).
    let s = symbol.trim();
    if s.is_empty() {
        return None;
    }
    let mut chars = s.chars();
    let first = chars.next()?.to_ascii_uppercase();
    let rest: String = chars.as_str().to_ascii_lowercase();
    let normalized = format!("{}{}", first, rest);

    match normalized.as_str() {
        "H" => Some(1.007_825_032_23),
        "He" => Some(4.002_603_254_13),
        "C" => Some(12.0),
        "N" => Some(14.003_074_004_43),
        "O" => Some(15.994_914_619_57),
        "F" => Some(18.998_403_162_73),
        "Ne" => Some(19.992_440_176_2),
        "P" => Some(30.973_761_998_42),
        "S" => Some(31.972_071_174_4),
        "Cl" => Some(34.968_852_682),
        "Ar" => Some(39.962_383_123_7),
        "Br" => Some(78.918_337_6),
        "Kr" => Some(83.911_497_728_2),
        "I" => Some(126.904_471_9),
        "Xe" => Some(131.904_155_085_6),
        _ => None,
    }
}

pub fn mass_vector_from_symbols_amu(symbols: &[String]) -> Result<Vec<f64>, String> {
    let mut out = Vec::with_capacity(symbols.len());
    for (idx, sym) in symbols.iter().enumerate() {
        let m = atomic_mass_amu(sym).ok_or_else(|| {
            format!(
                "Unknown element symbol '{}' at atom index {} (1-based {}).",
                sym,
                idx,
                idx + 1
            )
        })?;
        out.push(m);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_table_holds_the_masses_of_the_most_abundant_isotopes() {
        // Atomic masses of 1H, 12C, 14N and 16O (u), AME2020: Wang et al., Chin. Phys. C 45, 030003 (2021).
        for (symbol, mass) in [("H", 1.007_825_032_23), ("C", 12.0), ("N", 14.003_074_004_43), ("O", 15.994_914_619_57)] {
            assert!((atomic_mass_amu(symbol).unwrap() - mass).abs() < 1e-10, "{symbol}");
        }
        // Monoisotopic mass of H2O, 18.010 564 684 u.
        let water: f64 = mass_vector_from_symbols_amu(&["H".into(), "H".into(), "O".into()]).unwrap().iter().sum();
        assert!((water - 18.010_564_684).abs() < 1e-8, "{water}");
    }

    #[test]
    fn symbols_are_case_insensitive_and_unknown_symbols_are_errors() {
        assert_eq!(atomic_mass_amu("cl"), atomic_mass_amu("Cl"));
        assert!(mass_vector_from_symbols_amu(&["Xx".into()]).is_err());
    }
}
