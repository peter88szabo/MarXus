/// Minimal atomic-mass lookup (amu) for converting XYZ symbols to masses.
///
/// This is intentionally lightweight: it supports the elements we currently
/// encounter in examples and typical combustion/atmospheric chemistry inputs.
///
/// If you need broader coverage, extend this table (keep values in amu).
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

    // Values: conventional standard atomic weights (sufficient precision for inertia).
    match normalized.as_str() {
        "H" => Some(1.007_84),
        "C" => Some(12.0),
        "N" => Some(14.0067),
        "O" => Some(15.999),
        "F" => Some(18.998_403_163),
        "P" => Some(30.973_761_998),
        "S" => Some(32.06),
        "Cl" => Some(35.45),
        "Br" => Some(79.904),
        "I" => Some(126.904_47),
        "Ar" => Some(39.948),
        "He" => Some(4.002_602),
        "Ne" => Some(20.1797),
        "Kr" => Some(83.798),
        "Xe" => Some(131.293),
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
