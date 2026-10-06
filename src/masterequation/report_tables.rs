//! Human-readable tables of master-equation results, in the layout of the MESS output files: every
//! tabulated quantity (a rate coefficient, a yield, ...) is shown
//! - by temperature: for each T a table with the pressures as rows and the quantities as columns;
//! - by pressure: for each p a table with the temperatures as rows and the quantities as columns;
//! - as a temperature-pressure table: for each quantity a table with the pressures as rows and the
//!   temperatures as columns (P\T).
//!
//! Numbers are written with six significant digits, `1.04275e+06`; a missing value (a condition without
//! a result) is written as `***`. Tables wider than `MAX_COLUMNS` quantities are split into several tables
//! with the same rows.

use std::io::Write;

/// Largest number of quantity columns in one table.
pub const MAX_COLUMNS: usize = 8;

/// Width of a number column (11 characters of `1.04275e+06` and the separating spaces).
const COLUMN_WIDTH: usize = 13;

/// A named quantity at every temperature and pressure of a run: `values[t][p]`, None where the condition
/// has no result.
#[derive(Debug, Clone, PartialEq)]
pub struct Quantity {
    pub name: String,
    pub values: Vec<Vec<Option<f64>>>,
}

/// Quantities shown together, under one title and unit.
#[derive(Debug, Clone, PartialEq)]
pub struct QuantityGroup {
    pub title: String,
    pub quantities: Vec<Quantity>,
}

/// A number with six significant digits and a two-digit signed exponent, `1.04275e+06`; `***` for NaN
/// and infinities.
pub fn sci(x: f64) -> String {
    if !x.is_finite() {
        return "***".into();
    }
    let raw = format!("{x:.5e}");
    let (mantissa, exponent) = raw
        .split_once('e')
        .expect("scientific format has an exponent");
    let exponent: i32 = exponent.parse().expect("integer exponent");
    format!(
        "{mantissa}e{}{:02}",
        if exponent < 0 { '-' } else { '+' },
        exponent.abs()
    )
}

/// Temperature or pressure as a row or column label: `300`, `0.1`.
fn condition_label(x: f64) -> String {
    format!("{x}")
}

fn value_text(value: Option<f64>) -> String {
    value.map_or_else(|| "***".to_string(), sci)
}

/// Width of the column of a quantity: wide enough for a number and for its name.
fn column_width(name: &str) -> usize {
    COLUMN_WIDTH.max(name.chars().count() + 2)
}

/// One table: `corner` and the row labels in the first column (right-aligned if `numeric_rows`), the
/// columns right-aligned, followed by an empty line.
fn write_table<W: Write>(
    out: &mut W,
    corner: &str,
    numeric_rows: bool,
    columns: &[&str],
    rows: &[(String, Vec<String>)],
) -> std::io::Result<()> {
    let first = rows
        .iter()
        .map(|(label, _)| label.chars().count())
        .chain([corner.chars().count()])
        .max()
        .unwrap_or(0)
        + 2;
    let widths: Vec<usize> = columns.iter().map(|c| column_width(c)).collect();
    let cell = |text: &str, width: usize| format!("{text:>width$}");
    let label_cell = |text: &str| {
        if numeric_rows {
            format!("{text:>first$}")
        } else {
            format!("{text:<first$}")
        }
    };
    let mut line = label_cell(corner);
    for (column, &width) in columns.iter().zip(&widths) {
        line.push_str(&cell(column, width));
    }
    writeln!(out, "{line}")?;
    for (label, values) in rows {
        let mut line = label_cell(label);
        for (value, &width) in values.iter().zip(&widths) {
            line.push_str(&cell(value, width));
        }
        writeln!(out, "{line}")?;
    }
    writeln!(out)
}

/// For each temperature a table: pressures (Torr) as rows, the quantities of the group as columns.
pub fn write_tables_by_temperature<W: Write>(
    out: &mut W,
    temperatures: &[f64],
    pressures: &[f64],
    group: &QuantityGroup,
) -> std::io::Result<()> {
    writeln!(out, "{}, by temperature:\n", group.title)?;
    for (t, &temperature) in temperatures.iter().enumerate() {
        writeln!(out, "   Temperature = {} K\n", condition_label(temperature))?;
        for chunk in group.quantities.chunks(MAX_COLUMNS) {
            let columns: Vec<&str> = chunk.iter().map(|q| q.name.as_str()).collect();
            let rows: Vec<(String, Vec<String>)> = pressures
                .iter()
                .enumerate()
                .map(|(p, &pressure)| {
                    (
                        condition_label(pressure),
                        chunk.iter().map(|q| value_text(q.values[t][p])).collect(),
                    )
                })
                .collect();
            write_table(out, "P(torr)", true, &columns, &rows)?;
        }
    }
    Ok(())
}

/// For each pressure a table: temperatures (K) as rows, the quantities of the group as columns.
pub fn write_tables_by_pressure<W: Write>(
    out: &mut W,
    temperatures: &[f64],
    pressures: &[f64],
    group: &QuantityGroup,
) -> std::io::Result<()> {
    writeln!(out, "{}, by pressure:\n", group.title)?;
    for (p, &pressure) in pressures.iter().enumerate() {
        writeln!(out, "   Pressure = {} torr\n", condition_label(pressure))?;
        for chunk in group.quantities.chunks(MAX_COLUMNS) {
            let columns: Vec<&str> = chunk.iter().map(|q| q.name.as_str()).collect();
            let rows: Vec<(String, Vec<String>)> = temperatures
                .iter()
                .enumerate()
                .map(|(t, &temperature)| {
                    (
                        condition_label(temperature),
                        chunk.iter().map(|q| value_text(q.values[t][p])).collect(),
                    )
                })
                .collect();
            write_table(out, "T(K)", true, &columns, &rows)?;
        }
    }
    Ok(())
}

/// For each quantity of the group a table: pressures (Torr) as rows, temperatures (K) as columns.
pub fn write_temperature_pressure_tables<W: Write>(
    out: &mut W,
    temperatures: &[f64],
    pressures: &[f64],
    group: &QuantityGroup,
) -> std::io::Result<()> {
    writeln!(out, "{}, temperature-pressure tables:\n", group.title)?;
    let temperature_labels: Vec<String> =
        temperatures.iter().map(|&t| condition_label(t)).collect();
    let columns: Vec<&str> = temperature_labels.iter().map(|s| s.as_str()).collect();
    for quantity in &group.quantities {
        writeln!(out, "{}\n", quantity.name)?;
        let rows: Vec<(String, Vec<String>)> = pressures
            .iter()
            .enumerate()
            .map(|(p, &pressure)| {
                (
                    condition_label(pressure),
                    (0..temperatures.len())
                        .map(|t| value_text(quantity.values[t][p]))
                        .collect(),
                )
            })
            .collect();
        write_table(out, "P\\T", false, &columns, &rows)?;
    }
    Ok(())
}

/// A table with row labels and named columns, `corner` in the top-left cell (e.g. `From\To`).
pub fn write_labelled_table<W: Write>(
    out: &mut W,
    corner: &str,
    columns: &[String],
    rows: &[(String, Vec<Option<f64>>)],
) -> std::io::Result<()> {
    let columns: Vec<&str> = columns.iter().map(|c| c.as_str()).collect();
    let rows: Vec<(String, Vec<String>)> = rows
        .iter()
        .map(|(label, values)| {
            (
                label.clone(),
                values.iter().map(|&v| value_text(v)).collect(),
            )
        })
        .collect();
    write_table(out, corner, false, &columns, &rows)
}

/// Machine-readable form of the groups: for each group a block `# <title>`, a header
/// `T[K],P[Torr],<quantity names>` and one row per condition (temperatures outer, pressures inner), empty
/// fields for missing values, blocks separated by an empty line.
pub fn write_groups_csv<W: Write>(
    out: &mut W,
    temperatures: &[f64],
    pressures: &[f64],
    groups: &[QuantityGroup],
) -> std::io::Result<()> {
    for group in groups {
        writeln!(out, "# {}", group.title)?;
        let names: Vec<String> = group
            .quantities
            .iter()
            .map(|q| q.name.replace(',', ";"))
            .collect();
        writeln!(out, "T[K],P[Torr],{}", names.join(","))?;
        for (t, &temperature) in temperatures.iter().enumerate() {
            for (p, &pressure) in pressures.iter().enumerate() {
                let values: Vec<String> = group
                    .quantities
                    .iter()
                    .map(|q| q.values[t][p].map_or_else(String::new, |v| format!("{v:.6e}")))
                    .collect();
                writeln!(
                    out,
                    "{},{},{}",
                    condition_label(temperature),
                    condition_label(pressure),
                    values.join(",")
                )?;
            }
        }
        writeln!(out)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn group(count: usize) -> QuantityGroup {
        // values[t][p] = 1000 t + p + quantity index / 10, with one missing condition (t = 1, p = 0).
        QuantityGroup {
            title: "Test quantities (1/s)".into(),
            quantities: (0..count)
                .map(|q| Quantity {
                    name: format!("G{q}->P"),
                    values: (0..2)
                        .map(|t| {
                            (0..3)
                                .map(|p| {
                                    if t == 1 && p == 0 {
                                        None
                                    } else {
                                        Some(1000.0 * t as f64 + p as f64 + q as f64 / 10.0)
                                    }
                                })
                                .collect()
                        })
                        .collect(),
                })
                .collect(),
        }
    }

    fn text(write: impl FnOnce(&mut Vec<u8>) -> std::io::Result<()>) -> String {
        let mut out = Vec::new();
        write(&mut out).unwrap();
        String::from_utf8(out).unwrap()
    }

    /// Lines of the table that starts with a line beginning with `header`, up to the next empty line.
    fn table<'a>(text: &'a str, header: &str) -> Vec<&'a str> {
        let lines: Vec<&str> = text.lines().collect();
        let start = lines
            .iter()
            .position(|l| l.trim_start().starts_with(header))
            .expect(header);
        lines[start..]
            .iter()
            .take_while(|l| !l.trim().is_empty())
            .copied()
            .collect()
    }

    #[test]
    fn numbers_have_six_significant_digits_and_a_signed_two_digit_exponent() {
        assert_eq!(sci(1_042_750.0), "1.04275e+06");
        assert_eq!(sci(-2.91554e-10), "-2.91554e-10");
        assert_eq!(sci(0.5), "5.00000e-01");
        assert_eq!(sci(0.0), "0.00000e+00");
        assert_eq!(sci(9.999_996e5), "1.00000e+06");
        assert_eq!(sci(1.5e123), "1.50000e+123");
        assert_eq!(sci(f64::NAN), "***");
        assert_eq!(sci(f64::INFINITY), "***");
    }

    #[test]
    fn tables_by_temperature_have_the_pressures_as_rows() {
        let t = text(|o| {
            write_tables_by_temperature(o, &[300.0, 400.0], &[10.0, 100.0, 760.0], &group(2))
        });
        assert!(
            t.contains("Temperature = 300 K") && t.contains("Temperature = 400 K"),
            "{t}"
        );
        let rows = table(&t, "P(torr)");
        assert_eq!(rows.len(), 4, "{t}");
        let header: Vec<&str> = rows[0].split_whitespace().collect();
        assert_eq!(header, ["P(torr)", "G0->P", "G1->P"]);
        let first: Vec<&str> = rows[1].split_whitespace().collect();
        assert_eq!(first, ["10", "0.00000e+00", "1.00000e-01"]);
        // All lines of a table have the same width (right-aligned columns).
        assert!(rows.iter().all(|r| r.len() == rows[0].len()), "{t}");
    }

    #[test]
    fn tables_by_pressure_have_the_temperatures_as_rows_and_mark_missing_values() {
        let t = text(|o| {
            write_tables_by_pressure(o, &[300.0, 400.0], &[10.0, 100.0, 760.0], &group(2))
        });
        assert!(
            t.contains("Pressure = 10 torr") && t.contains("Pressure = 760 torr"),
            "{t}"
        );
        let rows = table(&t, "T(K)");
        assert_eq!(rows.len(), 3, "{t}");
        let missing: Vec<&str> = rows[2].split_whitespace().collect();
        assert_eq!(missing, ["400", "***", "***"], "{t}");
    }

    #[test]
    fn temperature_pressure_tables_have_one_table_per_quantity() {
        let t = text(|o| {
            write_temperature_pressure_tables(o, &[300.0, 400.0], &[10.0, 100.0, 760.0], &group(2))
        });
        assert_eq!(
            t.lines()
                .filter(|l| l.trim_start().starts_with("P\\T"))
                .count(),
            2,
            "{t}"
        );
        let rows = table(&t, "P\\T");
        let header: Vec<&str> = rows[0].split_whitespace().collect();
        assert_eq!(header, ["P\\T", "300", "400"]);
        let last: Vec<&str> = rows[3].split_whitespace().collect();
        assert_eq!(last, ["760", "2.00000e+00", "1.00200e+03"]);
        assert!(t.contains("G0->P") && t.contains("G1->P"), "{t}");
    }

    #[test]
    fn wide_groups_are_split_into_tables_of_at_most_max_columns() {
        let t = text(|o| {
            write_tables_by_pressure(o, &[300.0, 400.0], &[10.0], &group(MAX_COLUMNS + 3))
        });
        let headers: Vec<&str> = t
            .lines()
            .filter(|l| l.trim_start().starts_with("T(K)"))
            .collect();
        assert_eq!(headers.len(), 2, "{t}");
        assert_eq!(headers[0].split_whitespace().count(), MAX_COLUMNS + 1);
        assert_eq!(headers[1].split_whitespace().count(), 4);
    }

    #[test]
    fn long_names_widen_their_column() {
        let mut g = group(1);
        g.quantities[0].name = "stabilization(G4)".into();
        let t = text(|o| write_tables_by_pressure(o, &[300.0, 400.0], &[10.0], &g));
        let rows = table(&t, "T(K)");
        assert!(rows.iter().all(|r| r.len() == rows[0].len()), "{t}");
        assert!(rows[0].ends_with("stabilization(G4)"), "{t}");
    }

    #[test]
    fn groups_are_written_as_titled_csv_blocks() {
        let t = text(|o| {
            write_groups_csv(
                o,
                &[300.0, 400.0],
                &[10.0, 100.0, 760.0],
                &[group(2), group(1)],
            )
        });
        let blocks: Vec<&str> = t.split("\n\n").filter(|b| !b.trim().is_empty()).collect();
        assert_eq!(blocks.len(), 2, "{t}");
        let lines: Vec<&str> = blocks[0].lines().collect();
        assert_eq!(lines[0], "# Test quantities (1/s)");
        assert_eq!(lines[1], "T[K],P[Torr],G0->P,G1->P");
        assert_eq!(lines.len(), 2 + 6);
        assert_eq!(lines[2], "300,10,0.000000e0,1.000000e-1");
        // Missing condition (t = 1, p = 0): empty fields.
        assert_eq!(lines[5], "400,10,,");
    }

    #[test]
    fn labelled_tables_have_a_corner_label_and_aligned_rows() {
        let columns = vec!["G2".to_string(), "escape(G4)".to_string()];
        let rows = vec![
            ("G2".to_string(), vec![Some(1.0e6), None]),
            ("R".to_string(), vec![Some(8.0e-12), Some(-1.3e-13)]),
        ];
        let t = text(|o| write_labelled_table(o, "From\\To", &columns, &rows));
        let lines = table(&t, "From\\To");
        assert_eq!(lines.len(), 3, "{t}");
        assert!(lines.iter().all(|r| r.len() == lines[0].len()), "{t}");
        assert_eq!(
            lines[2].split_whitespace().collect::<Vec<_>>(),
            ["R", "8.00000e-12", "-1.30000e-13"]
        );
        assert_eq!(
            lines[1].split_whitespace().collect::<Vec<_>>(),
            ["G2", "1.00000e+06", "***"]
        );
    }
}
