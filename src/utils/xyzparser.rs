pub fn parse_xyz(xyz: &str) -> (usize, Vec<String>, Vec<Vec<f64>>) {
    // Split the XYZ string into lines and filter out any empty or whitespace-only lines
    let lines: Vec<&str> = xyz.lines().filter(|line| !line.trim().is_empty()).collect();

    // Count the number of atoms from the non-empty lines
    let natoms = lines.len();

    let mut atoms = Vec::new();
    let mut qxyz = Vec::with_capacity(natoms); // Initialize qxyz with N rows

    // Loop through each line to extract atomic coordinates
    for line in lines {
        let parts: Vec<&str> = line.split_whitespace().collect();

        // Parse the atomic symbol
        let atom = parts[0].to_string();
        atoms.push(atom);

        // Parse the x, y, z coordinates and store them in a sub-vector
        let x = parts[1].parse::<f64>().unwrap();
        let y = parts[2].parse::<f64>().unwrap();
        let z = parts[3].parse::<f64>().unwrap();
        qxyz.push(vec![x, y, z]);
    }

    (natoms, atoms, qxyz)
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_xyz() {
        let xyz_data = "
          C      -1.185385      1.500364     -0.174799
          C       0.057100      1.525641      0.290486
          C       1.182421      0.671307     -0.090281
          C       1.182420     -0.671306     -0.090282
          C       0.057100     -1.525640      0.290485
          C      -1.185388     -1.500364     -0.174796
          H       0.285347      2.226153      1.090336
          H       2.144030      1.164815     -0.199258
          H       2.144029     -1.164815     -0.199261
          H       0.285350     -2.226154      1.090332
          H      -1.972016     -2.076478      0.294389
          H      -1.462109     -0.924082     -1.042658
          H      -1.972017      2.076476      0.294382
          H      -1.462101      0.924083     -1.042665
        ";

        let (natoms, atoms, qxyz) = parse_xyz(xyz_data);

        assert_eq!(natoms, 14);
        assert_eq!(atoms.len(), 14);
        assert_eq!(qxyz.len(), 14);
        assert_eq!(atoms[0], "C");
        assert_eq!(atoms[13], "H");

        // Verify specific coordinates
        assert!((qxyz[0][0] - -1.185385).abs() < 1e-6);
        assert!((qxyz[0][1] - 1.500364).abs() < 1e-6);
        assert!((qxyz[0][2] - -0.174799).abs() < 1e-6);
        assert!((qxyz[13][0] - -1.462101).abs() < 1e-6);
        assert!((qxyz[13][1] - 0.924083).abs() < 1e-6);
        assert!((qxyz[13][2] - -1.042665).abs() < 1e-6);
    }
}
