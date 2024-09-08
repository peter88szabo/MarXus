use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;

fn main() -> io::Result<()> {
    let path = Path::new("Case2_starting_from_ZZ_allyl+O2_version-0.1.inp");
    let file = File::open(&path)?;
    let reader = io::BufReader::new(file);

    let mut species_data = Vec::new();
    let mut current_atoms: Vec<String> = Vec::new();
    let mut current_coordinates: Vec<f64> = Vec::new();
    let mut in_species_section = false;
    let mut in_geometry_section = false;

    for line in reader.lines() {
        let line = line?;

        // Check for start of a new species section
        if line.contains("Bimolecular") || line.contains("Well") || line.contains("Barrier") {
             if in_species_section {
                // Save previous species data
                species_data.push((current_atoms.clone(), current_coordinates.clone()));
                current_atoms.clear();
                current_coordinates.clear();
            }
            in_species_section = true;
        }

        // Check for geometry section within a species
        if line.contains("Geometry") {
            in_geometry_section = true;
            continue;
        } else if line.contains("GeometryEnd") && in_geometry_section {
            in_geometry_section = false;
            continue;
        }

        //if in_species_section && in_geometry_section {
        if in_geometry_section {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 4 {
                current_atoms.push(parts[0].to_string());
                for i in 1..4 {
                    let coord: f64 = parts[i].parse().unwrap_or(0.0); // Replace with proper error handling
                    current_coordinates.push(coord);
                }
            }
        }
    }

    // Add the last species data if it exists
    if !current_atoms.is_empty() || !current_coordinates.is_empty() {
        species_data.push((current_atoms, current_coordinates));
    }

    // Print parsed data for each species
    for (i, (atoms, coordinates)) in species_data.iter().enumerate() {
        println!();
        println!("Species {}: Atoms: {:?}, Coordinates: {:?}", i + 1, atoms, coordinates);
    }

    Ok(())
}



