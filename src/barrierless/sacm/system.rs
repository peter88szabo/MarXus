use super::phase_space::phase_space_states;
use super::spcoord::compute_spcoord;
use super::types::{
    SacmCaptureGeometry, SacmEnergyGrid, SacmJRange, SacmPhaseSpaceStates, SacmReactantPair,
    SacmReactantStates, SacmJResolved, SacmSpCoordInput, SacmSpCoordResult,
};

#[derive(Debug)]
pub struct SacmInput {
    pub reactants: SacmReactantPair,
    pub grid: SacmEnergyGrid,
    pub j_range: SacmJRange,
    pub spcoord: Option<SacmSpCoordInput>,
}

#[derive(Debug)]
pub struct SacmPrepared {
    pub reactant_a: SacmReactantStates,
    pub reactant_b: SacmReactantStates,
    pub reactant_a_j: Vec<SacmJResolved>,
    pub reactant_b_j: Vec<SacmJResolved>,
    pub phase_space: SacmPhaseSpaceStates,
    pub spcoord: Option<SacmSpCoordResult>,
    pub capture_geometry: Option<SacmCaptureGeometry>,
}

pub fn prepare_sacm(input: &SacmInput) -> SacmPrepared {
    let spcoord = input.spcoord.as_ref().map(compute_spcoord);
    let capture_geometry = spcoord.as_ref().map(|sp| SacmCaptureGeometry {
        distance_primary: sp.distance_xy,
        distance_secondary: sp.distance_zy,
        angle_deg: sp.angle_deg,
        mass_x: sp.mass_x,
        mass_y: sp.mass_y,
        mass_z: sp.mass_z,
    });
    let reactant_a = input.reactants.reactant_a.rovib_states(input.grid);
    let reactant_b = input.reactants.reactant_b.rovib_states(input.grid);
    let reactant_a_j = input.reactants.reactant_a.j_resolved_states(input.grid, input.j_range);
    let reactant_b_j = input.reactants.reactant_b.j_resolved_states(input.grid, input.j_range);
    let phase_space = phase_space_states(&input.reactants, input.grid);

    SacmPrepared {
        reactant_a,
        reactant_b,
        reactant_a_j,
        reactant_b_j,
        phase_space,
        spcoord,
        capture_geometry,
    }
}
