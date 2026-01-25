use crate::molecule::MoleculeStruct;
use crate::rrkm::sum_and_density::{JResolvedStates, RotorSymmetry};

use super::centrifugal::SacmCentrifugal;

#[derive(Debug)]
pub struct SacmReactant {
    pub molecule: MoleculeStruct,
    pub rotor: RotorSymmetry,
    pub centrifugal: Option<SacmCentrifugal>,
}

#[derive(Debug)]
pub struct SacmReactantPair {
    pub reactant_a: SacmReactant,
    pub reactant_b: SacmReactant,
}

#[derive(Debug)]
pub struct SacmFragments {
    pub fragment_a: MoleculeStruct,
    pub fragment_b: MoleculeStruct,
}

#[derive(Debug, Clone, Copy)]
pub struct SacmEnergyGrid {
    pub dE: f64,
    pub emax: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct SacmJRange {
    pub j_start: usize,
    pub j_end: usize,
    pub j_step: usize,
}

#[derive(Debug)]
pub struct SacmJResolved {
    pub j: usize,
    pub states: JResolvedStates,
}

#[derive(Debug)]
pub struct SacmReactantStates {
    pub rho_e: Vec<f64>,
    pub we: Vec<f64>,
}

#[derive(Debug)]
pub struct SacmPhaseSpaceStates {
    pub rho_ts: Vec<f64>,
    pub we_ts: Vec<f64>,
}

pub type SacmSpCoordInput = super::spcoord::SpCoordInput;
pub type SacmSpCoordResult = super::spcoord::SpCoordResult;

#[derive(Debug, Clone)]
pub struct SacmCaptureGeometry {
    pub distance_primary: f64,
    pub distance_secondary: Option<f64>,
    pub angle_deg: Option<f64>,
    pub mass_x: f64,
    pub mass_y: f64,
    pub mass_z: Option<f64>,
}

#[derive(Debug, Clone, Copy)]
pub enum ReactantRotorCase {
    Linear,
    SphericalTop,
    OblateSymmetricTop,
    ProlateSymmetricTop,
}

impl ReactantRotorCase {
    pub fn famc_index(self) -> usize {
        match self {
            ReactantRotorCase::Linear => 0,
            ReactantRotorCase::SphericalTop => 1,
            ReactantRotorCase::OblateSymmetricTop => 2,
            ReactantRotorCase::ProlateSymmetricTop => 3,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub enum ProductRotorCase {
    LinearAtom,
    SphericalAtom,
    LinearLinear,
    LinearSpherical,
    SphericalSpherical,
}

impl ProductRotorCase {
    pub fn famc_index(self) -> usize {
        match self {
            ProductRotorCase::LinearAtom => 0,
            ProductRotorCase::SphericalAtom => 1,
            ProductRotorCase::LinearLinear => 2,
            ProductRotorCase::LinearSpherical => 3,
            ProductRotorCase::SphericalSpherical => 4,
        }
    }
}
