#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReactantType {
    Atom,
    Linear,
    SphericalTop,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CollisionType {
    AtomAtom,
    AtomLinear,
    LinearLinear,
    AtomSphericalTop,
    LinearSphericalTop,
    SphericalTopSphericalTop,
}

pub fn classify_collision(a: ReactantType, b: ReactantType) -> CollisionType {
    use ReactantType::*;
    use CollisionType::*;

    match (a, b) {
        (Atom, Atom) => AtomAtom,

        (Atom, Linear) | (Linear, Atom) => AtomLinear,

        (Linear, Linear) => LinearLinear,

        (Atom, SphericalTop) | (SphericalTop, Atom) => AtomSphericalTop,

        (Linear, SphericalTop) | (SphericalTop, Linear) => LinearSphericalTop,

        (SphericalTop, SphericalTop) => SphericalTopSphericalTop,
    }
}
