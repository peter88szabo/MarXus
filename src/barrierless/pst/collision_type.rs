#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReactantType {
    Atom,
    Linear,
    Spherical,
    Prolate,
    Oblate,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CollisionType {
    AtomAtom,

    AtomLinear,
    LinearLinear,

    AtomSpherical,
    LinearSpherical,
    SphericalSpherical,

    AtomProlate,
    LinearProlate,
    SphericalProlate,
    ProlateProlate,

    AtomOblate,
    LinearOblate,
    SphericalOblate,
    OblateOblate,

    ProlateOblate,
}

pub fn classify_collision(a: ReactantType, b: ReactantType) -> CollisionType {
    use ReactantType::*;
    use CollisionType::*;

    match (a, b) {
        (Atom, Atom) => AtomAtom,

        (Atom, Linear) | (Linear, Atom) => AtomLinear,
        (Atom, Spherical) | (Spherical, Atom) => AtomSpherical,
        (Atom, Prolate) | (Prolate, Atom) => AtomProlate,
        (Atom, Oblate) | (Oblate, Atom) => AtomOblate,

        (Linear, Linear) => LinearLinear,
        (Linear, Spherical) | (Spherical, Linear) => LinearSpherical,
        (Linear, Prolate) | (Prolate, Linear) => LinearProlate,
        (Linear, Oblate) | (Oblate, Linear) => LinearOblate,

        (Spherical, Spherical) => SphericalSpherical,
        (Spherical, Prolate) | (Prolate, Spherical) => SphericalProlate,
        (Spherical, Oblate) | (Oblate, Spherical) => SphericalOblate,

        (Prolate, Prolate) => ProlateProlate,
        (Oblate, Oblate) => OblateOblate,
        (Prolate, Oblate) | (Oblate, Prolate) => ProlateOblate,
    }
}
