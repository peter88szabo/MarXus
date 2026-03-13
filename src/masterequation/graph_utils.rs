//! Small, general-purpose graph utilities used by master-equation network builders.
//!
//! These are intentionally physics-agnostic: they operate on node indices and edges.
//! Keeping them separate helps keep the chemistry/physics code readable.

#[derive(Clone, Debug)]
pub struct DisjointSetUnion {
    parent: Vec<usize>,
    rank: Vec<u8>,
}

impl DisjointSetUnion {
    pub fn new(node_count: usize) -> Self {
        Self {
            parent: (0..node_count).collect(),
            rank: vec![0; node_count],
        }
    }

    pub fn find(&mut self, x: usize) -> usize {
        let p = self.parent[x];
        if p == x {
            return x;
        }
        let root = self.find(p);
        self.parent[x] = root;
        root
    }

    pub fn union(&mut self, a: usize, b: usize) {
        let mut ra = self.find(a);
        let mut rb = self.find(b);
        if ra == rb {
            return;
        }
        if self.rank[ra] < self.rank[rb] {
            std::mem::swap(&mut ra, &mut rb);
        }
        self.parent[rb] = ra;
        if self.rank[ra] == self.rank[rb] {
            self.rank[ra] = self.rank[ra].saturating_add(1);
        }
    }

    pub fn component_roots(&mut self) -> Vec<usize> {
        let mut roots: Vec<usize> = (0..self.parent.len()).map(|i| self.find(i)).collect();
        roots.sort_unstable();
        roots.dedup();
        roots
    }
}
