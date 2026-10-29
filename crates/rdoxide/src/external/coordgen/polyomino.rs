//! The hexagon lattice shapes of `CoordgenMacrocycleBuilder`: `hexCoords`, `vertexCoords`, `Hex` and `Polyomino`
//! (shapes made of hexagons whose contour is the path the atoms of a macrocycle are placed on).

use std::collections::HashSet;
use std::ops::{Add, Sub};

/// `hexCoords`
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct HexCoords {
    pub x: i32,
    pub y: i32,
}

impl HexCoords {
    pub fn new(x: i32, y: i32) -> Self {
        HexCoords { x, y }
    }

    pub fn distance_from(&self, origin: HexCoords) -> i32 {
        let dx = (self.x - origin.x).abs();
        let dy = (self.y - origin.y).abs();
        let dz = ((-self.x - self.y) - (-origin.x - origin.y)).abs();
        let max = dx.max(dy);
        dz.max(max)
    }

    pub fn rotate_30_degrees(&self) -> HexCoords {
        let z = -self.x - self.y;
        HexCoords { x: -z, y: -self.x }
    }

    pub fn to_vertex_coords(&self) -> VertexCoords {
        VertexCoords {
            x: self.x,
            y: self.y,
            z: -self.x - self.y,
        }
    }
}

/// `vertexCoords`
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct VertexCoords {
    pub x: i32,
    pub y: i32,
    pub z: i32,
}

impl VertexCoords {
    pub fn new(x: i32, y: i32, z: i32) -> Self {
        VertexCoords { x, y, z }
    }
}

impl Add for VertexCoords {
    type Output = VertexCoords;
    fn add(self, o: VertexCoords) -> VertexCoords {
        VertexCoords::new(self.x + o.x, self.y + o.y, self.z + o.z)
    }
}

impl Sub for VertexCoords {
    type Output = VertexCoords;
    fn sub(self, o: VertexCoords) -> VertexCoords {
        VertexCoords::new(self.x - o.x, self.y - o.y, self.z - o.z)
    }
}

/// `Hex::neighboringPositions`
pub fn neighboring_positions(h: HexCoords) -> Vec<HexCoords> {
    let (x, y) = (h.x, h.y);
    vec![
        HexCoords::new(x + 1, y),
        HexCoords::new(x + 1, y - 1),
        HexCoords::new(x, y - 1),
        HexCoords::new(x - 1, y),
        HexCoords::new(x - 1, y + 1),
        HexCoords::new(x, y + 1),
    ]
}

/// `Hex::followingVertex`: the next vertex around the hexagon
pub fn following_vertex(hex: HexCoords, v: VertexCoords) -> VertexCoords {
    let (hx, hy, hz) = (hex.x, hex.y, -hex.x - hex.y);
    let (mut dx, mut dy, mut dz) = (v.x - hx, v.y - hy, v.z - hz);
    if dx == 0 && dy == 0 {
        dx = -dz;
        dy = 0;
        dz = 0;
    } else if dx == 0 && dz == 0 {
        dz = -dy;
        dx = 0;
        dy = 0;
    } else if dy == 0 && dz == 0 {
        dy = -dx;
        dx = 0;
        dz = 0;
    }
    VertexCoords::new(hx + dx, hy + dy, hz + dz)
}

/// `Polyomino`: hexagons without holes
#[derive(Debug, Clone, Default)]
pub struct Polyomino {
    /// the hexagons in insertion order
    pub list: Vec<HexCoords>,
    present: HashSet<HexCoords>,
    pub pentagon_vertices: Vec<VertexCoords>,
}

impl Polyomino {
    pub fn new() -> Self {
        Polyomino::default()
    }

    pub fn size(&self) -> usize {
        self.list.len()
    }

    pub fn clear(&mut self) {
        self.list.clear();
        self.present.clear();
    }

    pub fn has_hex(&self, c: HexCoords) -> bool {
        self.present.contains(&c)
    }

    pub fn add_hex(&mut self, c: HexCoords) {
        assert!(self.present.insert(c), "hexagon already present");
        self.list.push(c);
    }

    pub fn remove_hex(&mut self, c: HexCoords) {
        assert!(self.present.remove(&c), "no such hexagon");
        if let Some(i) = self.list.iter().position(|&h| h == c) {
            self.list.remove(i);
        }
    }

    /// `isTheSameAs`: same shape up to translation and rotation (mirror images differ)
    pub fn is_the_same_as(&self, p: &Polyomino) -> bool {
        if self.size() != p.size() {
            return false;
        }
        let mut target: Vec<HexCoords> = p.list.clone();
        if target.is_empty() {
            return true;
        }
        let mut lowestx = self.list[0].x;
        let mut lowesty = self.list[0].y;
        for h in &self.list {
            lowestx = lowestx.min(h.x);
            lowesty = lowesty.min(h.y);
        }
        for _ in 0..6 {
            let (mut lowest_target_x, mut lowest_target_y) = (0, 0);
            for (j, t) in target.iter().enumerate() {
                if j == 0 || t.x < lowest_target_x {
                    lowest_target_x = t.x;
                }
                if j == 0 || t.y < lowest_target_y {
                    lowest_target_y = t.y;
                }
            }
            for t in target.iter_mut() {
                *t = HexCoords::new(
                    t.x + lowestx - lowest_target_x,
                    t.y + lowesty - lowest_target_y,
                );
            }
            if target.iter().all(|&t| self.has_hex(t)) {
                return true;
            }
            for t in target.iter_mut() {
                *t = t.rotate_30_degrees();
            }
        }
        false
    }

    pub fn build_skewed_box_shape(&mut self, x: i32, y: i32, pentagon: bool) {
        self.clear();
        for yy in 0..y {
            for xx in 0..x {
                self.add_hex(HexCoords::new(xx, yy));
            }
        }
        if pentagon {
            self.mark_one_vertex_as_pentagon();
        }
    }

    /// rows alternating between `x` and `x - 1` hexagons
    pub fn build_ragged_smaller_box_shape(&mut self, x: i32, y: i32, pentagon: bool) {
        self.clear();
        let mut startx = 0;
        let mut yy = 0;
        while yy < y {
            for xx in 0..x {
                self.add_hex(HexCoords::new(startx + xx, yy));
            }
            yy += 1;
            if yy >= y {
                break;
            }
            for xx in 0..x - 1 {
                self.add_hex(HexCoords::new(startx + xx, yy));
            }
            startx -= 1;
            yy += 1;
        }
        if pentagon {
            self.mark_one_vertex_as_pentagon();
        }
    }

    /// rows alternating between `x` and `x + 1` hexagons
    pub fn build_ragged_bigger_box_shape(&mut self, x: i32, y: i32, pentagon: bool) {
        self.clear();
        let mut startx = 0;
        let mut yy = 0;
        while yy < y {
            for xx in 0..x {
                self.add_hex(HexCoords::new(startx + xx, yy));
            }
            yy += 1;
            if yy >= y {
                break;
            }
            for xx in 0..x + 1 {
                self.add_hex(HexCoords::new(startx + xx - 1, yy));
            }
            startx -= 1;
            yy += 1;
        }
        if pentagon {
            self.mark_one_vertex_as_pentagon();
        }
    }

    /// zig-zag box with rows of `x` hexagons
    pub fn build_ragged_box_shape(&mut self, x: i32, y: i32, pentagon: bool) {
        self.clear();
        let mut startx = 0;
        let mut yy = 0;
        while yy < y {
            for xx in 0..x {
                self.add_hex(HexCoords::new(startx + xx, yy));
            }
            yy += 1;
            if yy >= y {
                break;
            }
            for xx in 0..x {
                self.add_hex(HexCoords::new(startx + xx, yy));
            }
            startx -= 1;
            yy += 1;
        }
        if pentagon {
            self.mark_one_vertex_as_pentagon();
        }
    }

    /// a round-ish polyomino with the given number of vertices
    pub fn build_with_vertices_n(&mut self, tot_vertices: i32) {
        self.clear();
        self.add_hex(HexCoords::new(0, 0));
        self.add_hex(HexCoords::new(1, 0));
        let mut vertices = 10;
        while vertices < tot_vertices {
            let next_coords = self.all_free_neighbors();
            let mut lowest_distance = -1;
            let mut best_i = 0;
            for (i, &c) in next_coords.iter().enumerate() {
                if self.count_neighbors(c) != 2 {
                    continue;
                }
                let distance = c.distance_from(HexCoords::new(0, 0));
                if lowest_distance == -1 || distance < lowest_distance {
                    lowest_distance = distance;
                    best_i = i;
                }
            }
            assert!(lowest_distance != -1);
            self.add_hex(next_coords[best_i]);
            for (i, &c) in next_coords.iter().enumerate() {
                if i == best_i {
                    continue;
                }
                if self.count_neighbors(c) == 3 {
                    self.add_hex(c);
                }
            }
            vertices += 2;
        }
        if vertices - tot_vertices == 1 {
            self.mark_one_vertex_as_pentagon();
        }
    }

    /// `markOneVertexAsPentagon`: turns a vertex into a pentagon corner to get a path with an odd number of vertices
    pub fn mark_one_vertex_as_pentagon(&mut self) {
        let path = self.get_path();
        let n = path.len();
        for (prev_target, this_target, next_target) in [(2, 1, 2), (1, 2, 1)] {
            let mut previous = self.hexagons_at_vertex(path[n - 1]);
            let mut this = self.hexagons_at_vertex(path[0]);
            for i in 0..n {
                let next_i = if i >= n - 1 { 0 } else { i + 1 };
                let next = self.hexagons_at_vertex(path[next_i]);
                if previous == prev_target && this == this_target && next == next_target {
                    self.pentagon_vertices.push(path[i]);
                    return;
                }
                previous = this;
                this = next;
            }
        }
    }

    pub fn hexagons_at_vertex(&self, v: VertexCoords) -> usize {
        self.vertex_neighbors(v).len()
    }

    /// `vertexNeighbors`: the hexagons that share the vertex
    pub fn vertex_neighbors(&self, v: VertexCoords) -> Vec<HexCoords> {
        let direction = v.x + v.y + v.z;
        if direction != 1 && direction != -1 {
            return Vec::new();
        }
        let mut out = Vec::new();
        for c in [
            HexCoords::new(v.x - direction, v.y),
            HexCoords::new(v.x, v.y - direction),
            HexCoords::new(v.x, v.y),
        ] {
            if self.has_hex(c) {
                out.push(c);
            }
        }
        out
    }

    /// `freeVertexNeighborPositions`
    pub fn free_vertex_neighbor_positions(&self, v: VertexCoords) -> Vec<HexCoords> {
        let direction = v.x + v.y + v.z;
        if direction != 1 && direction != -1 {
            return Vec::new();
        }
        let mut out = Vec::new();
        for c in [
            HexCoords::new(v.x - direction, v.y),
            HexCoords::new(v.x, v.y - direction),
            HexCoords::new(v.x, v.y),
        ] {
            if !self.has_hex(c) {
                out.push(c);
            }
        }
        out
    }

    /// `findOuterVertex`: an outer vertex to start the path around the polyomino
    pub fn find_outer_vertex(&self) -> VertexCoords {
        for h in &self.list {
            let vert = VertexCoords::new(h.x + 1, h.y, -h.x - h.y);
            if self.hexagons_at_vertex(vert) == 1 {
                return vert;
            }
        }
        VertexCoords::new(0, 0, 0)
    }

    pub fn count_neighbors(&self, h: HexCoords) -> i32 {
        neighboring_positions(h)
            .into_iter()
            .filter(|&n| self.has_hex(n))
            .count() as i32
    }

    /// `allFreeNeighbors`: positions next to the polyomino, in the order they are met
    pub fn all_free_neighbors(&self) -> Vec<HexCoords> {
        let mut out = Vec::new();
        let mut visited: HashSet<HexCoords> = HashSet::new();
        for &h in &self.list {
            for n in neighboring_positions(h) {
                if self.has_hex(n) {
                    continue;
                }
                if !visited.insert(n) {
                    continue;
                }
                out.push(n);
            }
        }
        out
    }

    /// `isEquivalentWithout`: removing the hexagon keeps the number of vertices
    pub fn is_equivalent_without(&self, c: HexCoords) -> bool {
        if self.count_neighbors(c) != 3 {
            return false;
        }
        let neighs = neighboring_positions(c);
        for i in 0..neighs.len() {
            let i2 = (i + 6 - 1) % 6;
            let i3 = (i + 6 - 2) % 6;
            if self.has_hex(neighs[i]) && self.has_hex(neighs[i2]) && self.has_hex(neighs[i3]) {
                return true;
            }
        }
        false
    }

    /// `coordinatesOfSubstituent`: where a substituent of the atom at `pos` would be
    pub fn coordinates_of_substituent(&self, pos: VertexCoords) -> VertexCoords {
        let neighbors = self.vertex_neighbors(pos);
        assert!(!neighbors.is_empty());
        assert!(neighbors.len() < 3);
        let mut out = pos;
        if neighbors.len() == 1 {
            let parent = neighbors[0].to_vertex_coords();
            let mut v = pos - parent;
            let d = if v.x + v.y + v.z > 0 { 1 } else { -1 };
            if v.x == 0 {
                v.x -= d;
            }
            if v.y == 0 {
                v.y -= d;
            }
            if v.z == 0 {
                v.z -= d;
            }
            out = parent + v;
        } else if neighbors.len() == 2 {
            let parent1 = neighbors[0].to_vertex_coords();
            let v = pos - parent1;
            let parent2 = neighbors[1].to_vertex_coords();
            out = parent2 - v;
        }
        out
    }

    /// `getPath`: the contour of the polyomino
    pub fn get_path(&self) -> Vec<VertexCoords> {
        let mut out = Vec::new();
        let first_vertex = self.find_outer_vertex();
        let mut current_vertex = first_vertex;
        let neighbors = self.vertex_neighbors(first_vertex);
        assert_eq!(neighbors.len(), 1);
        let mut current_hex = neighbors[0];
        let mut next_vertex = following_vertex(current_hex, current_vertex);
        loop {
            if !self.pentagon_vertices.contains(&current_vertex) {
                out.push(current_vertex);
            }
            current_vertex = next_vertex;
            let neighbors = self.vertex_neighbors(current_vertex);
            assert!(neighbors.len() <= 2);
            if neighbors.len() == 2 {
                current_hex = if neighbors[0] != current_hex {
                    neighbors[0]
                } else {
                    neighbors[1]
                };
            }
            next_vertex = following_vertex(current_hex, current_vertex);
            if current_vertex == first_vertex {
                break;
            }
        }
        out
    }
}
