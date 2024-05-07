use nalgebra::Vector3;
use rand::distributions::Uniform;
use crate::FloatType;
use crate::PI;

/// A struct to store the particles in the simulation
pub(crate) struct Particles {
    /// The positions of the particles
    positions: Vec<Vector3<FloatType>>,
    /// The total size of the box
    pub box_size: Vector3<FloatType>,
    /// The cell in which the particle is located
    cell: Vec<usize>,
    /// The number of cells in each direction
    cells_num: Vector3<usize>,
}

impl Particles {
    fn start(n: usize, cell_length: FloatType, cells: Vector3<usize>, min_cell_size: FloatType) -> Particles {
        let box_size = Vector3::new(cell_length*(cells.x as FloatType),
                                    cell_length*(cells.y as FloatType),
                                    cell_length*(cells.z as FloatType));
        if box_size.x < min_cell_size || box_size.y < min_cell_size || box_size.z < min_cell_size {
            panic!("Box size is too small");
        }
        Particles {
            positions: Vec::with_capacity(n),
            box_size,
            cell: Vec::with_capacity(n),
            cells_num: Vector3::new((box_size.x/min_cell_size) as usize,
                                    (box_size.y/min_cell_size) as usize,
                                    (box_size.z/min_cell_size) as usize),
        }
    }

    fn push_simple(positions: &mut Vec<Vector3<FloatType>>, n: usize, start: Vector3<FloatType>,
                   cell_length: FloatType, cells: Vector3<usize>){
        let mut placed_particles = 0;
        'outer: for i_x in 0..cells.x {
            for i_y in 0..cells.y {
                for i_z in 0..cells.z {
                    if placed_particles >= n {
                        break 'outer;
                    }
                    positions.push(Vector3::new(start.x + (i_x as FloatType)*cell_length,
                                                start.y + (i_y as FloatType)*cell_length,
                                                start.z + (i_z as FloatType)*cell_length));
                    placed_particles += 1;
                }
            }
        }
    }

    fn put_cell(&mut self){
        for i in 0..self.particle_num() {
            self.cell[i] = self.cell_num(self.positions[i])
        }
    }

    pub fn particles_in_cell(&self, cell: usize) -> Vec<usize> {
        self.cell.iter()
            .filter(|cell_num| **cell_num == cell)
            .map(|&cell_num| cell_num)
            .collect()
    }

    pub fn volume(&self) -> FloatType {
        self.box_size.x * self.box_size.y * self.box_size.z
    }

    pub fn simple_unit_cell(n: usize, density: FloatType, cells: Vector3<usize>, min_cell_size: FloatType) -> Particles {
        let cell_length = Particles::density_to_len(1, density);
        let mut particles = Particles::start(n, cell_length, cells, min_cell_size);
        Particles::push_simple(&mut particles.positions, n, Vector3::zeros(), cell_length, cells);
        particles
    }

    pub fn body_centered_cubic_cell(n: usize, density: FloatType, cells: Vector3<usize>, min_cell_size: FloatType) -> Particles {
        let cell_length = Particles::density_to_len(2, density);
        let mut particles = Particles::simple_unit_cell(n/2 + n%2, cell_length, cells, min_cell_size);
        Particles::push_simple(&mut particles.positions, n/2, Vector3::new(cell_length/2.0, cell_length/2.0, cell_length/2.0), cell_length, cells);
        particles
    }

    pub fn face_centered_cubic_cell(n: usize, density: FloatType, cells: Vector3<usize>, min_cell_size: FloatType) -> Particles {
        let mut n = if n % 4 == 0 { n/4 } else { n/4 + 1 };
        let cell_length = Particles::density_to_len(4, density);
        let mut particles = Particles::simple_unit_cell(n, cell_length, cells, min_cell_size);
        n = if n % 4 >= 2 { n/4 } else { n/4 + 1 };
        Particles::push_simple(&mut particles.positions, n, Vector3::new(cell_length/2.0, cell_length/2.0, 0.0), cell_length, cells);
        n = if n % 4 >= 3 { n/4 } else { n/4 + 1 };
        Particles::push_simple(&mut particles.positions, n, Vector3::new(0.0, cell_length/2.0, cell_length/2.0), cell_length, cells);
        n = n/4;
        Particles::push_simple(&mut particles.positions, n, Vector3::new(cell_length/2.0, 0.0, cell_length/2.0), cell_length, cells);
        particles
    }

    fn density_to_len(n: usize, density: FloatType) -> FloatType {
        let ball_volume = FloatType::PI() / 6.0;
        (n as FloatType * ball_volume / density).powf(1.0 / 3.0)
    }

    pub fn particle_num(&self) -> usize {
        self.positions.len()
    }

    pub fn position(&self, index: usize) -> Vector3<FloatType> {
        self.positions[index]
    }

    pub fn density(&self) -> FloatType {
        (self.particle_num() as FloatType)*FloatType::PI() / (6.0*self.volume())
    }

    pub fn get_positions(&self) -> &Vec<Vector3<FloatType>> {
        &self.positions
    }

    pub fn particle_distance(&self, i: usize, j: Vector3<FloatType>) -> FloatType {
        self.least_image_distance(self.positions[i] - j)
    }

    pub fn particles_distance(&self, i: usize, j: usize) -> FloatType {
        self.least_image_distance(self.positions[i] - self.positions[j])
    }

    fn least_image_distance(&self, mut vec: Vector3<FloatType>) -> FloatType {
        for i in 0..3 {
            if vec[i] > self.box_size[i]/2.0 {
                vec[i] -= self.box_size[i];
            } else if vec[i] < -self.box_size[i]/2.0 {
                vec[i] += self.box_size[i];
            }
        }
        vec.norm()
    }

    pub fn cell_num(&self, loc: Vector3<FloatType>) -> usize {
        let cell_x = (self.cells_num.x as FloatType * loc.x / self.box_size.x) as usize;
        let cell_y = (self.cells_num.y as FloatType * loc.y / self.box_size.y) as usize;
        let cell_z = (self.cells_num.z as FloatType * loc.z / self.box_size.z) as usize;
        cell_x + cell_y*self.cells_num.x + cell_z*self.cells_num.x*self.cells_num.y
    }

    pub fn cell_coords(&self, cell_index: usize) -> Vector3<usize> {
        Vector3::new(cell_index % self.cells_num.x,
                     (cell_index / self.cells_num.x) % self.cells_num.y,
                     cell_index / (self.cells_num.x * self.cells_num.y))
    }

    pub fn neighbor_cells_particle(&self, particle_index: usize) -> Vec<usize> {
        self.neighbor_cells_cell(self.cell[particle_index])
    }

    /// Includes cell itself
    pub fn neighbor_cells_cell(&self, cell_index: usize) -> Vec<usize> {
        let cell_indexes = self.cell_coords(cell_index);
        let mut neighbors = Vec::new();
        for x_neigh in -1..1 {
            for y_neigh in -1..1 {
                for z_neigh in -1..1 {
                    let x = if (x_neigh != -1 || cell_indexes.x > 0) {(cell_indexes.x as i16 + x_neigh) as usize}
                    else { self.cells_num.x - 1 };
                    let y = if (y_neigh != -1 || cell_indexes.y > 0) {(cell_indexes.y as i16 + y_neigh) as usize}
                    else { self.cells_num.y - 1 };
                    let z = if (z_neigh != -1 || cell_indexes.z > 0) {(cell_indexes.z as i16 + z_neigh) as usize}
                    else { self.cells_num.z - 1 };
                    let cell = x + y*self.cells_num.x + z*self.cells_num.x*self.cells_num.y;
                    neighbors.push(cell);
                }
            }
        }
        neighbors.sort();
        neighbors.dedup();
        neighbors
    }

    pub fn neighbor_cells_loc(&self, loc: Vector3<FloatType>) -> Vec<usize> {
        self.neighbor_cells_cell(self.cell_num(loc))
    }

    /// Excludes particle itself
    pub fn particle_close_particles(&self, particle_index: usize) -> Vec<usize> {
        self.neighbor_cells_particle(particle_index)
            .iter()
            .map(|cell| self.particles_in_cell(*cell))
            .flatten()
            .filter(|&i| i != particle_index)
            .collect()
    }

    pub fn loc_close_particles(&self, loc: Vector3<FloatType>) -> Vec<usize> {
        self.neighbor_cells_loc(loc)
            .iter()
            .map(|cell| self.particles_in_cell(*cell))
            .flatten()
            .collect()
    }

    pub fn cell_close_particles(&self, cell_index: usize) -> Vec<usize> {
        self.neighbor_cells_cell(cell_index)
            .iter()
            .map(|cell| self.particles_in_cell(*cell))
            .flatten()
            .collect()
    }

    pub fn len(&self) -> usize {
        self.positions.len()
    }

    pub fn set_position(&mut self, index: usize, position: Vector3<FloatType>) {
        self.positions[index] = position;
        self.cell[index] = self.cell_num(position);
    }

}