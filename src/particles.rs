use nalgebra::Vector3;
use crate::FloatType;


pub(crate) struct Particles {
    positions: Vec<Vector3<FloatType>>,
    pub box_size: Vector3<FloatType>,
    cell: Vec<usize>,
    cells_num: Vector3<usize>,
}

impl Particles {
    fn _start(n: usize, cell_length: FloatType, cells: Vector3<usize>, min_cell_size: FloatType) -> Particles {
        let box_size = Vector3::new(cell_length*cells.x, cell_length*cells.y, cell_length*cells.z);
        if box_size.x < min_cell_size || box_size.y < min_cell_size || box_size.z < min_cell_size {
            panic!("Box size is too small");
        }
        let mut particles = Particles {
            positions: Vec::with_capacity(n),
            box_size,
            cell: Vec::with_capacity(n),
            cells_num: Vector3::new((box_size.x/min_cell_size) as usize,
                                    (box_size.y/min_cell_size) as usize,
                                    (box_size.z/min_cell_size) as usize),
        };
        particles
    }

    fn _push_simple(positions: &mut Vec<Vector3<FloatType>>, n: usize, start: Vector3<FloatType>,
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

    fn put_cell(&self){
        for i in 0..self.particle_num() {
            self.cell[i] = self.cell(self.positions[i])
        }
    }

    pub fn particles_in_cell(&self, cell: usize) -> Vec<usize> {
        self.positions.iter().enumerate().filter(|(i, _)| self.cell[i] == cell).collect()
    }

    pub fn volume(&self) -> FloatType {
        self.box_size.x * self.box_size.y * self.box_size.z
    }

    pub fn simple_unit_cell(n: usize, density: FloatType, cells: Vector3<usize>, min_cell_size: FloatType) -> Particles {
        let cell_length = Particles::density_to_len(1, density);
        let mut particles = Particles::_start(n, cell_length, cells, min_cell_size);
        Particles::_push_simple(&mut particles.positions, n, Vector3::zeros(), cell_length, cells);
        particles
    }

    pub fn body_centered_cubic_cell(n: usize, density: FloatType, cells: Vector3<usize>, min_cell_size: FloatType) -> Particles {
        let cell_length = Particles::density_to_len(2, density);
        let mut particles = Particles::simple_unit_cell(n/2 + n%2, cell_length, cells, min_cell_size);
        Particles::_push_simple(&mut particles.positions, n/2, Vector3::new(cell_length/2.0, cell_length/2.0, cell_length/2.0), cell_length, cells);
        particles
    }

    pub fn face_centered_cubic_cell(n: usize, density: FloatType, cells: Vector3<usize>, min_cell_size: FloatType) -> Particles {
        let mut n = if n % 4 == 0 { n/4 } else { n/4 + 1 };
        let cell_length = Particles::density_to_len(4, density);
        let mut particles = Particles::simple_unit_cell(n, cell_length, cells, min_cell_size);
        n = if n % 4 >= 2 { n/4 } else { n/4 + 1 };
        Particles::_push_simple(&mut particles.positions, n, Vector3::new(cell_length/2.0, cell_length/2.0, 0.0), cell_length, cells);
        n = if n % 4 >= 3 { n/4 } else { n/4 + 1 };
        Particles::_push_simple(&mut particles.positions, n, Vector3::new(0.0, cell_length/2.0, cell_length/2.0), cell_length, cells);
        n = n/4;
        Particles::_push_simple(&mut particles.positions, n, Vector3::new(cell_length/2.0, 0.0, cell_length/2.0), cell_length, cells);
        particles
    }

    fn density_to_len(n: usize, density: FloatType) -> FloatType {
        let ball_volume = FloatType::consts::PI / 6;
        (n as FloatType * ball_volume / density).powf(1.0 / 3.0)
    }

    pub fn particle_num(&self) -> usize {
        self.positions.len()
    }

    pub fn position(&self, index: usize) -> Vector3<FloatType> {
        self.positions[index]
    }

    pub fn density(&self) -> FloatType {
        (self.particle_num() as FloatType)*FloatType::consts::PI / (6*self.volume())
    }

    pub fn get_positions(&self) -> &Vec<Vector3<FloatType>> {
        &self.positions
    }

    pub fn dist(&self, i: usize, j: usize) -> FloatType {
        let mut dist = self.positions[i] - self.positions[j];
        for i in 0..3 {
            if dist[i] > self.box_size[i]/2.0 {
                dist[i] -= self.box_size[i];
            } else if dist[i] < -self.box_size[i]/2.0 {
                dist[i] += self.box_size[i];
            }
        }
        dist.norm()
    }

    fn cell_num(&self, loc: Vector3<FloatType>) -> usize {
        let cell_x = (self.cells_num.x as FloatType * loc.x / self.box_size.x) as usize;
        let cell_y = (self.cells_num.y as FloatType * loc.y / self.box_size.y) as usize;
        let cell_z = (self.cells_num.z as FloatType * loc.z / self.box_size.z) as usize;
        cell_x + cell_y*self.cells_num.x + cell_z*self.cells_num.x*self.cells_num.y
    }
    
    fn cell_index(&self, i: usize) -> Vector3<usize> {
        Vector3::new(i % self.cells_num.x,
                     (i / self.cells_num.x) % self.cells_num.y,
                     i / (self.cells_num.x * self.cells_num.y))
    }
    
    fn neighbor_cells(&self, i: usize) -> Vec<usize> {
        let cell_indexes = self.cell_index(i);
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
                    neighbors.append(cell)
                }
            }
        }
        neighbors = neighbors.sort().dedup();
    }

    fn neighbor_particle_indexes(&self, i: usize) -> Vec<usize> {
        self.neighbor_cells(i).iter().map(|cell| self.particles_in_cell(*cell)).flatten().collect()
    }
}