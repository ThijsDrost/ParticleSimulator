use nalgebra::Vector3;
use crate::FloatType;


pub(crate) struct Particles {
    pub positions: Vec<Vector3<FloatType>>,
    pub box_size: Vector3<FloatType>,
    pub temperature: FloatType,
}

impl Particles {
    fn _start(n: usize, temperature: FloatType, cell_length: FloatType, cells: Vector3<usize>) -> Particles {
        let mut particles = Particles {
            positions: Vec::with_capacity(n),
            box_size: Vector3::new(cell_length*cells.x, cell_length*cells.y, cell_length*cells.z),
            temperature: temperature,
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

    pub fn volume(&self) -> FloatType {
        self.box_size.x * self.box_size.y * self.box_size.z
    }

    pub fn body_centered_cubic_cell(n: usize, temperature: FloatType, density: FloatType,
    cells: Vector3<usize>) -> Particles {
        let cell_length = Particles::density_to_len(2, density);
        let mut particles = Particles::simple_unit_cell(n/2 + n%2, temperature, cell_length, cells);
        Particles::_push_simple(&mut particles.positions, n/2, Vector3::new(cell_length/2.0, cell_length/2.0, cell_length/2.0), cell_length, cells);
        particles
    }

    pub fn simple_unit_cell(n: usize, temperature: FloatType, density: FloatType,
                            cells: Vector3<usize>) -> Particles {
        let cell_length = Particles::density_to_len(1, density);
        let mut particles = Particles::_start(n, temperature, cell_length, cells);
        Particles::_push_simple(&mut particles.positions, n, Vector3::zeros(), cell_length, cells);
        particles
    }

    pub fn face_centered_cubic_cell(n: usize, temperature: FloatType, density: FloatType,
                                    cells: Vector3<usize>) -> Particles {
        let mut n = if n % 4 == 0 { n/4 } else { n/4 + 1 };
        let cell_length = Particles::density_to_len(4, density);
        let mut particles = Particles::simple_unit_cell(n, temperature, cell_length, cells);
        n = if n % 4 >= 2 { n/4 } else { n/4 + 1 };
        Particles::_push_simple(&mut particles.positions, n, Vector3::new(cell_length/2.0, cell_length/2.0, 0.0), cell_length, cells);
        n = if n % 4 >= 3 { n/4 } else { n/4 + 1 };
        Particles::_push_simple(&mut particles.positions, n, Vector3::new(0.0, cell_length/2.0, cell_length/2.0), cell_length, cells);
        n = n/4;
        Particles::_push_simple(&mut particles.positions, n, Vector3::new(cell_length/2.0, 0.0, cell_length/2.0), cell_length, cells);
        particles
    }

    fn density_to_len(n: usize, density: FloatType) -> FloatType {
        let ball_volume = std::f64::consts::PI / 6;
        (n as FloatType * ball_volume / density).powf(1.0 / 3.0)
    }

    fn len(&self) -> usize {
        self.positions.len()
    }

    fn get(&self, index: usize) -> Vector3<FloatType> {
        self.positions[index]
    }
}