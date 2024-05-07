use nalgebra::Vector3;
use rand::Rng;
use crate::FloatType;
use crate::particles::Particles;
use crate::interactions::Potential;

struct MonteCarlo<T: Potential> {
    pub beta_epsilon: FloatType,
    particles: Particles,
    energy: FloatType,
    potential: T,
    accepted: u32,
    rejected: u32,
}

impl<T: Potential> MonteCarlo<T> {
    fn new(beta_epsilon: FloatType, particles: Particles, potential: T) -> MonteCarlo<T> {
        MonteCarlo {
            beta_epsilon,
            energy: MonteCarlo::_energy(&particles, &potential),
            particles,
            potential,
            accepted: 0,
            rejected: 0,
        }
    }

    pub fn crystal(beta_epsilon: FloatType, density: FloatType, potential: T, cells: Vector3<usize>, crystal_type: &str)-> MonteCarlo<T> {
        let particles = match crystal_type {
            "simple" => Particles::simple_unit_cell(0, density, cells, potential.min_cell_size()),
            "bcc" => Particles::body_centered_cubic_cell(0, density, cells, potential.min_cell_size()),
            "fcc" => Particles::face_centered_cubic_cell(0, density, cells, potential.min_cell_size()),
            _ => panic!("Unknown crystal type"),
        };
        let beta_epsilon = beta_epsilon;
        MonteCarlo::new(beta_epsilon, particles, potential)
    }

    pub fn get_particle(&self, i: usize) -> Vector3<FloatType> {
        self.particles.position(i)
    }

    pub fn particles(&self) -> &Particles {
        &self.particles
    }

    pub fn step(&mut self, delta: FloatType) {
        let mut rng = rand::thread_rng();
        let particle_index = rng.gen_range(0..self.particles.len());
        let mut new_position = self.particles.position(particle_index).clone_owned();
        new_position.x += rng.gen_range(-delta..delta);
        new_position.y += rng.gen_range(-delta..delta);
        new_position.z += rng.gen_range(-delta..delta);
        let delta_energy = self.delta_energy(particle_index, new_position);
        if delta_energy < 0.0 || rng.gen::<FloatType>() < (-self.beta_epsilon * delta_energy).exp() {
            self.particles.set_position(particle_index, new_position);
            self.energy += delta_energy;
            self.accepted += 1;
        }
        else {
            self.rejected += 1;
        }
    }

    fn _energy(particles: &Particles, potential: &T) -> FloatType {
        let mut energy = 0.0;
        for i in 0..particles.len() {
            for j in i+1..particles.len() {
                energy += potential.potential(particles.particles_distance(i, j));
            }
        }
        energy
    }

    fn energy(&self) -> FloatType {
        MonteCarlo::<T>::_energy(&self.particles, &self.potential)
    }

    fn particle_index_energy(&self, i: usize) -> FloatType {
        let neighbors = self.particles.particle_close_particles(i);
        neighbors.iter().map(|j| self.potential.potential(self.particles.particles_distance(i, *j))).sum()
    }

    fn particle_loc_energy(&self, loc: Vector3<FloatType>, exclude: usize) -> FloatType {
        self.particles.loc_close_particles(loc)
            .iter()
            .filter(|i| **i != exclude)
            .map(|j| self.potential.potential(self.particles.particle_distance(*j, loc)))
            .sum()
    }

    fn delta_energy(&self, i: usize, new_pos: Vector3<FloatType>) -> FloatType {
        self.particle_loc_energy(new_pos, i) - self.particle_index_energy(i)
    }


}
