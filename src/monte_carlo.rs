use std::ops::Add;
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
}

impl MonteCarlo<T: Potential> {
    fn new(beta_epsilon: FloatType, particles: Particles, potential: &impl Potential) -> MonteCarlo {
        MonteCarlo {
            beta_epsilon,
            particles,
            energy: 0.0,
            potential,
        }
    }

    pub fn crystal(beta_epsilon: FloatType, density: FloatType, potential: &impl Potential, cells: Vector3<usize>, crystal_type: &str)-> MonteCarlo {

        let particles = match crystal_type {
            "simple" => Particles::simple_unit_cell(0, density, cells, min_cell_size),
            "bcc" => Particles::body_centered_cubic_cell(0, density, cells, min_cell_size),
            "fcc" => Particles::face_centered_cubic_cell(0, density, cells, min_cell_size),
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
        let i = rng.gen_range(0..self.particles.len());
        let mut new_pos = self.particles.position(i);
        new_pos.x += rng.gen_range(-delta..delta);
        new_pos.y += rng.gen_range(-delta..delta);
        new_pos.z += rng.gen_range(-delta..delta);
        let delta_energy =
        if delta_energy < 0.0 || rng.gen::<FloatType>() < (-self.beta_epsilon * delta_energy).exp() {
            self.particles.set_position(i, new_pos);
        }
    }

    fn _energy(particles: &Particles, potential: fn(FloatType) -> FloatType ) -> FloatType {
        let mut energy = 0.0;
        for i in 0..particles.len() {
            for j in i+1..particles.len() {
                energy += potential(particles.dist(i, j));
            }
        }
        energy
    }

    fn energy(&self) -> FloatType {
        self._energy(&self.particles, self.potential)
    }

    fn particle_energy(&self, i: usize) -> FloatType {
        let neighbors = self.particles.neighbors(i);
        let mut energy = 0.0;
        for j in neighbors {
            if i == j {
                continue;
            }
            energy += self.potential(self.particles.dist(i, j));
        }
        energy
    }

    fn delta_energy(&self, i: usize, new_pos: Vector3<FloatType>) -> FloatType {
        let mut delta_energy = 0.0;

        delta_energy
    }


}
