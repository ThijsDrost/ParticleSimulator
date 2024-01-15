mod particles;
mod molecular_dynamics;
mod interactions;
mod monte_carlo;

use nalgebra::Vector3;

type FloatType = f32;


fn main() {
    let box_size = 2.0;
    let mut particles = particles::Particles::face_centered_cubic_cell(16, 0.5, Vector3::new(4, 4, 4));
    println!("{:?}", particles.positions);
}
