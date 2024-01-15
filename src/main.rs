mod monte_carlo;
mod molecular_dynamics;

use nalgebra::Vector3;

type FloatType = f32;


fn main() {
    let box_size = 2.0;
    let mut particles = monte_carlo::Particles::face_centered_cubic_cell(16, 1.0, 2.0, Vector3::new(4, 4, 4));
    println!("{:?}", particles.positions);
}
