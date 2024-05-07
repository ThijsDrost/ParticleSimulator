mod particles;
mod molecular_dynamics;
mod interactions;
mod monte_carlo;

use nalgebra::Vector3;

type FloatType = f32;

trait PI {
    fn PI() -> FloatType;
}

impl PI for f32 {
    fn PI() -> FloatType {
        std::f32::consts::PI as FloatType
    }
}

impl PI for f64 {
    fn PI() -> FloatType {
        std::f64::consts::PI as FloatType
    }
}


fn main() {
    let box_size = 2.0;
    let mut particles = particles::Particles::face_centered_cubic_cell(16, 0.5, Vector3::new(4, 4, 4), 0.0);
    println!("{:?}", particles.get_positions());
}
