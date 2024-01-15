use nalgebra::Vector3;
use crate::FloatType;

fn lennard_jones_12_6(p1: Particle, p2: Particle) -> Vector3<FloatType> {
    let r = p1.position - p2.position;
    let r2 = r.norm_squared();
    let r6 = r2 * r2 * r2;
    let r12 = r6 * r6;
    let force = 24.0 * (2.0 / r12 - 1.0 / r6) / r2;
    force * r
}


struct Particle {
    position: Vector3<FloatType>,
    velocity: Vector3<FloatType>,
    acceleration: Vector3<FloatType>,
    charge: FloatType,
    mass: FloatType,
}

struct ParticleSystem {
    particles: Vec<Particle>,
}


impl Particle {
    fn new() -> Particle {
        Particle {
            position: Vector3::zeros(),
            velocity: Vector3::zeros(),
            acceleration: Vector3::zeros(),
            charge: 0.0,
            mass: 1.0,
        }
    }

    fn step(&mut self, dt: FloatType) {
        self.position += self.velocity * dt;
        self.velocity += self.acceleration * dt;
    }

    fn apply_force(&mut self, force: Vector3<FloatType>) {
        self.acceleration += force / self.mass;
    }

    fn interact(&mut self, other: &mut Particle, force_func: fn(&Particle, &Particle) -> Vector3<FloatType>) {
        let force = force_func(self, other);
        self.apply_force(force);
        other.apply_force(-force);
    }
}