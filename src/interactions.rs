use crate::FloatType;

struct LennardJones_12_6 {
    pub r_cut: FloatType,
    r_cut_value: FloatType,
}

impl LennardJones_12_6 {
    fn new(r_cut: FloatType) -> LennardJones_12_6 {
        LennardJones_12_6 {
            r_cut: r_cut,
            r_cut_value: LennardJones_12_6::_potential(r_cut),
        }
    }

    fn force(r: FloatType) -> FloatType {
        let r2 = r * r;
        let r6 = r2 * r2 * r2;
        let r12 = r6 * r6;
        24.0 * (2.0 / (r12*r) - 1.0 / (r6*r))
    }

    fn _potential(r: FloatType) -> FloatType {
        let r2 = r * r;
        let r6 = r2 * r2 * r2;
        let r12 = r6 * r6;
        4.0 * (1.0 / r12 - 1.0 / r6)
    }

    fn potential_cutoff(&self, r: FloatType, r_cut: FloatType) -> FloatType {
        if r > r_cut {
            0.0
        } else {
            LennardJones_12_6::_potential(r) - self.r_cut_value
        }
    }

    fn potential(r: FloatType) -> FloatType {
        LennardJones_12_6::_potential(r)
    }
}