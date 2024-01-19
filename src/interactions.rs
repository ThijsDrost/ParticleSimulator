use crate::FloatType;

pub trait Potential{
    fn potential(&self, r: FloatType) -> FloatType;
    fn force(&self, r: FloatType) -> FloatType;
    fn min_cell_size(&self) -> FloatType;
}

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

    fn _potential(r: FloatType) -> FloatType {
        let r2 = r * r;
        let r6 = r2 * r2 * r2;
        let r12 = r6 * r6;
        4.0 * (1.0 / r12 - 1.0 / r6)
    }

    fn _force(r: FloatType) -> FloatType {
        let r2 = r * r;
        let r6 = r2 * r2 * r2;
        let r12 = r6 * r6;
        24.0 * (2.0 / (r12*r) - 1.0 / (r6*r))
    }
}

impl Potential for LennardJones_12_6 {
    fn potential(&self, r: FloatType) -> FloatType {
        if r < self.r_cut {
            LennardJones_12_6::_potential(r) - self.r_cut_value
        }
        else {
            0.0
        }

    }


    fn force(&self, r: FloatType) -> FloatType {
        if r < self.r_cut {
            let r2 = r * r;
            let r6 = r2 * r2 * r2;
            let r12 = r6 * r6;
            24.0 * (2.0 / (r12 * r) - 1.0 / (r6 * r))
        }
        else {
            0.0
        }
    }

    fn min_cell_size(&self) -> FloatType {
        self.r_cut
    }
}

pub static WCA: LennardJones_12_6 = LennardJones_12_6::new(2**(1/6));