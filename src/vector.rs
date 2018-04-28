use rand::Rng;
use std::ops::{Add, Div, Mul, Sub};

#[derive(Clone, Copy, Debug)]
pub struct Vector {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

impl Vector {
    pub fn new(x: f32, y: f32, z: f32) -> Vector {
        Vector { x, y, z }
    }

    pub fn from_f32(v: f32) -> Vector {
        Vector::new(v, v, v)
    }

    pub fn zero() -> Vector {
        Vector::from_f32(0.0)
    }

    fn length(&self) -> f32 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    pub fn normalize(&mut self) -> Vector {
        let len = self.length();
        self.x /= len;
        self.y /= len;
        self.z /= len;
        *self
    }

    pub fn saturate(&mut self) -> Vector {
        self.x = saturate(self.x);
        self.y = saturate(self.y);
        self.z = saturate(self.z);
        *self
    }

    pub fn dot(&self, other: Vector) -> f32 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    pub fn lerp(a: Vector, b: Vector, t: f32) -> Vector {
        (1.0 - t) * a + t * b
    }

    pub fn random_on_sphere(rng: &mut Rng) -> Vector {
        Vector::new(
            rng.next_f32() - 0.5,
            rng.next_f32() - 0.5,
            rng.next_f32() - 0.5,
        ).normalize()
    }
}

impl Add for Vector {
    type Output = Vector;
    fn add(self, other: Vector) -> Vector {
        Vector {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl Sub for Vector {
    type Output = Vector;
    fn sub(self, other: Vector) -> Vector {
        Vector {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl Mul<f32> for Vector {
    type Output = Vector;
    fn mul(self, other: f32) -> Vector {
        Vector {
            x: self.x * other,
            y: self.y * other,
            z: self.z * other,
        }
    }
}

impl Mul<Vector> for f32 {
    type Output = Vector;
    fn mul(self, other: Vector) -> Vector {
        Vector {
            x: self * other.x,
            y: self * other.y,
            z: self * other.z,
        }
    }
}

impl Div<f32> for Vector {
    type Output = Vector;
    fn div(self, other: f32) -> Vector {
        Vector {
            x: self.x / other,
            y: self.y / other,
            z: self.z / other,
        }
    }
}

pub fn saturate(v: f32) -> f32 {
    v.min(1.0).max(0.0)
}
