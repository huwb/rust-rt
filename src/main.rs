use std::fs::File;
use std::io::prelude::*;
use std::ops::{Add, Div, Mul, Sub};
use std::path::Path;

// cargo run ; start output.ppm

#[derive(Clone, Copy)]
struct Vector {
    x: f32,
    y: f32,
    z: f32,
}

impl Vector {
    fn new() -> Vector {
        Vector {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        }
    }

    fn from_f32(v: f32) -> Vector {
        Vector { x: v, y: v, z: v }
    }

    fn from_f32s(x: f32, y: f32, z: f32) -> Vector {
        Vector { x, y, z }
    }

    fn length(&self) -> f32 {
        f32::sqrt(self.x * self.x + self.y * self.y + self.z * self.z)
    }

    fn normalize(&mut self) {
        let len = self.length();
        self.x /= len;
        self.y /= len;
        self.z /= len;
    }

    fn saturate(&mut self) {
        self.x = saturate(self.x);
        self.y = saturate(self.y);
        self.z = saturate(self.z);
    }

    fn abs(&mut self) {
        self.x = f32::abs(self.x);
        self.y = f32::abs(self.y);
        self.z = f32::abs(self.z);
    }

    fn dot(&self, other: Vector) -> f32 {
        self.x * other.x + self.y * other.y + self.z * other.z
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

fn saturate(v: f32) -> f32 {
    f32::max(0.0, f32::min(v, 1.0))
}

fn sphere(center: Vector, radius: f32, point: Vector) -> f32 {
    (center - point).length() - radius
}

fn map(x: Vector) -> f32 {
    let mut d = 1000.0;
    d = f32::min(d, sphere(Vector::from_f32s(0.0, 0.0, 0.0), 5.0, x));
    d = f32::min(d, sphere(Vector::from_f32s(0.0, -1005.0, 0.0), 1000.0, x));
    d
}

fn normal(pos: Vector) -> Vector {
    let dx = 0.01;
    let c = map(pos);
    let mut result = Vector::from_f32s(
        map(pos + Vector::from_f32s(0.01, 0.0, 0.0)) - c,
        map(pos + Vector::from_f32s(0.0, 0.01, 0.0)) - c,
        map(pos + Vector::from_f32s(0.0, 0.0, 0.01)) - c,
    );
    result.normalize();
    result
}

fn raymarch(ro: Vector, rd: Vector) -> f32 {
    let mut t = 0.0;
    for _ in 0..32 {
        let pt = ro + t * rd;
        let d = map(pt);
        t += d;
        let pt = ro + t * rd;
        let d = map(pt);
        t += d;
        if d < 0.001 || d > 1000.0 {
            break;
        }
    }

    t
}

fn lighting(pt: Vector) -> Vector {
    let mut n = normal(pt);
    n.normalize();
    let mut l = Vector::from_f32s(-1.0, 1.0, 2.0);
    l.normalize();
    Vector::from_f32(n.dot(l))
}

// u, v are in NDC - [-1,1]
fn trace(u: f32, v: f32) -> Vector {
    let fov = 1.0;
    let ro = Vector::from_f32s(0.0, 0.0, 10.0);
    let mut rd = Vector::from_f32s(fov * u, fov * v, -1.0);
    rd.normalize();

    let t = raymarch(ro, rd);

    if t < 1000. {
        lighting(ro + t * rd)
    } else {
        Vector::new()
    }
}

fn main() {
    let mut image = vec![vec![Vector::new(); 200]; 100];

    let width = image[0].len();
    let height = image.len();
    let inv_width = 1.0 / width as f32;
    let inv_height = 1.0 / height as f32;
    let aspect = inv_height / inv_width;
    let inv_width = aspect * inv_width;
    let uoff = 0.5 * aspect;
    let voff = 0.5;

    for y in 0..height {
        for x in 0..width {
            let u = 2.0 * (x as f32 * inv_width - uoff);
            let v = -2.0 * (y as f32 * inv_height - voff);
            let mut i = trace(u, v);
            i.saturate();
            image[y][x] = i;
        }
    }

    write_output(image, "output.ppm");
}

fn write_output(image: Vec<Vec<Vector>>, filename: &str) {
    let filename = Path::new(filename);
    let mut file = File::create(filename).expect(&format!(
        "Could not open for writing: {}",
        filename.display()
    ));

    writeln!(&mut file, "P3").expect("Could not write to file");

    let width = image[0].len();
    let height = image.len();

    writeln!(&mut file, "{} {}", width, height).expect("Could not write to file");
    writeln!(&mut file, "{}", 255).expect("Could not write to file");

    for y in 0..height {
        for x in 0..width {
            let r = (255.99 * f32::sqrt(image[y][x].x)) as i32;
            let g = (255.99 * f32::sqrt(image[y][x].y)) as i32;
            let b = (255.99 * f32::sqrt(image[y][x].z)) as i32;
            write!(&mut file, "{}\t{}\t{}\t", r, g, b).expect("Could not write to file");
        }
        file.write_all(b"\n").expect("Could not write to file");
    }
}
