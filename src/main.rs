extern crate rand;
extern crate simple_stopwatch;

pub mod threadpool;

use rand::{Rng, ThreadRng};
use simple_stopwatch::Stopwatch;
use std::fs::File;
use std::io::prelude::*;
use std::ops::{Add, Div, Mul, Sub};
use std::path::Path;
use std::sync::{Arc, Mutex};
use threadpool::ThreadPool;

// cargo run ; start output.ppm

#[derive(Clone, Copy, Debug)]
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

    // fn abs(&mut self) {
    //     self.x = f32::abs(self.x);
    //     self.y = f32::abs(self.y);
    //     self.z = f32::abs(self.z);
    // }

    fn dot(&self, other: Vector) -> f32 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    fn lerp(a: Vector, b: Vector, t: f32) -> Vector {
        (1.0 - t) * a + t * b
    }

    fn random_on_sphere(rng: &mut ThreadRng) -> Vector {
        let mut v = Vector::from_f32s(
            rng.next_f32() - 0.5,
            rng.next_f32() - 0.5,
            rng.next_f32() - 0.5,
        );
        v.normalize();
        v
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

struct Ray {
    ro: Vector,
    rd: Vector,
}

impl Ray {
    fn new(ro: Vector, rd: Vector) -> Ray {
        Ray { ro, rd }
    }

    fn get_point(&self, t: f32) -> Vector {
        self.ro + t * self.rd
    }
}

struct HitRecord {
    t: f32,
    normal: Vector,
}

struct Sphere {
    center: Vector,
    radius: f32,
}

impl Sphere {
    fn new(center: Vector, radius: f32) -> Sphere {
        Sphere { center, radius }
    }

    fn intersect(&self, r: &Ray, tmin: f32, tmax: f32) -> Option<HitRecord> {
        let oc = r.ro - self.center;
        let b = oc.dot(r.rd);
        let c = oc.dot(oc) - self.radius * self.radius;
        // det
        let h = b * b - c;
        if h >= 0.0 {
            let sqrt = f32::sqrt(h);

            let t = -b - sqrt;
            if t >= tmin && t <= tmax {
                let mut normal = r.get_point(t) - self.center;
                normal.normalize();
                return Some(HitRecord { t, normal });
            }

            let t = -b + sqrt;
            if t >= tmin && t <= tmax {
                let mut normal = r.get_point(t) - self.center;
                normal.normalize();
                return Some(HitRecord { t, normal });
            }
        }

        None
    }
}

fn saturate(v: f32) -> f32 {
    f32::max(0.0, f32::min(v, 1.0))
}

fn intersect(r: &Ray, tmin: f32, tmax: f32) -> Option<HitRecord> {
    let mut result = None;

    let scene = vec![
        Sphere::new(Vector::from_f32s(0.0, 0.0, -1.0), 0.5),
        Sphere::new(Vector::from_f32s(0.0, -100.5, -1.0), 100.0),
        Sphere::new(Vector::from_f32s(1.0, 0.0, -1.0), 0.5),
        Sphere::new(Vector::from_f32s(-1.0, 0.0, -1.0), 0.5),
    ];

    let mut tmax = tmax;
    for prim in scene {
        if let Some(hr) = prim.intersect(&r, tmin, tmax) {
            tmax = hr.t;
            result = Some(hr);
        }
    }

    result
}

fn render(r: Ray, bounces_rem: i32, rng: &mut ThreadRng) -> Vector {
    if bounces_rem <= 0 {
        return Vector::new();
    }

    if let Some(hr) = intersect(&r, 0.0001, 1000.0) {
        let ro = r.get_point(hr.t);
        let mut rd = hr.normal + 0.999 * Vector::random_on_sphere(rng);
        rd.normalize();

        0.5 * render(Ray { ro, rd }, bounces_rem - 1, rng)
    } else {
        let t = saturate(r.rd.y / 2.0 + 0.5);
        Vector::lerp(Vector::from_f32(1.0), Vector::from_f32s(0.5, 0.7, 1.0), t)
    }
}

// u, v are in NDC - [-1,1]
fn color(u: f32, v: f32, rng: &mut ThreadRng) -> Vector {
    let fov = 1.0;
    let ro = Vector::from_f32s(0.0, 0.0, 0.0);
    let mut rd = Vector::from_f32s(fov * u, fov * v, -1.0);
    rd.normalize();

    let r = Ray::new(ro, rd);
    render(r, 10, rng)
}

fn main() {
    let width = 800;
    let height = 450;
    let spp = 256;
    let thread_count = 8;

    // create a pool of worker threads
    let tp = ThreadPool::new(thread_count);

    // this creates a vector of references to the same data :(
    // let mut image = vec![Arc::new(Mutex::new(vec![Vector::new(); width])); height];
    let mut image = vec![];
    for _ in 0..height {
        image.push(Arc::new(Mutex::new(vec![Vector::new(); width])));
    }

    let inv_width = 1.0 / width as f32;
    let inv_height = 1.0 / height as f32;
    let aspect = inv_height / inv_width;
    let inv_width = aspect * inv_width;
    let uoff = 0.5 * aspect;
    let voff = 0.5;

    let sw = Stopwatch::start_new();

    for y in 0..height {
        let mut yf = y as f32;

        let mut row = Arc::clone(&&image[y]);

        // send this row off for processing by a worker thread
        tp.dispatch(Box::new(move || {
            let mut row = row.lock().unwrap();
            let mut rng = rand::thread_rng();

            for x in 0..width {
                let mut xf = x as f32;
                let mut col = Vector::new();

                for _s in 0..spp {
                    let u = 2.0 * ((xf + rng.next_f32()) * inv_width - uoff);
                    let v = -2.0 * ((yf + rng.next_f32()) * inv_height - voff);

                    let mut i = color(u, v, &mut rng);

                    col = col + i;
                }

                col = col / spp as f32;
                col.saturate();
                row[x] = col;
            }
        }));
    }

    // this will call the little threadies home
    drop(tp);

    let elapsed_s = sw.s();
    println!("Render time: {}s", elapsed_s);

    let sw = Stopwatch::start_new();
    write_output(image, "output.ppm");
    let elapsed_s = sw.s();
    println!("Write image time: {}s", elapsed_s);
}

fn write_output(image: Vec<Arc<Mutex<Vec<Vector>>>>, filename: &str) {
    let filename = Path::new(filename);
    let mut file = File::create(filename).expect(&format!(
        "Could not open for writing: {}",
        filename.display()
    ));

    writeln!(&mut file, "P3").expect("Could not write to file");

    let width = { image[0].lock().unwrap().len() };
    let height = image.len();

    writeln!(&mut file, "{} {}", width, height).expect("Could not write to file");
    writeln!(&mut file, "{}", 255).expect("Could not write to file");

    for y in 0..height {
        let row = &image[y].lock().unwrap();
        for x in 0..width {
            let r = (255.99 * f32::sqrt(row[x].x)) as i32;
            let g = (255.99 * f32::sqrt(row[x].y)) as i32;
            let b = (255.99 * f32::sqrt(row[x].z)) as i32;
            write!(&mut file, "{}\t{}\t{}\t", r, g, b).expect("Could not write to file");
        }
        file.write_all(b"\n").expect("Could not write to file");
    }
}
