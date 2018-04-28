extern crate rand;
extern crate simple_stopwatch;

pub mod threadpool;
pub mod vector;

use rand::Rng;
use simple_stopwatch::Stopwatch;
use std::fs::File;
use std::io::prelude::*;
use std::path::Path;
use std::sync::{Arc, Mutex};
use threadpool::ThreadPool;
use vector::Vector;

// cargo run ; start output.ppm

static MAX_BOUNCES: usize = 10;

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
            let h = h.sqrt();

            let t = -b - h;
            if t >= tmin && t <= tmax {
                let normal = (r.get_point(t) - self.center).normalize();
                return Some(HitRecord { t, normal });
            }

            let t = -b + h;
            if t >= tmin && t <= tmax {
                let normal = (r.get_point(t) - self.center).normalize();
                return Some(HitRecord { t, normal });
            }
        }

        None
    }
}

fn intersect(r: &Ray, tmin: f32, tmax: f32) -> Option<HitRecord> {
    let mut result = None;

    let scene = vec![
        Sphere::new(Vector::new(0.0, 0.0, -1.0), 0.5),
        Sphere::new(Vector::new(0.0, -100.5, -1.0), 100.0),
        Sphere::new(Vector::new(1.0, 0.0, -1.0), 0.5),
        Sphere::new(Vector::new(-1.0, 0.0, -1.0), 0.5),
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

fn render(r: Ray, bounces_rem: usize, rng: &mut Rng) -> Vector {
    if bounces_rem <= 0 {
        return Vector::zero();
    }

    if let Some(hr) = intersect(&r, 0.0001, 1000.0) {
        let ro = r.get_point(hr.t);
        let rd = (hr.normal + 0.999 * Vector::random_on_sphere(rng)).normalize();
        let albedo = 0.5;
        albedo * render(Ray { ro, rd }, bounces_rem - 1, rng)
    } else {
        let t = vector::saturate(r.rd.y / 2.0 + 0.5);
        Vector::lerp(Vector::from_f32(1.0), Vector::new(0.5, 0.7, 1.0), t)
    }
}

// u, v are in NDC - [-1,1]
fn color(u: f32, v: f32, rng: &mut Rng) -> Vector {
    let fov = 1.0;
    let ro = Vector::new(0.0, 0.0, 0.0);
    let rd = Vector::new(fov * u, fov * v, -1.0).normalize();

    let r = Ray::new(ro, rd);
    render(r, MAX_BOUNCES, rng)
}

fn main() {
    let width = 800;
    let height = 450;
    let spp = 256;
    let thread_count = 8;

    println!(
        "{} x {}, {} spp, {} threads",
        width, height, spp, thread_count
    );

    // create a pool of worker threads
    let tp = ThreadPool::new(thread_count);

    // this creates a vector of references to the same data :(
    // let mut image = vec![Arc::new(Mutex::new(vec![Vector::new(); width])); height];
    let mut image = vec![];
    for _ in 0..height {
        image.push(Arc::new(Mutex::new(vec![Vector::zero(); width])));
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
            let mut rng = rand::weak_rng();

            for x in 0..width {
                let mut xf = x as f32;
                let mut col = Vector::zero();

                for _s in 0..spp {
                    let u = 2.0 * ((xf + rng.next_f32()) * inv_width - uoff);
                    let v = -2.0 * ((yf + rng.next_f32()) * inv_height - voff);

                    let mut i = color(u, v, &mut rng);

                    col = col + i;
                }

                row[x] = (col / spp as f32).saturate();
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

    writeln!(&mut file, "P3").unwrap();

    let width = { image[0].lock().unwrap().len() };
    let height = image.len();

    writeln!(&mut file, "{} {}", width, height).unwrap();
    writeln!(&mut file, "{}", 255).unwrap();

    for y in 0..height {
        let row = &image[y].lock().unwrap();
        for x in 0..width {
            let r = (255.99 * row[x].x.sqrt()) as i32;
            let g = (255.99 * row[x].y.sqrt()) as i32;
            let b = (255.99 * row[x].z.sqrt()) as i32;
            write!(&mut file, "{}\t{}\t{}\t", r, g, b).unwrap();
        }
        file.write_all(b"\n").unwrap();
    }
}
