#![allow(non_snake_case)]
#![allow(non_upper_case)]

mod grid;
use grid::get_idx_from_pos_vec as IX;


fn add_source(N: usize, x: &mut Vec<f64>, s: &Vec<f64>, dt: f64) {
    for idx in 0..(N+2)*(N+2) { x[idx] += dt * s[idx] }
}

fn diffuse(N: usize, b: usize, x: &mut Vec<f64>, x0: &Vec<f64>, diff: f64, dt: f64) {
    let a = dt * diff * (N*N) as f64;
    let dim = Vec::from([N, N]);

    for k in 0..20 {
        for i in 1..N {
            for j in 1..N {
                let idx = IX(&Vec::from([i, j]), &dim);
                let jdx = IX(&Vec::from([i-1, j]), &dim);
                let kdx = IX(&Vec::from([i+1, j]), &dim);
                let ldx = IX(&Vec::from([i, j-1]), &dim);
                let mdx = IX(&Vec::from([i, j+1]), &dim);
                x[idx] = x0[idx] + a*( x[jdx] + x[kdx] + x[ldx] + x[mdx] ) / (1.+4.*a);
            }
        }
        set_bnd(N, b, x);
    }

}

fn advect(N: usize, b: usize, d: &mut Vec<f64>, d0: &Vec<f64>, u: &Vec<f64>, v: &Vec<f64>, dt: f64) {
    let dt0 = dt * N as f64;
    let dim = Vec::from([N, N]);
    for i in 1..N {
        for j in 1..N {
            let N = N as f64;
            let idx = IX(&Vec::from([i, j]), &dim);

            let mut x = i as f64 - dt0*u[idx];
            let mut y = j as f64 - dt0*v[idx];
            if x < 0.5 { x = 0.5; }
            if x > N+0.5 { x = N+0.5; }
            if y < 0.5 { y = 0.5; }
            if y > N+0.5 { y = N+0.5; }

            let i0 = x as usize;
            let j0 = y as usize;
            let i1 = i0 + 1;
            let j1 = j0 + 1;

            let s1 = x - i0 as f64;
            let s0 = 1. - s1;
            let t1 = y - j0 as f64;
            let t0 = 1. - t1;

            let jdx = IX(&Vec::from([i0, j0]), &dim);
            let kdx = IX(&Vec::from([i0, j1]), &dim);
            let ldx = IX(&Vec::from([i1, j0]), &dim);
            let mdx = IX(&Vec::from([i1, j1]), &dim);

            d[idx] = s0*(t0*d0[jdx] + t1*d0[kdx]) + s1*(t0*d0[ldx] + t0*d0[mdx]);
        }
    }
    set_bnd(N, b, d);
}

fn dens_step(N: usize, x: &mut Vec<f64>, x0: &Vec<f64>, u: &Vec<f64>, v: &Vec<f64>, diff: f64, dt: f64) {
    add_source(N, x, x0, dt);
    // TODO swap x & x0 ???
    diffuse(N, 0, x, x0, diff, dt);  
    // TODO swap x & x0 ???    (again)
    advect(N, 0, x, x0, u, v, dt);
}

fn set_bnd(N: usize, b: usize, x: &mut Vec<f64>) {
    let dim = Vec::from([N, N]);
    for i in 1..N {
        x[IX(&Vec::from([0, i]), &dim)] = match b {
            1 => -x[IX(&Vec::from([1, i]), &dim)],
            _ => x[IX(&Vec::from([1, i]), &dim)],
        };
        x[IX(&Vec::from([N+1, i]), &dim)] = match b {
            1 => -x[IX(&Vec::from([N, i]), &dim)],
            _ => x[IX(&Vec::from([N, i]), &dim)],
        };
        x[IX(&Vec::from([i, 0]), &dim)] = match b {
            1 => -x[IX(&Vec::from([i, 1]), &dim)],
            _ => x[IX(&Vec::from([i, 1]), &dim)],
        };
        x[IX(&Vec::from([i, N+1]), &dim)] = match b {
            1 => -x[IX(&Vec::from([i, N]), &dim)],
            _ => x[IX(&Vec::from([i, N]), &dim)],
        };
    }
    x[IX(&Vec::from([0, 0]), &dim)] = 0.5*( x[IX(&Vec::from([1,0]), &dim)] + x[IX(&Vec::from([0,1]), &dim)] );
    x[IX(&Vec::from([0, N+1]), &dim)] = 0.5*( x[IX(&Vec::from([1,N+1]), &dim)] + x[IX(&Vec::from([0,N]), &dim)] );
    x[IX(&Vec::from([N+1, 0]), &dim)] = 0.5*( x[IX(&Vec::from([N,0]), &dim)] + x[IX(&Vec::from([N+1,1]), &dim)] );
    x[IX(&Vec::from([N+1, N+1]), &dim)] = 0.5*( x[IX(&Vec::from([N,N+1]), &dim)] + x[IX(&Vec::from([N+1,N]), &dim)] );
}

// TODO u0 and v0 as mutable references?
fn vel_step(N: usize, u: &mut Vec<f64>, v: &mut Vec<f64>, u0: &mut Vec<f64>, v0: &mut Vec<f64>, visc: f64, dt: f64) {
    add_source(N, u, u0, dt);
    add_source(N, v, v0, dt);
    // TODO swap u & u0 ???
    // TODO swap v & v0 ???
    diffuse(N, 1, u, u0, visc, dt);
    diffuse(N, 2, v, v0, visc, dt);
    project(N, u, v, u0, v0);
    // TODO swap u & u0 ???
    // TODO swap v & v0 ???
    advect(N, 1, u, u0, u0, v0, dt);
    advect(N, 2, v, v0, u0, v0, dt);
    project(N, u, v, u0, v0);
}

fn project(N: usize, u: &mut Vec<f64>, v: &mut Vec<f64>, p: &mut Vec<f64>, div: &mut Vec<f64>) {
    let h = 1. / N as f64;
    let dim = Vec::from([N, N]);

    for i in 1..N {
        for j in 1..N {
            let idx = IX(&Vec::from([i, j]), &dim);
            let jdx = IX(&Vec::from([i+1, j]), &dim);
            let kdx = IX(&Vec::from([i-1, j]), &dim);
            let ldx = IX(&Vec::from([i, j+1]), &dim);
            let mdx = IX(&Vec::from([i, j-1]), &dim);
            div[idx] = -0.5*h*( u[jdx]-u[kdx] + v[ldx]-v[mdx] );
            p[idx] = 0.;
        }
    }
    set_bnd(N, 0, div);
    set_bnd(N, 0, p);

    for k in 0..20 {
        for i in 1..N {
            for j in 1..N {
                let idx = IX(&Vec::from([i, j]), &dim);
                let jdx = IX(&Vec::from([i+1, j]), &dim);
                let kdx = IX(&Vec::from([i-1, j]), &dim);
                let ldx = IX(&Vec::from([i, j+1]), &dim);
                let mdx = IX(&Vec::from([i, j-1]), &dim);
                p[idx] = ( div[idx] + p[kdx] + p[jdx] + p[mdx] + p[ldx] ) / 4.;
            }
        }
        set_bnd(N, 0, p);
    }

    for i in 1..N {
        for j in 1..N {
            let idx = IX(&Vec::from([i, j]), &dim);
            let jdx = IX(&Vec::from([i+1, j]), &dim);
            let kdx = IX(&Vec::from([i-1, j]), &dim);
            let ldx = IX(&Vec::from([i, j+1]), &dim);
            let mdx = IX(&Vec::from([i, j-1]), &dim);
            u[idx] -= 0.5*( p[jdx] - p[kdx] ) / h;
            v[idx] -= 0.5*( p[ldx] - p[mdx] ) / h;
        }
    }
    set_bnd(N, 1, u);
    set_bnd(N, 2, v);
}

fn main() {
    grid::test();

    const dt: f64 = 0.01;
    const visc: f64 = 0.01;
    const diff: f64 = 1.;
    // grid properties
    const N: usize = 30;            // nr. of cells per dimension
    const h: f64 = 1. / N as f64;   // grid spacing
    let dim: Vec<usize> = Vec::from([N+2, N+2]);
    // current 
    let mut dens_c = grid::initialize_grid(&dim);
    let mut u_c = grid::initialize_grid(&dim);
    let mut v_c = grid::initialize_grid(&dim);
    // previous
    let mut dens_p = grid::initialize_grid(&dim);
    let mut u_p = grid::initialize_grid(&dim);
    let mut v_p = grid::initialize_grid(&dim);

    // add_source(N, &mut dens_p, &velx_p, dt);  // TODO

    let nr_of_iterations = 1000;
    for iter_idx in 0..nr_of_iterations {
        // get_from_UI(dens_p, u_p, v_p);  // TODO get density/forces from UI
        vel_step(N, &mut u_c, &mut v_c, &mut u_p, &mut v_p, visc, dt);
        dens_step(N, &mut dens_c, &mut dens_p, &mut u_c, &mut v_c, diff, dt);
        // draw_dens(N, dens_c);  // TODO draw density
        println!("{}", iter_idx);
        save_to_file("d", iter_idx, &dens_c);
        save_to_file("u", iter_idx, &u_c);
        save_to_file("v", iter_idx, &v_c);
    }
}

use std::path::PathBuf;

fn save_to_file(quantity: &str, iter_idx: usize, x: &Vec<f64>) {
    let path_to_file = PathBuf::from(&format!("data/{}/{}.txt", quantity, iter_idx));
    // let file = std::fs::File::open(path_to_file).unwrap();
    let mut content = String::from("");
    for i in x.iter() {
        content += &format!("{}\n", i); 
    }
    // println!("{:?}", path_to_file);
    std::fs::write(&path_to_file, content).unwrap();
}

