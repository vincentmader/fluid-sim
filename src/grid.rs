
pub fn initialize_grid(
    dim: &Vec<usize>,
) -> Vec<f64> {
    let mut grid = Vec::new();
    let mut nr_of_cells = 1;
    for d in dim.iter() { nr_of_cells *= d; }
    let mut c = 0;
    for _ in 0..nr_of_cells {
        let cell = c as f64;
        // let cell = 0.;
        grid.push(cell);
        c += 1;
    }
    grid
}

pub fn get_idx_from_pos_vec(
    pos: &Vec<usize>,
    dim: &Vec<usize>,
) -> usize {
    if pos.len() != dim.len() { panic!(); }
    if pos.len() <= 0 { panic!(); }
    // let mut res = 0;
    // for d in 0..dim.len() {
    //     res += pos[d] * prod(&dim, d);
    // }
    (0..dim.len()).map(|d| pos[d] * prod(&dim, d)).sum()
}

fn prod(dim: &Vec<usize>, n: usize) -> usize {  // product of dimensions of to nth entry
    let mut res = 1;
    for i in 0..n { res *= dim[i]; }
    res
}

fn floor(x: usize, y: usize) -> usize { 
    ( x - (x % y) ) / y 
}

fn get_pos_vec_from_idx(
    idx: usize,
    dims: &Vec<usize>,
) -> Vec<usize> {

    let mut res = Vec::new();
    for dim in (0..dims.len()).rev() {
        let mut coord = idx;
        for dim2 in 0..res.len() {
            coord -= res[dim2] * prod(dims, dims.len()-dim2-1);
        }
        if dim != 0 { coord = floor(coord, prod(dims, dim)); }
        res.push(coord);
    }
    assert_eq!(dims.len(), res.len());
    res.into_iter().rev().collect()
}

pub fn test() {
    let dim: Vec<usize> = Vec::from([10, 10, 18, 20, 19]);

    let pos = Vec::from([8, 9, 9, 2, 16]);
    let idx = get_idx_from_pos_vec(&pos, &dim);
    let pos_new = get_pos_vec_from_idx(idx, &dim);
    // println!("{:?}", pos);
    // println!("{:?}", pos_new);
    assert_eq!(pos, pos_new);
    println!("\nIt works! (nD-vec -> idx -> nD-vec)");
}

