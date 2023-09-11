use rayon::prelude::*;

pub fn parallel_fn<T: Clone + Send, F: Fn(&T, &T) -> T + Sync>(values: &[T], f: F) -> T {
    match values.len() {
        0 => panic!("Empty slices"),
        1 => values[0].clone(),
        _ => parallel_fn(
            &values
                .to_vec()
                .into_par_iter()
                .chunks(2)
                .map(|es| match es.len() {
                    1 => es[0].clone(),
                    _ => f(&es[0], &es[1]),
                })
                .collect::<Vec<_>>(),
            f,
        ),
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn correct_parallel_fn() {
        let values = vec![1, 2, 3, 4, 5, 6, 7, 8, 9];
        let result = super::parallel_fn(&values, |a, b| a + b);
        assert_eq!(result, 45);
    }
}
