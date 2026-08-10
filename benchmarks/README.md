# Archived benchmark data

`frontera_size_sweep.csv` is a transcription of Table 1 in
`Parallel_Computing_Project.pdf`. All displayed digits are preserved. The
all-pairs rows leave `theta` blank because the opening angle is not used by that
algorithm.

These values are an archival record, not fresh CI benchmarks. The original
Frontera scheduler output and scaling-run logs were not retained, so this
repository does not claim that every committed figure can be regenerated from
raw measurements. New runs emit machine-readable `SUMMARY` lines that
`plots/plot_results.py` can parse.
