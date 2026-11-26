# Example data

These files are sample inputs for testing and demonstrating the algorithms
in this package.

Source: World Athletics 100m all-time lists, retrieved January 2024.

Sample Data Attribution: The sample 100m timing datasets (top 20 performances and 2379 all time bests) contain historical athletic performance records compiled from publicly available statistics published by World Athletics (https://worldathletics.org/). These datasets are provided strictly as static benchmark/example inputs for testing and demonstrating the algorithms in this package. This project is not affiliated with, endorsed by, or sponsored by World Athletics.

Each file is a one-column CSV with a header row. Only the time column is kept. Athlete names and other fields from the original lists are not included.

| File | N | Notes |
|---|---|---|
| `athletics20.csv` | 20 | Prefix of the 2,379-row list. 12 distinct times. |
| `worldathletics-mens-100m-2379.csv` | 2,379 | Full kept list. 58 distinct times. |

The same files are shipped in the Python package (`lex1dx.datasets.load`) and in the R package (`inst/extdata`, `load_dataset()`).
