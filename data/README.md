# Third-party data

Small published tables, transcribed for the analysis; cite the original papers when using them.

- `chae2021/chae2021_table2_fits.csv`, `chae2021/chae2021_table3_eNenv.csv` -- Tables 2 and 3 of Chae, K.-H., et al. 2021, ApJ 921, 104 (arXiv:2109.04745), transcribed from the arXiv source.
- `mistele2024_table1_mrt.txt` -- machine-readable Table 1 of Mistele, T., et al. 2024, ApJL 969, L3 (arXiv:2406.09685).

SPARC (Lelli, McGaugh & Schombert 2016, AJ 152, 157) is not included; `python3 tools/fetch_sparc.py` downloads it.
- `heracles2009/leroy2009_table3.csv` -- molecular-gas masses and CO scale lengths, Table 3 (with Table 1 distances) of Leroy, A. K., et al. 2009, AJ 137, 4670 (arXiv:0905.4742), transcribed; used by `notebooks/h2_refit_20261103.py` (v2.0.1).
- DiskMass Survey (Martinsson et al. 2013a,b, A&A 557, A130/A131) is not included; `bash tools/fetch_diskmass.sh` downloads the arXiv sources whose tables `notebooks/dms_fb_20261105.py` reads (v2.0.1).
