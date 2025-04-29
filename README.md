# 2025-gtdb-rs226

Build sourmash databases for the new GTDB release.

First, clone this repo and cd in. Then:

1. Install and activate conda env

```
mamba env create -f environment.yml
mamba activate 2025-gtdb-rs226
```

2. Build the zipfiles
```
snakemake build -c 1 -j 1 -n
```

3. Double check the zips:
```
snakemake check -c1 -j1 -n
```

4. Index zips to rocksdb and tar
```
snakemake index -c1 -j1 -n
```

*In each case, `-n` is the dry run. Remove `-n` to actually run the workflow.*
