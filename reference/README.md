# Reference data

`manifest.v1.json` locks the two NCBI assemblies, annotation identity, source URL,
file size, SHA-256, sequence classification, coordinate convention, license link,
and known limitations used by Slice 0.

Large reference files are intentionally excluded from source control. Download and
verify them with:

```powershell
python -m pip install -e .
python -m pichia_safe_harbor.cli fetch strain-b
python -m pichia_safe_harbor.cli fetch strain-c
python -m pichia_safe_harbor.cli validate strain-b
python -m pichia_safe_harbor.cli validate strain-c
```

If a complete NCBI Datasets zip is already available, install it through the same
whole-bundle atomic validation path:

```powershell
python -m pichia_safe_harbor.cli install strain-b --archive path/to/strain-b.zip
python -m pichia_safe_harbor.cli install strain-c --archive path/to/strain-c.zip
```

To use pre-downloaded files, place them under `reference/data/<key>/` using the
`local_name` values in the manifest, then run `validate`. Validation checks both
the exact byte size and SHA-256 before any analysis starts.
