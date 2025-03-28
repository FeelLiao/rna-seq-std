# Changelog

## release/1.0

The first version of RNA-Seq Standard.

### v1.0.1

- Change `get_fq` function to automatically recognize the right extension name and paired ends name.
- Add `*.gz` format support for reference genome. The compressed genome will automatically decompress and press it to `hisat2-build`.
- `sample_pre.py` update to v0.1.1, fixed the extension name and paired ends name issue.