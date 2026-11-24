# lex1dx (Python)

Python bindings for lex1dx. Install from this directory (`pip install .`) in a full repository checkout; the extension is built from the C++ core at the repository root, with the same strict floating-point flags, so results are bit-identical to the C++ library.

```python
import lex1dx
r = lex1dx.cluster(x, 64)
s = lex1dx.cluster(x, 64, mode="sig5")
```
