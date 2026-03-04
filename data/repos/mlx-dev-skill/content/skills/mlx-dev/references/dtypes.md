# MLX Data Types Reference

## Contents

1. [Supported Data Types](#supported-data-types)
2. [Critical Gotchas](#critical-gotchas)
3. [Conversion Patterns](#conversion-patterns)
4. [Decision Tree](#decision-tree)

---

## Supported Data Types

| Type | GPU Support | Notes |
|------|-------------|-------|
| `bool_` | Yes | |
| `uint8` | Yes | |
| `uint16` | Yes | |
| `uint32` | Yes | |
| `uint64` | Yes | |
| `int8` | Yes | |
| `int16` | Yes | |
| `int32` | Yes | Default int type |
| `int64` | Yes | |
| `float16` | Yes | |
| `bfloat16` | Yes | M3+ recommended |
| `float32` | Yes | Default float type |
| `float64` | **CPU ONLY** | GPU operations throw! |
| `complex64` | Partial | No matmul, limited ops |

---

## Critical Gotchas

### float64 Is CPU-Only

```python
a = mx.array([1.0], dtype=mx.float64)
mx.exp(a, stream=mx.gpu)  # RuntimeError!

# Solutions:
mx.exp(a, stream=mx.cpu)           # Force CPU
mx.exp(a.astype(mx.float32))       # Cast to float32
```

### bfloat16 From External Sources Gets Misinterpreted

```python
from ml_dtypes import bfloat16
import numpy as np
x = np.array(1., dtype=bfloat16)
mx.array(x)  # Returns complex64, not bfloat16!

# Correct approach:
mx_arr = mx.array(x.astype(np.float32), dtype=mx.bfloat16)
```

**Why this happens:** NumPy doesn't natively support bfloat16. The `ml_dtypes` package stores bfloat16 as uint16 with custom dtype metadata. When MLX reads this array, it sees the raw uint16 bit pattern and misinterprets it as a different type.

### Integer Statistics Can Overflow

```python
a = mx.array([0, 100, 200, 123], dtype=mx.uint8)
mx.mean(a)                      # Returns 41.75 - WRONG!
mx.mean(a.astype(mx.float32))   # Returns 105.75 - CORRECT
```

### Default Types Differ From NumPy

- MLX defaults: float32, int32
- NumPy defaults: float64, int64 (on 64-bit systems)
- NumPy float64 arrays auto-convert to float32 when creating MLX arrays

### bfloat16 Cannot Convert to NumPy

```python
a = mx.array([1.0], dtype=mx.bfloat16)
np.array(a)  # Error!
np.array(a.astype(mx.float32))  # OK
```

---

## Conversion Patterns

### From NumPy

```python
# Standard types
mx.array(np_arr)

# bfloat16 (indirect)
mx.array(np_arr.astype(np.float32), dtype=mx.bfloat16)

# float64 (avoid or force CPU)
mx.array(np_arr.astype(np.float32))
```

### To NumPy

```python
# Standard types
np.array(mx_arr)

# bfloat16 (indirect)
np.array(mx_arr.astype(mx.float32))
```

### Type Casting

```python
x.astype(mx.float16)
x.astype(mx.bfloat16)
x.astype(mx.float32)
```

---

## Decision Tree

```
Need float64 precision?
  → Use stream=mx.cpu or restructure for float32

Working with bfloat16 from external source?
  → Convert through float32: mx.array(x.astype(np.float32), dtype=mx.bfloat16)

Computing statistics on integers?
  → Cast first: arr.astype(mx.float32)

Need complex numbers?
  → Limited support; avoid matmul with complex

Converting to NumPy?
  → bfloat16 must cast to float32 first
```
