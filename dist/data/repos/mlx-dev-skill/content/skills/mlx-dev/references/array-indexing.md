# MLX Array Indexing Reference

## Contents

1. [Basic Reading](#basic-reading)
2. [Critical Gotchas](#critical-gotchas)
3. [The at[] Syntax](#the-at-syntax)
4. [Gather and Scatter Operations](#gather-and-scatter-operations)

---

## Basic Reading

Standard slicing works as expected:

```python
arr = mx.arange(10)
arr[3]              # array(3, dtype=int32)
arr[-2]             # array(8, dtype=int32)
arr[2:8:2]          # array([2, 4, 6], dtype=int32)

arr = mx.arange(8).reshape(2, 2, 2)
arr[:, :, 0]        # Multi-dimensional slicing
arr[..., 0]         # Ellipsis works
arr[None].shape     # [1, 8] - None inserts new axis
```

---

## Critical Gotchas

### Array-Based Indexing Requires mx.array

```python
a = mx.array([1, 2, 3])
a[[0, 1]]              # ValueError: Cannot index mlx array using the given type
a[mx.array([0, 1])]    # Works: array([1, 2], dtype=int32)
```

### Slice Indices Must Be Python Integers

```python
x = mx.array([18, 47, 56, 57, 58])
i = mx.array(2)
x[i:i+2]               # ValueError: Slice indices must be integers or None
x[i.item():i.item()+2] # Works, but forces evaluation
```

### Slicing Creates Copies, Not Views

This is **opposite to NumPy**:

```python
# NumPy behavior (views)
a = np.array([1, 2, 3])
b = a[:]
b[2] = 0            # a becomes [1, 2, 0] - original modified!

# MLX behavior (copies)
a = mx.array([1, 2, 3])
b = a[:]            # Independent copy
b[2] = 0            # b is [1, 2, 0], a stays [1, 2, 3]
```

**Why MLX differs:** Views require tracking aliasing relationships between arrays, which complicates lazy evaluation and GPU memory management. Copies are simpler and enable more aggressive kernel fusion and optimization.

### Duplicate Index Updates Are Nondeterministic

```python
a = mx.array([1, 2, 3])
a[[0, 0]] = mx.array([4, 5])  # a[0] could be 4 OR 5 - undefined!
```

### Boolean Mask Reads Are NOT Supported

```python
a = mx.array([1, 2, 3])
mask = mx.array([True, False, True])
a[mask]             # NOT SUPPORTED - would return variable-length array
```

**Boolean mask assignment works:**
```python
a[mask] = mx.array([5.0, 6.0])   # OK
a[mask] = 0.0                    # Scalar broadcasts to True positions
```

**Use mx.where() for conditional selection:**
```python
# Instead of: arr[arr > 0]
result = mx.where(arr > 0, arr, mx.zeros_like(arr))
```

### No Bounds Checking

Out-of-bounds indexing is undefined behavior—no error, just garbage:

```python
a = mx.ones([4, 4])
a[0, 100]           # Returns garbage - no IndexError!
```

---

## The at[] Syntax

Regular `x[idx] += y` only applies one update per index. Use `at[]` for accumulation:

```python
a = mx.array([0, 0])
idx = mx.array([0, 1, 0, 1])  # Each index appears twice

# WRONG - each incremented only once
a[idx] += 1                    # array([1, 1])

# CORRECT - properly accumulated
a = a.at[idx].add(1)          # array([2, 2])
```

### Available at[] Operations

| Operation | Syntax | Use Case |
|-----------|--------|----------|
| Accumulating add | `x.at[idx].add(y)` | Histogram binning, gradient accumulation |
| Accumulating subtract | `x.at[idx].subtract(y)` | Sparse updates |
| Accumulating multiply | `x.at[idx].multiply(y)` | Scaling at indices |
| Accumulating divide | `x.at[idx].divide(y)` | Normalization |
| Element-wise max | `x.at[idx].maximum(y)` | Scatter-max pooling |
| Element-wise min | `x.at[idx].minimum(y)` | Scatter-min |

---

## Gather and Scatter Operations

| NumPy/PyTorch | MLX Equivalent |
|---------------|----------------|
| `np.take()` | `mx.take()` |
| `np.take_along_axis()` | `mx.take_along_axis()` |
| `np.put_along_axis()` | `mx.put_along_axis()` |
| `torch.gather()` | `mx.take_along_axis()` |
| `torch.scatter_add_()` | `arr.at[idx].add()` |
| `np.add.at()` | `arr.at[idx].add()` |
| `np.put()` | Use indexed assignment |
| `np.place()` | Use boolean mask assignment |
| `np.nonzero()` | **NOT AVAILABLE** |
| `np.unique()` | **NOT AVAILABLE** |
