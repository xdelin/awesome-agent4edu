# MLX Automatic Differentiation Reference

## Contents

1. [Basic Patterns](#basic-patterns)
2. [Available Transforms](#available-transforms)
3. [Critical Gotchas](#critical-gotchas)
4. [Control Flow](#control-flow)

---

## Basic Patterns

### Simple Gradient

```python
mx.grad(mx.sin)(mx.array(0.0))  # Returns 1.0
```

### Training Pattern with value_and_grad

```python
def loss_fn(model, x, y):
    return nn.losses.cross_entropy(model(x), y, reduction="mean")

loss_and_grad_fn = nn.value_and_grad(model, loss_fn)
loss, grads = loss_and_grad_fn(model, x, y)
```

### Full Training Step

```python
def loss_fn(model, x, y):
    return nn.losses.cross_entropy(model(x), y, reduction="mean")

def train_step(model, optimizer, x, y):
    loss, grads = nn.value_and_grad(model, loss_fn)(model, x, y)
    optimizer.update(model, grads)
    return loss
```

**Note:** Always define `loss_fn` with explicit parameters `(model, x, y)` rather than using closures. This makes the function more reusable and works correctly with `mx.compile`.

---

## Available Transforms

| Transform | Purpose |
|-----------|---------|
| `mx.grad()` | Reverse-mode autodiff |
| `mx.value_and_grad()` | Loss and gradients together |
| `nn.value_and_grad()` | Model-aware version |
| `mx.vjp()` | Vector-Jacobian product |
| `mx.jvp()` | Jacobian-vector product |
| `mx.vmap()` | Vectorization |
| `mx.checkpoint()` | Gradient checkpointing for memory |

### Transforms Compose

```python
mx.grad(mx.vmap(mx.grad(fn)))  # Valid!
```

---

## Critical Gotchas

### Scatter Operations Have No Gradient

```python
def fn(x):
    output = mx.zeros([4])
    output[mx.array([0, 1])] = x  # No VJP implemented!
    return output.sum()

mx.grad(fn)(mx.array([2, 3]))  # ValueError!

# Workaround: use concatenate or where
def fn(x):
    return mx.concatenate([x, mx.zeros([2])]).sum()
```

### Not All Operations Are Differentiable

Operations without gradients:
- Indexed assignment (scatter)
- Some integer operations
- Operations that break graph tracing

### Gradient Checkpointing for Memory

For very deep models, use checkpointing to trade compute for memory:

```python
from functools import partial

@partial(mx.checkpoint, interval=4)
def forward(model, x):
    for layer in model.layers:
        x = layer(x)
    return x
```

---

## Control Flow

### Control Flow Breaks Lazy Graphs

```python
def fn(x):
    y = expensive_op(x)
    if y > 0:  # Forces evaluation, breaks tracing!
        return positive_path(y)
    return negative_path(y)

# Better: branchless with mx.where
def fn(x):
    y = expensive_op(x)
    return mx.where(y > 0, positive_path(y), negative_path(y))
```

### mx.where for Conditional Operations

```python
# Conditional selection
result = mx.where(condition, if_true, if_false)

# Gradient flows through both branches
# but only the selected values contribute
```
