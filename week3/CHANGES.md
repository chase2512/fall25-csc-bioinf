# Codon Compatibility Changes Applied

## Summary of Professor's Fixes Applied

### 1. TreeError Exception Class
**Changed:** `class TreeError(Exception)` → `class TreeError(Static[Exception])`
- All Codon exceptions must inherit from `Static[Exception]`

### 2. Set Hashing Support
**Added:** Custom `__hash__` method for `set` class using `@extend`
- Codon doesn't have `frozenset`, so we use regular `set` with custom hash implementation
- This allows TreeNode to use `hash(set(self._children))` instead of `hash(self._children)`

### 3. TreeNode._children Type
**Changed:** `Optional[Tuple[TreeNode, ...]]` → `List[TreeNode]`
- Codon tuples must have compile-time known sizes
- Changed initialization:
  - `self._children = None` → `self._children = []` (for leaf nodes)
  - `self._children = tuple(children)` → `self._children = [i for i in children]` (for intermediate nodes)
- Updated `is_leaf()` method: `return self._children is None` → `return len(self._children) == 0`
- Updated `children` property to return `None` if empty list

### 4. Float Types
**Changed:** All `np.float32` → `np.float64`
- Files affected: `upgma_codon.py`, `nj_codon.py`, `test_phylo_codon.py`
- Changed `MAX_FLOAT = float(np.finfo(np.float32).max)` → `MAX_FLOAT = np.finfo(np.float64).max`
- Changed array type conversions: `.astype(np.float32, ...)` → `.astype(np.float64, ...)`

### 5. Distance Initialization
**Already correct:** All `dist` variables initialized as `0.0` instead of `0`
- Ensures float/int compatibility in Codon

### 6. Import Path Fix
**Added:** `sys.path.insert(0, '../code')` in test files
- Allows Codon to find modules in the code directory
- Simplified from `os.path` usage (which Codon doesn't fully support)

## Files Modified

1. **week3/code/tree_codon.py**
   - Exception class fix
   - Added `@extend class set` with custom hash
   - Changed `_children` from `Optional[Tuple]` to `List`
   - Updated all related methods
   - Changed `__hash__` to use `set` instead of `tuple`

2. **week3/code/upgma_codon.py**
   - Changed to `np.float64`
   - Updated `MAX_FLOAT` calculation

3. **week3/code/nj_codon.py**
   - Changed to `np.float64`
   - Updated `MAX_FLOAT` calculation

4. **week3/test/test_phylo_codon.py**
   - Added `sys.path.insert` for imports
   - Changed test data to use `np.float64`

5. **week3/test/evaluate.sh**
   - Added `PYTHONPATH` handling
   - Fixed unbound variable error

## Note About Loading Data
If you need to load `distances.txt` (not currently used in our tests), use:
```python
from python import numpy as pnp
import numpy.pybridge
distances: np.ndarray[int,2] = pnp.loadtxt("tests/sequence/data/distances.txt", dtype=pnp.int64)
tree = upgma(distances)
```

This is due to a bug in Codon's NumPy parser.
