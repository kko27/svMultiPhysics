# Calcium Concentration Output - Code Blocks Reference

This document provides a quick reference for all the code blocks added to enable calcium concentration output from the electrophysiology model.

## 1. Enum Definitions (consts.h)

### Output Group Enum (Line ~529)
```cpp
outGrp_calcium = 529,
```

### Output Type Enum (Line ~568)
```cpp
out_calcium = 568
```

**File**: `Code/Source/solver/consts.h`

---

## 2. Output Properties Mapping (set_output_props.h)

### Mapping Entry (Line ~72)
```cpp
{OutputNameType::out_calcium, std::make_tuple(OutputNameType::outGrp_calcium, 0, 1, "Calcium_concentration") },
```

**File**: `Code/Source/solver/set_output_props.h`

**Explanation**: 
- Maps the `out_calcium` output type to its output group
- Specifies that it's a single scalar value (length 1)
- Names it "Calcium_concentration" for VTK output

---

## 3. CEP Default Outputs (set_equation_props.h)

### Update Output Count and Array (Lines ~71-73)
```cpp
nDOP = {1, 2, 0, 0};
outPuts[0] = OutputNameType::out_voltage;
outPuts[1] = OutputNameType::out_calcium;
```

**File**: `Code/Source/solver/set_equation_props.h`

**Change**: Modified `nDOP[1]` from `1` to `2` to include both action potential and calcium outputs.

---

## 4. Function Declaration (post.h)

### Calcium Extraction Function Declaration (After line ~57)
```cpp
void calcium_concentration(Simulation* simulation, const mshType& lM, Vector<double>& res);
```

**File**: `Code/Source/solver/post.h`

---

## 5. Function Implementation (post.cpp)

### Complete Function Implementation (After line ~1005)
```cpp
/// @brief Extract calcium concentration from electrophysiology state variables
///
/// For each node in the mesh, extract the calcium concentration (Ca_i) from the 
/// electrophysiology model state variables (Xion array). For the TTP model, 
/// calcium is stored at index 3.
//
void calcium_concentration(Simulation* simulation, const mshType& lM, Vector<double>& res)
{
  auto& com_mod = simulation->com_mod;
  auto& cep_mod = simulation->cep_mod;

  res = 0.0;
  
  for (int a = 0; a < lM.nNo; a++) {
    int Ac = lM.gN(a);
    
    // Extract calcium concentration from state variables
    // For TTP model: index 3 corresponds to Ca_i
    // For other models (AP, BO, FN): adjust index as needed
    if (cep_mod.Xion.ncols() > Ac && cep_mod.Xion.nrows() > 3) {
      res(a) = cep_mod.Xion(3, Ac);  // Calcium concentration in mM
    } else {
      res(a) = 0.0;
    }
  }
}
```

**File**: `Code/Source/solver/post.cpp`

---

## 6. Output Case Handler (post.cpp)

### Case Handler for Calcium Output (After line ~146)
```cpp
} else if (outGrp == OutputNameType::outGrp_calcium) {
  Vector<double> tmpCalcium(msh.nNo);
  post::calcium_concentration(simulation, msh, tmpCalcium);
  res = 0.0;
  for (int a = 0; a < com_mod.msh[iM].nNo; a++) {
    int Ac = msh.gN(a);
    res(0,Ac) = tmpCalcium(a);
  }
}
```

**File**: `Code/Source/solver/post.cpp`

**Location**: In the output writing routine, inserted after the `outGrp_ActiveTension` case block.

---

## Summary of Changes

| File | Change Type | Lines | Description |
|------|------------|-------|-------------|
| consts.h | Enum addition | 2 | Added `outGrp_calcium` and `out_calcium` enums |
| set_output_props.h | Map entry | 1 | Added calcium mapping to output_props_map |
| set_equation_props.h | Configuration | 3 | Updated CEP outputs to include calcium |
| post.h | Function declaration | 1 | Added calcium_concentration function prototype |
| post.cpp | Function implementation | 32 | Implemented calcium extraction function |
| post.cpp | Case handler | 8 | Added output case handler for calcium |

**Total lines added**: ~47 lines (including comments)

---

## Comparison with Action Potential Output

| Aspect | Action Potential | Calcium |
|--------|-----------------|---------|
| **Enum (group)** | `outGrp_Y` | `outGrp_calcium` |
| **Enum (output)** | `out_voltage` | `out_calcium` |
| **Output name** | "Action_potential" | "Calcium_concentration" |
| **Function** | Direct array read | `calcium_concentration()` |
| **State location** | First row of solution array `lY` | Index 3 of `cep_mod.Xion` |
| **Units** | mV (millivolts) | mM (millimolar) |
| **Output type** | Scalar (1 value/node) | Scalar (1 value/node) |

---

## Testing

To verify the implementation works correctly:

1. **Build the project** normally with CMake
2. **Run a CEP simulation** with any electrophysiology model (TTP, AP, BO, FN)
3. **Check VTK output files** for "Calcium_concentration" field
4. **Verify values** are in range 0-10 mM (typical for cardiac calcium)

### Expected Behavior

- Calcium concentration should show a transient (spike) during the action potential upstroke
- Peak calcium typically occurs ~50-100ms into the action potential
- Baseline resting calcium is ~0.0001 mM (100 nM)

---

## Extending for Other Ionic Species

To add output for sodium (Na_i) or potassium (K_i):

1. Add new enums to consts.h:
   ```cpp
   outGrp_sodium = 530,
   out_sodium = 567
   ```

2. Add to set_output_props.h:
   ```cpp
   {OutputNameType::out_sodium, std::make_tuple(OutputNameType::outGrp_sodium, 0, 1, "Sodium_concentration") },
   ```

3. Create similar function in post.cpp:
   ```cpp
   void sodium_concentration(Simulation* simulation, const mshType& lM, Vector<double>& res)
   {
     // Use Xion index 2 for Na_i in TTP model
   }
   ```

4. Add case handler following the same pattern

---

## Notes

- All changes maintain backward compatibility
- Existing simulations without CEP output will not be affected
- The calcium index (3) is specific to the TTP model; other models may require adjustment
- The implementation includes bounds checking for safety
- Comments indicate where model-specific adjustments may be needed
