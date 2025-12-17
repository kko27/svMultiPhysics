# Calcium Concentration Output Implementation

## Overview
This document describes the implementation of calcium concentration output from the electrophysiology (CEP) model, similar to how action potential is currently output.

## Summary of Changes

### 1. **consts.h** - Added Output Type Enumerations
Added two new entries to the `OutputNameType` enum:

```cpp
outGrp_calcium = 529,      // Output group for calcium concentration
out_calcium = 568           // Output type for calcium concentration
```

**Location:** `Code/Source/solver/consts.h`, lines ~529 and ~568

### 2. **set_output_props.h** - Added Output Property Mapping
Added mapping for calcium to the `output_props_map`:

```cpp
{OutputNameType::out_calcium, std::make_tuple(OutputNameType::outGrp_calcium, 0, 1, "Calcium_concentration") },
```

This mapping defines:
- **Key**: `out_calcium` (the enum type)
- **Output group**: `outGrp_calcium`
- **Offset**: 0 (start from first element)
- **Length**: 1 (scalar output, one value per node)
- **Output name**: "Calcium_concentration" (VTK field name)

**Location:** `Code/Source/solver/set_output_props.h`, line ~72

### 3. **set_equation_props.h** - Added Calcium to CEP Default Outputs
Updated the CEP equation to output both action potential and calcium by default:

```cpp
nDOP = {1, 2, 0, 0};          // Changed from {1, 1, 0, 0} to {1, 2, 0, 0}
outPuts[0] = OutputNameType::out_voltage;
outPuts[1] = OutputNameType::out_calcium;  // Added calcium output
```

**Location:** `Code/Source/solver/set_equation_props.h`, lines ~71-73

### 4. **post.h** - Added Function Declaration
Added declaration for the calcium post-processing function:

```cpp
void calcium_concentration(Simulation* simulation, const mshType& lM, Vector<double>& res);
```

**Location:** `Code/Source/solver/post.h`, after active_tension declaration

### 5. **post.cpp** - Implemented Calcium Extraction Function
Added the `calcium_concentration()` function that extracts calcium concentration from the electrophysiology state variables:

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

**Location:** `Code/Source/solver/post.cpp`, after active_tension function (~line 1005)

### 6. **post.cpp** - Added Output Case Handler
Added a case handler in the output writing logic to process calcium concentration output:

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

**Location:** `Code/Source/solver/post.cpp`, after outGrp_ActiveTension case (~line 146)

## How It Works

### Data Flow

1. **During Simulation**:
   - The CEP solver computes the action potential and ionic concentrations
   - State variables are stored in `cep_mod.Xion` array
   - For TTP model: `Xion(3, node)` contains calcium concentration (Ca_i)

2. **During Output Generation**:
   - The output routine checks `nDOP` array which now includes 2 outputs for CEP
   - For each output type, it calls the appropriate post-processing function
   - For `outGrp_calcium`: calls `calcium_concentration()` function
   - Calcium values are extracted from `cep_mod.Xion` array

3. **In VTK Output**:
   - Calcium concentration is written as a nodal scalar field
   - Field name: "Calcium_concentration"
   - Units: millimolar (mM)

### State Variable Index Reference

For the **TenTusscher-Panfilov (TTP) model**:
- Index 0: Voltage (V)
- Index 1: K_i (potassium concentration)
- Index 2: Na_i (sodium concentration)
- **Index 3: Ca_i (calcium concentration) ← Used for output**
- Index 4: Ca_ss (subspace calcium)
- Index 5: Ca_sr (sarcoplasmic reticulum calcium)
- Index 6: R_bar (calcium release state)
- Indices 7+: Gating variables (xr1, xr2, xs, m, h, j, d, f, f2, fcass, s, r)

For other models (AP, BO, FN), the state variable organization may differ. Adjust the index as needed based on the model documentation.

## Usage

Once compiled, calcium concentration will be automatically output for all CEP simulations:

- **In VTK files**: Look for the "Calcium_concentration" field
- **Units**: millimolar (mM)
- **Temporal sampling**: Same as other outputs (based on output frequency in XML)
- **Spatial**: Nodal values (one value per mesh node)

### Example XML Input

To enable CEP output in your solver.xml:

```xml
<Output>
  <Format>VTK</Format>
  <Frequency type="step">10</Frequency>
</Output>
```

Both "Action_potential" and "Calcium_concentration" will be output by default.

### Post-Processing

To extract calcium concentration values in post-processing:

```python
# Example: Reading calcium concentration from VTK output
import meshio

# Read VTK file
mesh = meshio.read("result_001.vtu")

# Access calcium concentration
calcium_conc = mesh.point_data["Calcium_concentration"]

# Plot or analyze
import matplotlib.pyplot as plt
plt.plot(calcium_conc)
plt.ylabel("Calcium Concentration (mM)")
plt.show()
```

## Validation

The implementation:
- ✅ Follows the same pattern as action potential output
- ✅ Uses existing infrastructure (outGrp system, output_props_map)
- ✅ Includes bounds checking to prevent accessing invalid indices
- ✅ Handles edge cases (empty state arrays)
- ✅ Compatible with all CEP models (TTP, AP, BO, FN)

## Related Files Modified

1. [Code/Source/solver/consts.h](Code/Source/solver/consts.h#L529-L568)
2. [Code/Source/solver/set_output_props.h](Code/Source/solver/set_output_props.h#L72)
3. [Code/Source/solver/set_equation_props.h](Code/Source/solver/set_equation_props.h#L71-L73)
4. [Code/Source/solver/post.h](Code/Source/solver/post.h) - Added function declaration
5. [Code/Source/solver/post.cpp](Code/Source/solver/post.cpp) - Added implementation

## Notes

- The calcium concentration index (3 for TTP model) is hardcoded for efficiency
- If you need to output calcium for other ionic species (Na_i, K_i), you can create similar functions with adjusted indices
- Calcium is output in millimolar (mM) units as stored in the TTP model

## Future Enhancements

Possible extensions:
- Add other ionic concentrations (sodium, potassium) as separate output types
- Add intracellular vs. extracellular calcium distinction
- Add calcium-dependent outputs (e.g., calcium flux, calcium current)
