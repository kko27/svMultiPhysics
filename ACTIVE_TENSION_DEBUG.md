# Debugging Zero Active Tension Output

## Summary
The active tension output reads from `cem.Ya()` array, which is populated during `cep_integ()`. The computation happens in `actv_strs_land()`, but only if certain conditions are met.

## Code Flow for Active Tension
1. **During solve**: `cep_integ()` → `cep_integ_l()` → checks conditions → calls `actv_strs_land()` → stores result in `yl` → copies to `cem.Ya(Ac)`
2. **During output**: `active_tension()` reads from `cem.Ya(Ac)` and writes to VTK

## Checklist of Conditions That Must Be True

### 1. ✅ Active Stress Enabled in XML
**Location**: `solver.xml` line 69
```xml
<EM_Active_stress>true</EM_Active_stress>
```
**Status**: ✅ Confirmed from your XML

### 2. Electromechanics Coupling Enabled
**Location**: `solver.xml` line 68  
```xml
<EM_Coupled>true</EM_Coupled>
```
**Status**: ✅ Confirmed from your XML

### 3. Land Model States Initialized  
**Code**: `cep_ion.cpp:367` and `464`
```cpp
Vector<double> Y_land_node;
if (dmn.cep.cepType == ElectrophysiologyModelType::TTP) {
  Y_land_node = cep_mod.Y_land.col(Ac);
}
```
**Check**: Is `cep.cepType` actually set to TTP?

### 4. Y_land_node Has Non-Zero Size
**Code**: `cep_ion.cpp:785, 812, 840`
```cpp
if (Y_land_node.size() > 0) {
  cep_mod.ttp.actv_strs_land(Y_land_node, X(3), lambda_old, lambda_new, dlambda_dt, cep.dt, yl);
}
```
**Issue**: If `Y_land_node` is empty, active tension is NOT computed!

### 5. Calcium Concentration is Non-Zero
**Code**: `X(3)` is calcium concentration passed to `actv_strs_land`
**Check**: Is the TTP model properly computing calcium transients?

### 6. Fiber Stretch is Computed
**Code**: `lambda_new` and `lambda_old` passed to `actv_strs_land`
**Check**: Are fiber directions defined and is I4f computed?

## Diagnostic Steps

### Step 1: Add Debug Output to Check Conditions
Add this to `cep_ion.cpp` around line 367 (domain-based) and 464 (non-domain):

```cpp
// After extracting Y_land_node
if (Ac == 0) {  // Only print for node 0 to avoid spam
  std::cout << "[DEBUG cep_integ] Node " << Ac << std::endl;
  std::cout << "  cepType == TTP? " << (dmn.cep.cepType == ElectrophysiologyModelType::TTP) << std::endl;
  std::cout << "  Y_land_node.size() = " << Y_land_node.size() << std::endl;
  std::cout << "  cem.aStress = " << cem.aStress << std::endl;
  std::cout << "  cem.cpld = " << cem.cpld << std::endl;
}
```

### Step 2: Check Active Tension Computation
Add this to `cep_ion.cpp` around line 786 (inside FE case, similar for RK4/CN2):

```cpp
if (cem.aStress) {
  if (Y_land_node.size() > 0) {
    double yl_before = yl;
    cep_mod.ttp.actv_strs_land(Y_land_node, X(3), lambda_old, lambda_new, dlambda_dt, cep.dt, yl);
    
    if (Ac == 0) {  // Debug print for node 0
      std::cout << "[DEBUG actv_strs_land] Node " << Ac << std::endl;
      std::cout << "  c_Ca = " << X(3) << " mM" << std::endl;
      std::cout << "  lambda_old = " << lambda_old << std::endl;
      std::cout << "  lambda_new = " << lambda_new << std::endl;
      std::cout << "  dlambda_dt = " << dlambda_dt << std::endl;
      std::cout << "  yl_before = " << yl_before << std::endl;
      std::cout << "  yl_after = " << yl << std::endl;
    }
  } else {
    if (Ac == 0) {
      std::cout << "[DEBUG] Node " << Ac << ": Y_land_node is EMPTY - active tension NOT computed!" << std::endl;
    }
  }
}
```

### Step 3: Check `actv_strs_land` Inputs
Add this at the start of `CepModTtp::actv_strs_land()` in `CepModTtp.cpp`:

```cpp
// Debug first call
static bool first_call = true;
if (first_call) {
  std::cout << "[actv_strs_land] First call:" << std::endl;
  std::cout << "  c_Ca = " << c_Ca << " mM" << std::endl;
  std::cout << "  lambda_old = " << lambda_old << std::endl;
  std::cout << "  lambda_new = " << lambda_new << std::endl;
  std::cout << "  dlambda_dt = " << dlambda_dt << std::endl;
  std::cout << "  dt = " << dt << std::endl;
  std::cout << "  Y_land_node = [" << Y_land_node << "]" << std::endl;
  first_call = false;
}
```

## Most Likely Causes

### Cause 1: Y_land_node is Empty (Most Likely!)
**Symptom**: `Y_land_node.size() == 0`
**Reason**: 
- `cep.cepType` is not set to TTP
- `Y_land` global array was not properly initialized
- Column extraction `cep_mod.Y_land.col(Ac)` is returning empty

**Fix**: Check that `cep.cepType == ElectrophysiologyModelType::TTP` is true

### Cause 2: Calcium is Zero
**Symptom**: `X(3) == 0` (no calcium transient)
**Reason**: TTP model is not solving properly or initial conditions are wrong
**Fix**: Check TTP model parameters and initial conditions

### Cause 3: Lambda is 1.0 Everywhere
**Symptom**: `lambda_old == lambda_new == 1.0` (no deformation)
**Reason**: Structural solver hasn't run yet or fiber directions not defined
**Fix**: Check that FSI/CEM coupling is working

### Cause 4: Active Tension Computed But Not Stored
**Symptom**: `yl` is non-zero but `cem.Ya(Ac)` is still zero
**Reason**: Something wrong with storage in lines 394 or 488 of `cep_ion.cpp`

## Quick Test

Add this simple check right before the active tension output in `post.cpp:998`:

```cpp
if (a == 0) {  // Node 0
  std::cout << "[post::active_tension] Node 0: cem.Ya(0) = " << cem.Ya(0) << std::endl;
  std::cout << "  cem.cpld = " << cem.cpld << std::endl;
}
```

This will immediately tell you if `cem.Ya` has data at output time.

## Expected Values
- `c_Ca` (calcium): Should vary between ~0.0001 mM (resting) to ~0.001 mM (peak)
- `lambda`: Should vary between 0.9 - 1.2 depending on deformation
- `yl` (active tension): Should be 0 at rest, ~10^4-10^6 dyne/cm² during contraction
- After conversion in `actv_strs_land`: `Tact` should be similar magnitude

