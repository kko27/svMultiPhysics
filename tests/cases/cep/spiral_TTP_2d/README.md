# Spiral Wave for TTP Model 

Another approach to spiral wave initialization is to 'rig' the simulation domain such that there are 4 sub-domains which are all initialized to different sets of ionic and gating variable states. More details are provided in this reference. A visualization of this simulation case is shown below

<p align="center">
   <img src="./SpiralWave.jpg" width="600">
</p>

Initial conditions for the Ten-tusscher-Panfilov (TTP) EP model to be specified in the solver.xml to override the default variables that are hard-coded inside the solver. 

## Overview

The user-defined initial conditions override system allows users to specify custom initial states and gating variables directly in the XML configuration file. When `Initial_conditions` is specified in the XML, it automatically overrides the default initial conditions that would normally be set based on the myocardial zone (epicardium, endocardium, or midmyocardium).

## How It Works

1. **XML Specification**: Users specify `Initial_conditions` in the domain configuration
2. **Flag Activation**: The system automatically sets a flag indicating user-defined initial conditions are being used
3. **Override Logic**: During initialization, the system uses the user-defined values instead of defaults
4. **Fallback**: If no `Initial_conditions` are specified, the system uses default initialization based on myocardial zone

## XML Configuration

### Full Override Example

```xml
<Domain id="1">
    <Electrophysiology_model> TTP </Electrophysiology_model>
    <Myocardial_zone> epicardium </Myocardial_zone>
    <!-- User-defined initial conditions - these override the default epicardium values -->
    <Initial_conditions>
        <Initial_States>
            <V>-82.0</V>
            <K_i>135.0</K_i>
            <Na_i>9.5</Na_i>
            <Ca_i>0.00015</Ca_i>
            <Ca_ss>0.00015</Ca_ss>
            <Ca_sr>0.18</Ca_sr>
            <R_bar>0.0</R_bar>
        </Initial_States>
        <Gating_variables>
            <Rectifier_current>
                <x_r1>0.008</x_r1>
                <x_r2>0.995</x_r2>
                <x_s>0.012</x_s>
            </Rectifier_current>
            <Fast_sodium_current>
                <m>0.001</m>
                <h>0.995</h>
                <j>0.995</j>
            </Fast_sodium_current>
            <Slow_inward_current>
                <d>0.001</d>
                <f>0.995</f>
                <f2>0.995</f2>
                <fcass>0.995</fcass>
            </Slow_inward_current>
            <Transient_outward_current>
                <s>0.995</s>
                <r>0.001</r>
            </Transient_outward_current>
        </Gating_variables>
    </Initial_conditions>
</Domain>
```

### Partial Override Example

```xml
<Domain id="1">
    <Electrophysiology_model> TTP </Electrophysiology_model>
    <Myocardial_zone> epicardium </Myocardial_zone>
    <!-- Partial override - only specified values override defaults -->
    <Initial_conditions>
        <Initial_States>
            <V>-82.0</V>
            <K_i>135.0</K_i>
            <!-- Other values (Na_i, Ca_i, Ca_ss, Ca_sr, R_bar) use epicardium defaults -->
        </Initial_States>
        <Gating_variables>
            <Rectifier_current>
                <x_r1>0.008</x_r1>
                <!-- x_r2 and x_s use epicardium defaults -->
            </Rectifier_current>
            <!-- Other gating variable groups use epicardium defaults -->
        </Gating_variables>
    </Initial_conditions>
</Domain>
```

### Default Initialization (No Override)

```xml
<Domain id="1">
    <Electrophysiology_model> TTP </Electrophysiology_model>
    <Myocardial_zone> epicardium </Myocardial_zone>
    <!-- No Initial_conditions specified - uses epicardium defaults -->
</Domain>
```

## Available Parameters

### Initial States (7 parameters)
1. `V` - Membrane potential [mV]
2. `K_i` - Intracellular potassium concentration [mM]
3. `Na_i` - Intracellular sodium concentration [mM]
4. `Ca_i` - Intracellular calcium concentration [mM]
5. `Ca_ss` - Subspace calcium concentration [mM]
6. `Ca_sr` - Sarcoplasmic reticulum calcium concentration [mM]
7. `R_bar` - Ryanodine receptor state variable [-]

### Gating Variables (12 parameters)

#### Rectifier Current
- `x_r1` - Rapid delayed rectifier activation gate [-]
- `x_r2` - Rapid delayed rectifier inactivation gate [-]
- `x_s` - Slow delayed rectifier activation gate [-]

#### Fast Sodium Current
- `m` - Fast sodium activation gate [-]
- `h` - Fast sodium fast inactivation gate [-]
- `j` - Fast sodium slow inactivation gate [-]

#### Slow Inward Current
- `d` - L-type calcium activation gate [-]
- `f` - L-type calcium inactivation gate [-]
- `f2` - L-type calcium inactivation gate 2 [-]
- `fcass` - L-type calcium inactivation gate (Ca-dependent) [-]

#### Transient Outward Current
- `s` - Transient outward activation gate [-]
- `r` - Transient outward inactivation gate [-]

## Implementation Details

### Key Changes Made

1. **Added Flag**: `user_initial_conditions_defined` flag in `CepModTtp` class
2. **Modified Initialization**: `init()` method checks the flag and uses user values when set
3. **XML Integration**: `read_cep_domain()` sets the flag when `Initial_conditions` is specified

### Logic Flow

1. **XML Parsing**: `read_cep_domain()` reads initial conditions from XML
2. **Flag Setting**: Sets `user_initial_conditions_defined = true` when initial conditions are found
3. **Value Assignment**: Assigns user values to TTP model member variables
4. **Initialization**: `init()` method uses user values instead of defaults when flag is set

### Code Structure

```cpp
// In read_files.cpp
if (domain_params->initial_conditions.defined()) {
    // Read and set user-defined values
    cep_mod.ttp.set_user_initial_conditions_flag(true);
}

// In CepModTtp.cpp
void CepModTtp::init(...) {
    if (user_initial_conditions_defined) {
        // Use current member variable values (set by XML)
        copyStateToVectors(X, Xg);
    } else {
        // Use default initialization based on myocardial zone
        // ... default logic ...
    }
}
```