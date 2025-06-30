# TTP Model User-Defined Initial Conditions Override

This document explains how to use the user-defined initial conditions override system for the Ten-tusscher-Panfilov (TTP) cardiac electrophysiology model in svMultiPhysics.

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

## Default Values by Myocardial Zone

### Epicardium (imyo = 1)
- `V`: -85.23 mV
- `K_i`: 138.3 mM
- `Na_i`: 10.0 mM
- `Ca_i`: 0.0002 mM
- `Ca_ss`: 0.0002 mM
- `Ca_sr`: 0.2 mM
- `R_bar`: 0.0
- Gating variables: Specific values for epicardium

### Endocardium (imyo = 2)
- `V`: -85.23 mV
- `K_i`: 138.3 mM
- `Na_i`: 10.0 mM
- `Ca_i`: 0.0002 mM
- `Ca_ss`: 0.0002 mM
- `Ca_sr`: 0.2 mM
- `R_bar`: 0.0
- Gating variables: Specific values for endocardium

### Mid-myocardium (imyo = 3)
- `V`: -85.23 mV
- `K_i`: 138.3 mM
- `Na_i`: 10.0 mM
- `Ca_i`: 0.0002 mM
- `Ca_ss`: 0.0002 mM
- `Ca_sr`: 0.2 mM
- `R_bar`: 0.0
- Gating variables: Specific values for mid-myocardium

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

## Advantages

1. **Simple Integration**: Works seamlessly with existing XML structure
2. **Flexible Override**: Can override all or just some initial conditions
3. **No External Files**: All configuration is contained in the XML file
4. **Backward Compatible**: Existing simulations without initial conditions work unchanged
5. **Clear Logic**: Easy to understand when overrides are active

## Usage Scenarios

### Scenario 1: Custom Initial Conditions
```xml
<!-- Use custom initial conditions for research
<Initial_conditions>
    <Initial_States>
        <V>-80.0</V>
        <K_i>140.0</K_i>
        <!-- ... other custom values ... -->
    </Initial_States>
</Initial_conditions>
```

### Scenario 2: Continuation from Previous State
```xml
<!-- Continue simulation from specific state
<Initial_conditions>
    <Initial_States>
        <V>-45.2</V>
        <K_i>135.8</K_i>
        <!-- ... values from previous simulation end ... -->
    </Initial_States>
</Initial_conditions>
```

### Scenario 3: Parameter Study
```xml
<!-- Study effect of specific initial conditions
<Initial_conditions>
    <Initial_States>
        <V>-85.0</V>
        <Ca_i>0.0003</Ca_i>
        <!-- Vary Ca_i while keeping other values constant -->
    </Initial_States>
</Initial_conditions>
```

## Best Practices

1. **Document Values**: Keep track of what initial conditions you're using
2. **Validate Physiologically**: Ensure values are within reasonable ranges
3. **Test Thoroughly**: Verify that your custom initial conditions work correctly
4. **Use Meaningful Names**: Name your XML files to reflect the initial conditions used
5. **Version Control**: Track changes to initial conditions in version control

## Troubleshooting

### Common Issues

1. **Simulation Instability**: Check that initial conditions are physiologically reasonable
2. **Unexpected Behavior**: Verify that the correct myocardial zone is specified
3. **Parameter Conflicts**: Ensure no conflicting parameters are set elsewhere

### Debug Tips

1. **Print Values**: Add debug output to verify loaded initial conditions
2. **Compare with Defaults**: Check against known good default values
3. **Start Simple**: Test with values close to defaults first
4. **Check XML Syntax**: Ensure XML is properly formatted 