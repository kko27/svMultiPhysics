# TTP Electrophysiology Model Initial Conditions: Flow Control Documentation

This document describes the complete logic flow of how the Ten-Tusscher-Panfilov (TTP) electrophysiology model initial conditions (states) are read from the XML configuration file, stored in memory, and utilized by the cardiac electrophysiology solver.

## Overview

The flow can be divided into four main phases:
1. **XML Parsing Phase**: Reading and parsing initial conditions from XML
2. **Storage Phase**: Storing parsed values in domain-specific data structures
3. **Initialization Phase**: Setting up the TTP model with initial conditions
4. **Utilization Phase**: Using initial conditions during solver execution

---

## Phase 1: XML Parsing Phase

### Entry Point
- **Location**: `main.cpp` → `read_files()` → `read_files_ns::read_files()`
- **Trigger**: Simulation starts and reads the solver XML file (e.g., `solver.xml`)

### Step 1.1: XML File Loading
- **Function**: `Parameters::read_xml(std::string file_name)`
- **Location**: `Code/Source/solver/Parameters.cpp:150`
- **Action**: 
  - Uses TinyXML2 library to load and parse the XML file
  - Locates root element `<svMultiPhysicsFile>`
  - Extracts general simulation parameters

```cpp
tinyxml2::XMLDocument doc;
auto error = doc.LoadFile(file_name.c_str());
auto root_element = doc.FirstChildElement(FSI_FILE.c_str());
```

### Step 1.2: Equation Parameters Extraction
- **Function**: `Parameters::set_equation_values(tinyxml2::XMLElement* root_element)`
- **Location**: `Code/Source/solver/Parameters.cpp:198`
- **Action**:
  - Iterates through all `<Add_equation>` elements
  - For each equation with `type="CEP"`, creates an `EquationParameters` object
  - Calls `EquationParameters::set_values()` to parse equation-level parameters

### Step 1.3: Domain Parameters Extraction
- **Function**: `EquationParameters::set_values(tinyxml2::XMLElement* eq_elem)`
- **Location**: `Code/Source/solver/Parameters.cpp:2068`
- **Action**:
  - Iterates through child elements of `<Add_equation>`
  - For each `<Domain id="...">` element:
    - Creates a `DomainParameters` object
    - Calls `DomainParameters::set_values()` to parse domain-specific parameters

### Step 1.4: Initial Conditions XML Element Detection
- **Function**: `DomainParameters::set_values(tinyxml2::XMLElement* domain_elem)`
- **Location**: `Code/Source/solver/Parameters.cpp:1673`
- **Action**:
  - Iterates through child elements of `<Domain>`
  - When encountering `<TTP_initial_conditions>` element:
    - Creates/calls `InitialConditionsParameters::set_values()`

### Step 1.5: Initial States Parsing
- **Function**: `InitialConditionsParameters::set_values(tinyxml2::XMLElement* xml_elem)`
- **Location**: `Code/Source/solver/Parameters.cpp:2986`
- **Action**:
  - Checks for `<Initial_states>` child element
  - If found, calls `InitialStatesParameters::set_values()`

### Step 1.6: Individual State Variable Parsing
- **Function**: `InitialStatesParameters::set_values(tinyxml2::XMLElement* xml_elem)`
- **Location**: `Code/Source/solver/Parameters.cpp:3043`
- **Action**:
  - Iterates through child elements of `<Initial_states>`
  - For each element (e.g., `<V>`, `<K_i>`, `<Na_i>`, etc.):
    - Extracts text value using `item->GetText()`
    - Calls `set_parameter_value(name, value)` to store the value
    - Sets `value_set = true` flag
- **Supported Initial State Variables**:
  - `V`: Membrane potential [mV]
  - `K_i`: Intracellular potassium concentration [mM]
  - `Na_i`: Intracellular sodium concentration [mM]
  - `Ca_i`: Intracellular calcium concentration [mM]
  - `Ca_ss`: Subspace calcium concentration [mM]
  - `Ca_sr`: Sarcoplasmic reticulum calcium concentration [mM]
  - `R_bar`: Ryanodine receptor state variable [-]

### Step 1.7: Gating Variables Parsing (Optional)
- **Function**: `GatingVariablesParameters::set_values(tinyxml2::XMLElement* xml_elem)`
- **Location**: `Code/Source/solver/Parameters.cpp:3113`
- **Action**:
  - Similar to Step 1.6, but for gating variables
  - Parses gating variable elements (e.g., `<x_r1_rectifier>`, `<m_fast_Na>`, etc.)

---

## Phase 2: Storage Phase

### Step 2.1: Domain Setup
- **Function**: `read_cep_domain(CepMod* cep_mod, Simulation* simulation, DomainParameters* domain_params, cepModelType& lDmn)`
- **Location**: `Code/Source/solver/read_files.cpp:1058`
- **Action**:
  - Checks if `domain_params->ttp_initial_conditions.defined()` is true
  - If initial conditions are defined in XML, proceeds to extract and store values

### Step 2.2: Initial States Storage
- **Location**: `Code/Source/solver/read_files.cpp:1061-1071`
- **Action**:
  - Extracts values from `InitialStatesParameters` structure
  - Stores each value into `lDmn.cep.ttp_initial_state` structure
  - Sets `lDmn.cep.ttp_user_initial_state = true` if any state variable is set
- **Code Snippet**:
```cpp
if (initial_states_params.V.defined()) { 
    lDmn.cep.ttp_initial_state.V = initial_states_params.V.value(); 
    any_set = true; 
}
// ... (similar for K_i, Na_i, Ca_i, Ca_ss, Ca_sr, R_bar)
lDmn.cep.ttp_user_initial_state = any_set;
```

### Step 2.3: Gating Variables Storage
- **Location**: `Code/Source/solver/read_files.cpp:1074-1094`
- **Action**:
  - Extracts values from `GatingVariablesParameters` structure
  - Maps XML parameter names to internal state structure fields
  - Stores into `lDmn.cep.ttp_initial_state` structure
- **Mapping Examples**:
  - `x_r1_rectifier` → `ttp_initial_state.x_r1`
  - `m_fast_Na` → `ttp_initial_state.m`
  - `d_slow_in` → `ttp_initial_state.d`
  - etc.

### Step 2.4: Data Structure Location
- **Data Structure**: `cepModelType::ttp_initial_state`
- **Type**: `TenTusscherPanfilovState`
- **Location in Memory**: `eq.dmn[iDmn].cep.ttp_initial_state`
- **Flag**: `eq.dmn[iDmn].cep.ttp_user_initial_state` (boolean)
- **Header File**: `Code/Source/solver/CepMod.h:197-198`

---

## Phase 3: Initialization Phase

### Step 3.1: Solver Initialization Entry
- **Function**: `cep_init(Simulation* simulation)`
- **Location**: `Code/Source/solver/cep_ion.cpp:45`
- **Action**:
  - Called during solver initialization
  - Loops through all equations to find CEP (Cardiac Electrophysiology) equations
  - For each domain node, calls `cep_init_l()` to initialize local state vectors

### Step 3.2: Local Initialization
- **Function**: `cep_init_l(CepMod& cep_mod, cepModelType& cep, int nX, int nG, Vector<double>& X, Vector<double>& Xg)`
- **Location**: `Code/Source/solver/cep_ion.cpp:144`
- **Action**:
  - Checks the electrophysiology model type
  - For TTP model (`ElectrophysiologyModelType::TTP`):
    - Checks if `cep.ttp_user_initial_state` flag is true
    - If true: calls `cep.ttp.init()` with pointer to `cep.ttp_initial_state`
    - If false: calls `cep.ttp.init()` with `nullptr` to use default values

### Step 3.3: TTP Model Initialization
- **Function**: `CepModTtp::init(const int imyo, const int nX, const int nG, Vector<double>& X, Vector<double>& Xg, const TenTusscherPanfilovState* user_state)`
- **Location**: `Code/Source/solver/CepModTtp.cpp:549`
- **Action**:
  - **If `user_state != nullptr`** (user-defined initial conditions):
    - Copies user state: `initial_state = *user_state;`
  - **Else** (default initialization):
    - Uses myocardial zone (`imyo`) to select default state:
      - `imyo = 1`: Epicardium defaults
      - `imyo = 2`: Endocardium defaults
      - `imyo = 3`: Midmyocardium defaults
  - Calls `copyStateToVectors(X, Xg)` to populate state vectors

### Step 3.4: State Vector Population
- **Function**: `CepModTtp::copyStateToVectors(Vector<double>& X, Vector<double>& Xg) const`
- **Location**: `Code/Source/solver/CepModTtp.cpp:842`
- **Action**:
  - Copies initial state structure to state variable vectors:
    - `X(0)` = `V` (membrane potential)
    - `X(1)` = `K_i`
    - `X(2)` = `Na_i`
    - `X(3)` = `Ca_i`
    - `X(4)` = `Ca_ss`
    - `X(5)` = `Ca_sr`
    - `X(6)` = `R_bar`
  - Copies gating variables to gating vector:
    - `Xg(0)` = `x_r1`
    - `Xg(1)` = `x_r2`
    - `Xg(2)` = `x_s`
    - `Xg(3)` = `m`
    - `Xg(4)` = `h`
    - `Xg(5)` = `j`
    - `Xg(6)` = `d`
    - `Xg(7)` = `f`
    - `Xg(8)` = `f2`
    - `Xg(9)` = `fcass`
    - `Xg(10)` = `s`
    - `Xg(11)` = `r`

### Step 3.5: Global State Array Assignment
- **Location**: `Code/Source/solver/cep_ion.cpp:129-134`
- **Action**:
  - Copies local state vectors to global state array `cep_mod.Xion`
  - State variables (`X`) go to first `nX` rows
  - Gating variables (`Xg`) go to rows `nX` through `nX+nG-1`
  - Applied to all mesh nodes within the domain

---

## Phase 4: Utilization Phase

### Step 4.1: Time Integration Loop
- **Function**: `cep_integ(Simulation* simulation, const int iEq, const int iDof, const Array<double>& Dg)`
- **Location**: `Code/Source/solver/cep_ion.cpp:175`
- **Action**:
  - Called at each time step during simulation
  - Extracts current state from `cep_mod.Xion` array
  - Separates into state variables (`X`) and gating variables (`Xg`)
  - Calls domain-specific integration function `cep_integ_l()`

### Step 4.2: Domain-Specific Integration
- **Function**: `cep_integ_l(CepMod& cep_mod, cepModelType& cep, int nX, int nG, Vector<double>& X, Vector<double>& Xg, const double t1, double& yl, const double I4f, const double dt)`
- **Location**: `Code/Source/solver/cep_ion.cpp:358`
- **Action**:
  - For TTP model, selects appropriate time integration scheme
  - Calls TTP-specific integration method:
    - Forward Euler: `cep.ttp.integ_fe()`
    - Runge-Kutta 4: `cep.ttp.integ_rk()`
    - Crank-Nicholson: `cep.ttp.integ_cn2()`

### Step 4.3: State Evolution
- **Functions**: `CepModTtp::integ_fe()`, `CepModTtp::integ_rk()`, `CepModTtp::integ_cn2()`
- **Location**: `Code/Source/solver/CepModTtp.cpp`
- **Action**:
  - Updates state variables (`X`) and gating variables (`Xg`) over time
  - Uses initial conditions as starting point for first time step
  - Subsequent time steps use previous step's state as input

### Step 4.4: State Update Back to Global Array
- **Location**: `Code/Source/solver/cep_ion.cpp:307-312` (multi-domain) or `332-338` (single-domain)
- **Action**:
  - After integration, updates `cep_mod.Xion` array with new state values
  - Copies updated `X` and `Xg` vectors back to global state array
  - Ensures all mesh nodes have consistent state information

---

## Flow Diagram Summary

```
XML File (solver.xml)
    │
    ├─► Parameters::read_xml()
    │       │
    │       ├─► set_equation_values()
    │       │       │
    │       │       └─► EquationParameters::set_values()
    │       │               │
    │       │               └─► DomainParameters::set_values()
    │       │                       │
    │       │                       ├─► InitialConditionsParameters::set_values()
    │       │                       │       │
    │       │                       │       ├─► InitialStatesParameters::set_values()
    │       │                       │       │       └─► Stores: V, K_i, Na_i, Ca_i, Ca_ss, Ca_sr, R_bar
    │       │                       │       │
    │       │                       │       └─► GatingVariablesParameters::set_values()
    │       │                       │               └─► Stores: x_r1, x_r2, x_s, m, h, j, d, f, f2, fcass, s, r
    │       │                       │
    │       └─► [Other parameters...]
    │
    ├─► read_files_ns::read_files()
    │       │
    │       └─► read_cep_domain()
    │               │
    │               └─► Stores into: lDmn.cep.ttp_initial_state
    │                       Sets flag: lDmn.cep.ttp_user_initial_state = true
    │
    ├─► cep_init() [During solver initialization]
    │       │
    │       └─► cep_init_l()
    │               │
    │               └─► CepModTtp::init()
    │                       │
    │                       ├─► If ttp_user_initial_state == true:
    │                       │       initial_state = *ttp_initial_state (user values)
    │                       │
    │                       └─► Else:
    │                               initial_state = defaults[imyo] (zone-based defaults)
    │
    │                       └─► copyStateToVectors(X, Xg)
    │                               │
    │                               └─► Populates X[0..6] and Xg[0..11]
    │                                       │
    │                                       └─► Assigned to cep_mod.Xion(:, node)
    │
    └─► cep_integ() [During each time step]
            │
            └─► cep_integ_l()
                    │
                    └─► CepModTtp::integ_fe/rk/cn2()
                            │
                            └─► Updates X and Xg over time
                                    │
                                    └─► Back to cep_mod.Xion(:, node)
```

---

## Key Data Structures

### Parameter Structures (XML → Memory)
- `InitialStatesParameters`: Holds parsed XML values for initial states
- `GatingVariablesParameters`: Holds parsed XML values for gating variables
- `InitialConditionsParameters`: Container for both above structures

### Domain Structures (Per-Domain Storage)
- `cepModelType::ttp_initial_state`: `TenTusscherPanfilovState` structure storing initial conditions
- `cepModelType::ttp_user_initial_state`: Boolean flag indicating user-defined initial conditions
- `cepModelType::ttp`: `CepModTtp` object containing the TTP model instance

### Global State Arrays
- `CepMod::Xion`: Global array storing state and gating variables for all mesh nodes
  - Dimensions: `[nXion, tnNo]` where `nXion = nX + nG`
  - Rows `[0, nX-1]`: State variables (V, K_i, Na_i, Ca_i, Ca_ss, Ca_sr, R_bar)
  - Rows `[nX, nX+nG-1]`: Gating variables (x_r1, x_r2, x_s, m, h, j, d, f, f2, fcass, s, r)

---

## Decision Points

1. **User-Defined vs. Default Initial Conditions**
   - **Check**: `cep.ttp_user_initial_state` flag
   - **If true**: Use `cep.ttp_initial_state` (from XML)
   - **If false**: Use zone-based defaults based on `cep.imyo`

2. **Multi-Domain vs. Single-Domain Handling**
   - **Check**: `com_mod.dmnId.size() != 0`
   - **If multi-domain**: Average initial conditions across overlapping domains
   - **If single-domain**: Use domain 0's initial conditions directly

3. **Myocardial Zone Selection**
   - **Epicardium** (`imyo = 1`): Default epicardium state
   - **Endocardium** (`imyo = 2`): Default endocardium state
   - **Midmyocardium** (`imyo = 3`): Default midmyocardium state

---

## Important Notes

1. **Partial Override**: Users can specify only a subset of initial conditions. Unspecified values will use defaults.

2. **Per-Domain Initialization**: Each domain can have its own initial conditions, allowing for multi-domain setups with different starting states.

3. **State Persistence**: Initial conditions are used only for the first time step. Subsequent time steps use the previous step's state as input.

4. **XML Parameter Names**: The XML uses different naming conventions than internal C++ structures (e.g., `x_r1_rectifier` in XML maps to `x_r1` internally).

5. **Validation**: The parser will throw exceptions if unknown XML elements are encountered, ensuring configuration file correctness.

---

## File Locations Reference

| Phase | File | Key Functions |
|-------|------|---------------|
| XML Parsing | `Code/Source/solver/Parameters.cpp` | `Parameters::read_xml()`, `DomainParameters::set_values()`, `InitialConditionsParameters::set_values()`, `InitialStatesParameters::set_values()` |
| Storage | `Code/Source/solver/read_files.cpp` | `read_cep_domain()` |
| Initialization | `Code/Source/solver/cep_ion.cpp` | `cep_init()`, `cep_init_l()` |
| TTP Model | `Code/Source/solver/CepModTtp.cpp` | `CepModTtp::init()`, `CepModTtp::copyStateToVectors()` |
| Data Structures | `Code/Source/solver/CepMod.h` | `cepModelType`, `TenTusscherPanfilovState` |
| Time Integration | `Code/Source/solver/cep_ion.cpp` | `cep_integ()`, `cep_integ_l()` |

---

## Example XML Structure

```xml
<Domain id="1">
    <Electrophysiology_model> TTP </Electrophysiology_model>
    <Myocardial_zone> epi </Myocardial_zone>
    
    <TTP_initial_conditions>
        <Initial_states>
            <V>18.167</V>
            <K_i>136.897</K_i>
            <Na_i>8.6105</Na_i>
            <Ca_i>1.2592E-4</Ca_i>
            <Ca_ss>3.6181E-4</Ca_ss>
            <Ca_sr>3.6399</Ca_sr>
            <R_bar>0.9078</R_bar>
        </Initial_states>
        
        <Gating_variables>
            <x_r1_rectifier>6.368E-3</x_r1_rectifier>
            <x_r2_rectifier>0.36774</x_r2_rectifier>
            <x_s_rectifier>9.77E-3</x_s_rectifier>
            <m_fast_Na>0.9338</m_fast_Na>
            <h_fast_Na>0.5276</h_fast_Na>
            <j_fast_Na>0.6837</j_fast_Na>
            <d_slow_in>1.9520E-2</d_slow_in>
            <f_slow_in>0.7896</f_slow_in>
            <f2_slow_in>0.9752</f2_slow_in>
            <fcass_slow_in>0.9954</fcass_slow_in>
            <s_out>0.9941</s_out>
            <r_out>4.2757E-5</r_out>
        </Gating_variables>
    </TTP_initial_conditions>
</Domain>
```

---

*Document generated based on codebase analysis of svMultiPhysics solver*

