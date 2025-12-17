# Architecture & Data Flow Diagrams

## System Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                    SVMultiPhysics Solver                        │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  ┌──────────────────┐         ┌──────────────────┐             │
│  │  CEP Solver      │         │  Other Physics   │             │
│  │  (cep_integ)     │         │  (fluid, struct) │             │
│  └────────┬─────────┘         └──────────────────┘             │
│           │                                                     │
│           ├─→ Compute V_m, gating variables                    │
│           ├─→ Compute ion currents                             │
│           ├─→ Compute ionic concentrations:                    │
│           │   - Xion(0, node) = V (voltage)                   │
│           │   - Xion(1, node) = K_i (potassium)               │
│           │   - Xion(2, node) = Na_i (sodium)                 │
│           │   - Xion(3, node) = Ca_i (CALCIUM) ← We use this  │
│           │   - Xion(4, node) = Ca_ss (subspace Ca)           │
│           │   - Xion(5, node) = Ca_sr (SR calcium)            │
│           │   - Xion(6+, node) = gating variables             │
│           │                                                     │
│           └─→ Store in cep_mod.Xion array                      │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
                            │
                            ↓
┌─────────────────────────────────────────────────────────────────┐
│               Output Configuration                              │
│                                                                 │
│  Read from: set_equation_props.h                               │
│  ┌──────────────────────────────────┐                          │
│  │ CEP Equation Output Configuration │                          │
│  │                                  │                          │
│  │  nDOP = {1, 2, 0, 0}             │                          │
│  │  outPuts[0] = out_voltage        │ → Action_potential      │
│  │  outPuts[1] = out_calcium ← NEW! │ → Calcium_concentration │
│  │                                  │                          │
│  └──────────────────────────────────┘                          │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
                            │
                            ↓
┌─────────────────────────────────────────────────────────────────┐
│           Output Properties Lookup                              │
│                                                                 │
│  From: output_props_map (set_output_props.h)                   │
│                                                                 │
│  out_calcium → outGrp_calcium (OUTPUT GROUP)                   │
│               → offset: 0, length: 1                           │
│               → field name: "Calcium_concentration"             │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
                            │
                            ↓
┌─────────────────────────────────────────────────────────────────┐
│         Post-Processing / Output Generation                     │
│                                                                 │
│  Match: outGrp == outGrp_calcium                               │
│  ┌────────────────────────────────────────┐                    │
│  │ Call: post::calcium_concentration()    │                    │
│  │                                        │                    │
│  │ For each node in mesh:                 │                    │
│  │   calcium[node] = Xion(3, Ac)          │                    │
│  │   (where Ac = global node ID)          │                    │
│  │                                        │                    │
│  └────────────────────────────────────────┘                    │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
                            │
                            ↓
┌─────────────────────────────────────────────────────────────────┐
│                   VTK Output File                               │
│                                                                 │
│  result_001.vtu:                                               │
│  ├─ Points (node coordinates)                                  │
│  ├─ Cells (mesh topology)                                      │
│  ├─ PointData:                                                │
│  │  ├─ Action_potential (float, 1 component)                   │
│  │  ├─ Calcium_concentration (float, 1 component) ← NEW!      │
│  │  ├─ Displacement (if struct)                               │
│  │  └─ ... other fields                                       │
│  └─ CellData: ... (if any)                                    │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

## Execution Flow Diagram

```
START CEP SIMULATION
│
├─→ Initialize state variables
│   └─ Create Xion array with all ion states
│
├─→ FOR each time step:
│   │
│   ├─→ cep_integ()
│   │   └─→ cep_integ_l() for each node
│   │       ├─ Solve ion channel equations
│   │       ├─ Update Xion(3, node) = Ca_i
│   │       └─ Store results
│   │
│   ├─→ [If output_requested at this time step]
│   │   │
│   │   ├─→ FOR each equation:
│   │   │   │
│   │   │   ├─→ FOR each output type (outPuts[i]):
│   │   │   │   │
│   │   │   │   ├─→ Check if outPuts[i] == out_voltage
│   │   │   │   │   └─ Output action potential directly
│   │   │   │   │
│   │   │   │   ├─→ Check if outPuts[i] == out_calcium
│   │   │   │   │   │
│   │   │   │   │   └─→ Call post::calcium_concentration()
│   │   │   │   │       │
│   │   │   │   │       ├─→ FOR each node:
│   │   │   │   │       │   res(node) = Xion(3, global_id)
│   │   │   │   │       │
│   │   │   │   │       └─→ Return calcium array
│   │   │   │   │
│   │   │   │   └─→ Write result to VTK output
│   │   │   │       ├─ Field name: Calcium_concentration
│   │   │   │       └─ Format: Float32, 1 component/node
│   │   │   │
│   │   │   └─→ Write other outputs (stress, displacement, etc.)
│   │   │
│   │   └─→ Finalize VTK file
│   │       result_NNN.vtu ← Contains both Action_potential
│   │                           and Calcium_concentration
│   │
│   └─→ Next time step
│
└─→ END

VTK FILE CONTAINS:
• Action_potential field (mV)
• Calcium_concentration field (mM) ← NEW!
• Other physics outputs as configured
```

## State Variable Organization

```
╔════════════════════════════════════════════════════════╗
║               cep_mod.Xion Array                       ║
║  (nStates × nNodes matrix for TTP model)              ║
╠════════════════════════════════════════════════════════╣
║ Index │ Variable  │ TTP Description                   ║
╠═══════╪═══════════╪═══════════════════════════════════╣
║   0   │ V         │ Membrane voltage (mV)              ║
║   1   │ K_i       │ Intracellular K⁺ (mM)             ║
║   2   │ Na_i      │ Intracellular Na⁺ (mM)            ║
║   3   │ Ca_i      │ Intracellular Ca²⁺ (mM) ✓ USED    ║
║   4   │ Ca_ss     │ Subspace Ca²⁺ (mM)                ║
║   5   │ Ca_sr     │ SR Ca²⁺ (mM)                      ║
║   6   │ R_bar     │ Calcium release channel state     ║
║   7   │ xr1       │ Rapid K current gating (voltage)  ║
║   8   │ xr2       │ Rapid K current gating (time)     ║
║   9   │ xs        │ Slow K current gating             ║
║  10   │ m         │ Na current activation             ║
║  11   │ h         │ Na current inactivation (fast)    ║
║  12   │ j         │ Na current inactivation (slow)    ║
║  13   │ d         │ L-type Ca current activation      ║
║  14   │ f         │ L-type Ca current inactivation    ║
║  15   │ f2        │ L-type Ca current inactivation 2  ║
║  16   │ fcass     │ Ca-dependent inactivation         ║
║  17   │ s         │ Transient outward K curr gating   ║
║  18   │ r         │ Transient outward K curr gating   ║
╚═══════╧═══════════╧═══════════════════════════════════╝

Memory Layout:
┌─────┬─────┬─────┬─────┬─────┬─────┬─────┬─────┬─────┐
│Ca_i │Ca_i │Ca_i │Ca_i │Ca_i │Ca_i │Ca_i │Ca_i │Ca_i│ ← Xion row 3
│ n1  │ n2  │ n3  │ n4  │ n5  │ n6  │ n7  │ n8  │ n9 │
└─────┴─────┴─────┴─────┴─────┴─────┴─────┴─────┴─────┘
 node  node  node  node  node  node  node  node  node
   1     2     3     4     5     6     7     8     9

Access pattern in calcium_concentration():
for (int a = 0; a < lM.nNo; a++) {
  int Ac = lM.gN(a);           // Global node ID
  res(a) = cep_mod.Xion(3, Ac);  // Get Ca_i from row 3
}
```

## Comparison: Action Potential vs Calcium

```
╔══════════════════════════════════════════════════════════╗
║          Data Flow Comparison                            ║
╠══════════════════════════════════════════════════════════╣
║                                                          ║
║  ACTION POTENTIAL (out_voltage)                         ║
║  ══════════════════════════════════════════════════════ ║
║                                                          ║
║  Step 1: Computation in CEP solver                       ║
║          └─→ Xion, gating variables computed            ║
║                                                          ║
║  Step 2: Output request (outGrp_Y)                      ║
║          └─→ Direct array read: lY(0, Ac)               ║
║          └─→ No post-processing function needed         ║
║                                                          ║
║  Step 3: VTK output                                     ║
║          └─→ Action_potential field                     ║
║                                                          ║
╠══════════════════════════════════════════════════════════╣
║                                                          ║
║  CALCIUM CONCENTRATION (out_calcium) ← NEW              ║
║  ══════════════════════════════════════════════════════ ║
║                                                          ║
║  Step 1: Computation in CEP solver                       ║
║          └─→ Xion(3, Ac) = Ca_i computed               ║
║                                                          ║
║  Step 2: Output request (outGrp_calcium)                ║
║          └─→ Call post::calcium_concentration()         ║
║          └─→ Extract Xion(3, Ac) for each node          ║
║                                                          ║
║  Step 3: VTK output                                     ║
║          └─→ Calcium_concentration field                ║
║                                                          ║
╚══════════════════════════════════════════════════════════╝
```

## Implementation Components

```
┌─────────────────────────────────────────────────────────┐
│              CONSTS.H (Enumerations)                    │
├─────────────────────────────────────────────────────────┤
│                                                         │
│  enum class OutputNameType {                           │
│    ...                                                 │
│    outGrp_calcium = 529,  ← Output group identifier   │
│    ...                                                 │
│    out_calcium = 568       ← Output type identifier   │
│    ...                                                 │
│  }                                                    │
│                                                         │
└─────────────────────────────────────────────────────────┘
            │                         │
            ↓                         ↓
    ┌───────────────┐        ┌──────────────────┐
    │SET_OUTPUT_    │        │SET_EQUATION_     │
    │PROPS.H        │        │PROPS.H           │
    │               │        │                  │
    │output_props_  │        │CEP Equation:     │
    │map =          │        │                  │
    │{out_calcium→  │        │nDOP = {1, 2, ..}│
    │ outGrp_calc...│        │outPuts[1]=       │
    │ 0,1,         │        │ out_calcium      │
    │"Calcium_..."}│        │                  │
    └───────────────┘        └──────────────────┘
            │                         │
            └────────┬────────────────┘
                     ↓
            ┌──────────────────┐
            │   POST.H & CPP   │
            │                  │
            │ calcium_          │
            │concentration()    │
            │ {                │
            │  res(a) =       │
            │  Xion(3, Ac)    │
            │ }               │
            │                  │
            │ + case handler   │
            │   for           │
            │   outGrp_calcium │
            └──────────────────┘
                     ↓
            ┌──────────────────┐
            │   VTK OUTPUT     │
            │                  │
            │ Calcium_         │
            │concentration     │
            │ [mM per node]    │
            └──────────────────┘
```

## Time Step Output Sequence

```
Time step N:
├─ Solve CEP equations
│  └─ Update Xion array with all state variables
│     └─ Xion(3, :) = new calcium values
│
├─ Check output frequency
│  └─ Is this an output step?
│
├─ If YES:
│  ├─ FOR each output in outPuts array:
│  │  │
│  │  ├─ If outPuts[i] == out_voltage
│  │  │  └─ ReadVoltage from lY directly
│  │  │     └─ Output as "Action_potential"
│  │  │
│  │  ├─ If outPuts[i] == out_calcium ← HERE
│  │  │  └─ Call calcium_concentration()
│  │  │     ├─ FOR each node:
│  │  │     │  res(a) = Xion(3, gN(a))
│  │  │     └─ Return result array
│  │  │     └─ Output as "Calcium_concentration"
│  │  │
│  │  └─ ... other outputs ...
│  │
│  └─ Write VTK file with all results
│     └─ result_NNN.vtu contains both:
│        ├─ Action_potential
│        └─ Calcium_concentration
│
└─ Continue to next time step
```

## Summary

The implementation adds a **calcium extraction layer** between the CEP solver's internal state arrays and the VTK output. It:

1. **Reads** from `cep_mod.Xion(3, node)` 
2. **Processes** via `calcium_concentration()` function
3. **Outputs** as "Calcium_concentration" field in VTK

All integrated seamlessly using the existing output framework!

---

**Key Insight**: Calcium is already being computed by the CEP solver - we're just extracting it and making it available for output.
