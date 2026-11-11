#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
# University of California, and others. SPDX-License-Identifier: BSD-3-Clause

"""Generate TP06 + ultra-slow I_Na gate (TTP_slowNa) FE/Rush--Larsen CSVs.

The CellML TP06 cores and the protocol of generate_ttp.py are reused
unchanged. The driver adds a 13th gate k, following the openCARP model
tenTusscherPanfilov_slow_Na:

    k_inf = 1 / (1 + exp((V - V_half_k) / slope_k))
    tau_k = 5 + usi_param / (1 + exp((V + 15) / 4))
    I_Na  = G_Na m^3 h j k^3 (V - E_Na)

The generated computeRates returns I_Na without k (algebraic[50]). The
driver replaces that term in dV/dt and dNai/dt with k^3 * I_Na.
"""

from __future__ import annotations

import argparse
import csv
import sys
from math import exp, isfinite
from pathlib import Path
from typing import Iterator, TextIO

import generate_ttp as ttp


USI_PARAM = 1050.0
V_HALF_K = 10.0
SLOPE_K = 5.0
K_INITIAL = 1.0

STATE_COLUMNS = (*ttp.STATE_COLUMNS, "k")

# Indices into the generated CellML arrays.
V_INDEX = 0
NAI_INDEX = 2
I_NA_ALGEBRAIC_INDEX = 50
F_CONSTANT_INDEX = 2
CM_CONSTANT_INDEX = 3
VC_CONSTANT_INDEX = 4


def k_kinetics(voltage: float) -> tuple[float, float]:
    """Return (k_inf, tau_k) for the ultra-slow I_Na gate."""

    steady_state = 1.0 / (1.0 + exp((voltage - V_HALF_K) / SLOPE_K))
    time_constant = 5.0 + USI_PARAM / (1.0 + exp((voltage + 15.0) / 4.0))
    return steady_state, time_constant


def iter_full_trajectory(
    profile: str = "epi",
) -> Iterator[tuple[int, tuple[float, ...]]]:
    """Yield every completed-update state from the canonical protocol."""

    configuration = ttp.PROFILE_CONFIGURATIONS[profile]
    states, constants = configuration.model.initConsts()
    if len(states) != 19:
        raise AssertionError(f"expected 19 states, found {len(states)}")
    k = K_INITIAL

    sodium_concentration_factor = constants[CM_CONSTANT_INDEX] / (
        constants[VC_CONSTANT_INDEX] * constants[F_CONSTANT_INDEX]
    )

    for step in range(ttp.UPDATE_COUNT + 1):
        values = (*(states[index] for index in ttp.OUTPUT_STATE_INDICES), k)
        if not all(isfinite(value) for value in values):
            raise AssertionError(f"non-finite state after update {step}")
        yield step, values

        if step == ttp.UPDATE_COUNT:
            continue

        old_time_ms = step * ttp.DT_MS
        rates, algebraic = configuration.model.computeRates(
            old_time_ms, states, constants, ttp.stimulus_current(old_time_ms)
        )

        # Replace I_Na with k^3 * I_Na in dV/dt and dNai/dt.
        sodium_current_change = algebraic[I_NA_ALGEBRAIC_INDEX] * (k**3 - 1.0)
        rates[V_INDEX] -= sodium_current_change
        rates[NAI_INDEX] -= sodium_current_change * sodium_concentration_factor

        updated_states = list(states)
        for state_index in ttp.MAIN_STATE_INDICES:
            updated_states[state_index] = (
                states[state_index] + ttp.DT_MS * rates[state_index]
            )
        for state_index, steady_state_index, time_constant_index in ttp.GATE_KINETICS:
            steady_state = algebraic[steady_state_index]
            time_constant = algebraic[time_constant_index]
            updated_states[state_index] = steady_state - (
                steady_state - states[state_index]
            ) * exp(-ttp.DT_MS / time_constant)

        k_steady_state, k_time_constant = k_kinetics(states[V_INDEX])
        k = k_steady_state - (k_steady_state - k) * exp(-ttp.DT_MS / k_time_constant)

        states = updated_states


def generate_canonical_trajectory(profile: str = "epi") -> ttp.CanonicalTrajectory:
    """Run the 600 ms protocol and collect the profile's checkpoints."""

    checkpoints = ttp.PROFILE_CONFIGURATIONS[profile].checkpoints
    checkpoint_set = set(checkpoints)
    rows = tuple(
        (step, values)
        for step, values in iter_full_trajectory(profile)
        if step in checkpoint_set
    )
    if tuple(step for step, _ in rows) != checkpoints:
        raise AssertionError("not every canonical checkpoint was recorded")
    return ttp.CanonicalTrajectory(rows)


def write_csv(trajectory: ttp.CanonicalTrajectory, stream: TextIO) -> None:
    writer = csv.writer(stream, lineterminator="\n")
    writer.writerow(("step", *STATE_COLUMNS))
    for step, values in trajectory.rows:
        writer.writerow((step, *values))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, help="write canonical CSV to this path instead of stdout"
    )
    parser.add_argument(
        "--profile", choices=tuple(ttp.PROFILE_CONFIGURATIONS), default="epi"
    )
    args = parser.parse_args()

    trajectory = generate_canonical_trajectory(args.profile)
    if args.output is None:
        write_csv(trajectory, sys.stdout)
    else:
        with args.output.open("w", encoding="utf-8", newline="") as stream:
            write_csv(trajectory, stream)


if __name__ == "__main__":
    main()
