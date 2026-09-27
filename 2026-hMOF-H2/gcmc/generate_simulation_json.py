#!/usr/bin/env python3
"""Create a RASPA3 simulation.json for one temperature/pressure state.

The template is an actual production input. Only the framework name and the
external state are changed; all sampling and charge settings are preserved.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


ALLOWED_TEMPERATURES = (77, 120, 160, 200, 233, 253, 273, 298)
ALLOWED_PRESSURES_BAR = (5, 100)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--template", type=Path, required=True)
    parser.add_argument("--framework-name", required=True)
    parser.add_argument("--temperature", type=int, required=True)
    parser.add_argument("--pressure-bar", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    if args.temperature not in ALLOWED_TEMPERATURES:
        parser.error(f"temperature must be one of {ALLOWED_TEMPERATURES}")
    if args.pressure_bar not in ALLOWED_PRESSURES_BAR:
        parser.error(f"pressure-bar must be one of {ALLOWED_PRESSURES_BAR}")

    data = json.loads(args.template.read_text(encoding="utf-8"))
    systems = data.get("Systems")
    if not isinstance(systems, list) or len(systems) != 1:
        parser.error("template must contain exactly one RASPA3 framework system")
    system = systems[0]
    system["Name"] = args.framework_name
    system["ExternalTemperature"] = args.temperature
    system["ExternalPressure"] = args.pressure_bar * 100000

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(data, indent=4) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
