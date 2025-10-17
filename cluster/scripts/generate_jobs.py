#!/usr/bin/env python3
"""
Generate job parameter file from YAML configuration.

This script reads a YAML configuration file defining a parametric sweep
and generates a JSON file containing all parameter combinations for each job.

Usage:
    python generate_jobs.py configs/spherical_sweep.yaml
    python generate_jobs.py configs/my_sweep.yaml --validate-only
"""

import argparse
import sys
import json
import yaml
import numpy as np
from pathlib import Path
from itertools import product
from typing import Dict, List, Any


def load_config(config_path: Path) -> Dict[str, Any]:
    """Load YAML configuration file."""
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config


def expand_parameter_values(param_def: Dict[str, Any]) -> List[Any]:
    """
    Expand parameter definition to list of values.

    Supports:
    - values: [list of discrete values]
    - range: {start, stop, num} for linspace
    """
    if 'values' in param_def:
        return param_def['values']
    elif 'range' in param_def:
        r = param_def['range']
        return np.linspace(r['start'], r['stop'], r['num']).tolist()
    else:
        raise ValueError(f"Parameter definition must have 'values' or 'range': {param_def}")


def generate_parameter_grid(config: Dict[str, Any]) -> List[Dict[str, Any]]:
    """
    Generate all parameter combinations from config.

    Returns:
        List of parameter dictionaries, one per job
    """
    parameters = config.get('parameters', {})

    # Extract parameter names and their values
    param_names = []
    param_values = []

    for name, definition in parameters.items():
        param_names.append(name)
        param_values.append(expand_parameter_values(definition))

    # Generate Cartesian product of all parameter combinations
    combinations = list(product(*param_values))

    # Create job parameter dictionaries
    jobs = []
    for job_id, combo in enumerate(combinations):
        job_params = {
            'job_id': job_id,
            'parameters': dict(zip(param_names, combo)),
            'geometry': config.get('geometry', {}),
            'solver': config.get('solver', {}),
            'output': config.get('output', {}),
        }
        jobs.append(job_params)

    return jobs


def save_job_list(jobs: List[Dict[str, Any]], output_path: Path):
    """Save job list as JSON file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(jobs, f, indent=2)


def print_summary(config: Dict[str, Any], jobs: List[Dict[str, Any]]):
    """Print summary of generated job list."""
    print("=" * 70)
    print(f"Sweep: {config['sweep_name']}")
    print("=" * 70)

    # Count values for each parameter
    parameters = config.get('parameters', {})
    print("\nParameter ranges:")
    total_combinations = 1
    for name, definition in parameters.items():
        values = expand_parameter_values(definition)
        count = len(values)
        total_combinations *= count
        print(f"  {name:20s}: {count:4d} values")

    print(f"\nTotal combinations: {total_combinations}")
    print(f"Jobs generated:     {len(jobs)}")

    # Resource summary
    resources = config.get('resources', {})
    print(f"\nResources per job:")
    print(f"  CPUs:        {resources.get('cpus', 'N/A')}")
    print(f"  Memory:      {resources.get('mem_per_cpu', 'N/A')} per CPU")
    print(f"  Time limit:  {resources.get('time', 'N/A')}")
    print(f"  Partition:   {resources.get('partition', 'N/A')}")

    # Output location
    print(f"\nOutput directory: {config.get('output_dir', 'N/A')}")
    print("=" * 70)


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Generate job parameter file from YAML configuration"
    )
    parser.add_argument(
        'config',
        type=str,
        help='Path to YAML configuration file'
    )
    parser.add_argument(
        '-o', '--output',
        type=str,
        default=None,
        help='Output JSON file (default: auto-generated from config)'
    )
    parser.add_argument(
        '--validate-only',
        action='store_true',
        help='Only validate config and print summary, do not write output'
    )

    args = parser.parse_args()

    # Load configuration
    config_path = Path(args.config)
    if not config_path.exists():
        print(f"Error: Config file not found: {config_path}", file=sys.stderr)
        return 1

    try:
        config = load_config(config_path)
    except Exception as e:
        print(f"Error loading config: {e}", file=sys.stderr)
        return 1

    # Generate job list
    try:
        jobs = generate_parameter_grid(config)
    except Exception as e:
        print(f"Error generating jobs: {e}", file=sys.stderr)
        return 1

    # Print summary
    print_summary(config, jobs)

    # Determine output path
    if args.output:
        output_path = Path(args.output)
    else:
        sweep_name = config.get('sweep_name', 'sweep')
        output_path = config_path.parent / f"{sweep_name}_jobs.json"

    # Save job list (unless validate-only)
    if not args.validate_only:
        try:
            save_job_list(jobs, output_path)
            print(f"\n✓ Job list saved to: {output_path}")
            print(f"\nNext steps:")
            print(f"  LOCAL TEST:")
            print(f"    python scripts/run_local_test.py {sweep_name}")
            print(f"  CLUSTER:")
            print(f"    sbatch --array=0-{len(jobs)-1} slurm/array_job.slurm {sweep_name}")
        except Exception as e:
            print(f"Error saving job list: {e}", file=sys.stderr)
            return 1
    else:
        print("\n✓ Validation successful (no files written)")

    return 0


if __name__ == '__main__':
    sys.exit(main())
