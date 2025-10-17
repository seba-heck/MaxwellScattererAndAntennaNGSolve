#!/usr/bin/env python3
"""
Run local test of cluster workflow without SLURM.

This script runs a subset of jobs locally for testing and development.
It simulates what would happen on the cluster but uses local resources.

Usage:
    python run_local_test.py test_local
    python run_local_test.py my_sweep --jobs 0,1,2
    python run_local_test.py my_sweep --all  # Run all jobs locally
"""

import argparse
import sys
import subprocess
from pathlib import Path


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Run local test of cluster workflow"
    )
    parser.add_argument(
        'sweep_name',
        type=str,
        help='Sweep name (e.g., test_local)'
    )
    parser.add_argument(
        '--jobs',
        type=str,
        default='0',
        help='Comma-separated job IDs to run (default: 0)'
    )
    parser.add_argument(
        '--all',
        action='store_true',
        help='Run all jobs in the sweep'
    )
    parser.add_argument(
        '--threads',
        type=int,
        default=4,
        help='Number of threads per job (default: 4)'
    )

    args = parser.parse_args()

    # Paths
    script_dir = Path(__file__).parent
    cluster_dir = script_dir.parent
    configs_dir = cluster_dir / "configs"
    results_dir = cluster_dir / "results" / args.sweep_name

    # Job file
    job_file = configs_dir / f"{args.sweep_name}_jobs.json"

    if not job_file.exists():
        print(f"Error: Job file not found: {job_file}", file=sys.stderr)
        print(f"Did you run: python scripts/generate_jobs.py configs/{args.sweep_name}.yaml ?", file=sys.stderr)
        return 1

    # Determine which jobs to run
    if args.all:
        import json
        with open(job_file, 'r') as f:
            jobs_data = json.load(f)
        job_ids = list(range(len(jobs_data)))
        print(f"Running ALL {len(job_ids)} jobs locally...")
    else:
        job_ids = [int(x.strip()) for x in args.jobs.split(',')]
        print(f"Running {len(job_ids)} job(s) locally: {job_ids}")

    print("=" * 70)
    print(f"Local Test: {args.sweep_name}")
    print("=" * 70)
    print(f"Job file:      {job_file}")
    print(f"Results dir:   {results_dir}")
    print(f"Threads/job:   {args.threads}")
    print(f"Jobs to run:   {job_ids}")
    print("=" * 70)

    # Run each job
    failed_jobs = []
    for i, job_id in enumerate(job_ids):
        print(f"\n{'='*70}")
        print(f"Running job {job_id} ({i+1}/{len(job_ids)})...")
        print(f"{'='*70}")

        cmd = [
            sys.executable,
            str(script_dir / "run_simulation.py"),
            "--job-file", str(job_file),
            "--job-id", str(job_id),
            "--output-dir", str(results_dir),
            "--num-threads", str(args.threads),
            "--local"
        ]

        try:
            result = subprocess.run(cmd, check=True)
            print(f"\n✓ Job {job_id} completed successfully")
        except subprocess.CalledProcessError as e:
            print(f"\n✗ Job {job_id} failed with exit code {e.returncode}", file=sys.stderr)
            failed_jobs.append(job_id)
        except KeyboardInterrupt:
            print("\n\nInterrupted by user", file=sys.stderr)
            return 130

    # Summary
    print("\n" + "=" * 70)
    print("Local Test Summary")
    print("=" * 70)
    print(f"Total jobs:      {len(job_ids)}")
    print(f"Successful:      {len(job_ids) - len(failed_jobs)}")
    print(f"Failed:          {len(failed_jobs)}")
    if failed_jobs:
        print(f"Failed job IDs:  {failed_jobs}")
    print("=" * 70)

    # Post-process if any jobs completed
    if len(failed_jobs) < len(job_ids):
        print("\nRunning post-processing...")
        try:
            subprocess.run([
                sys.executable,
                str(script_dir / "post_process.py"),
                str(results_dir)
            ], check=True)
        except subprocess.CalledProcessError:
            print("Warning: Post-processing failed", file=sys.stderr)

    if failed_jobs:
        return 1

    print("\n✓ All local tests passed!")
    return 0


if __name__ == '__main__':
    sys.exit(main())
