#!/usr/bin/env python3
"""
Post-process and aggregate results from parametric sweep.

This script collects results from all completed jobs and generates
summary statistics, CSV files, and optional plots.

Usage:
    python post_process.py results/spherical_k_sweep
    python post_process.py results/my_sweep --no-plots
"""

import argparse
import sys
import json
import pandas as pd
from pathlib import Path
from typing import List, Dict, Any
import warnings


def find_completed_jobs(results_dir: Path) -> List[Path]:
    """Find all completed job directories with metadata."""
    job_dirs = []
    for job_dir in sorted(results_dir.glob("job_*")):
        metadata_file = job_dir / "metadata.json"
        if metadata_file.exists():
            job_dirs.append(job_dir)
    return job_dirs

def find_error_jobs(results_dir: Path) -> List[Path]:
    """Find all job directories that have an error log."""
    error_jobs = {"meshing_failed": [], "error": [], "segment": [], "other": []}
    for job_file in sorted(results_dir.glob("*.err")):
        job_nbr = int(job_file.stem.split("/")[-1])
        
        with open(job_file, 'r') as f:
            content = f.read()

            if len(content.strip()) > 131:
                if "Meshing failed!" in content:
                    # print(f"Meshing failed in {job_nbr}")
                    error_jobs["meshing_failed"].append(job_nbr)
                elif "Error" in content:
                    # print(f"Error in {job_nbr}")
                    error_jobs["error"].append(job_nbr)
                elif "Segmentation" in content:
                    error_jobs["segment"].append(job_nbr)
                else:
                    # print(f"Error Msg in {job_file}")
                    error_jobs["other"].append(job_nbr)
    
    return error_jobs


def load_job_metadata(job_dir: Path) -> Dict[str, Any]:
    """Load metadata from a job directory."""
    metadata_file = job_dir / "metadata.json"
    with open(metadata_file, 'r') as f:
        return json.load(f)


def flatten_dict(d: Dict[str, Any], parent_key: str = '', sep: str = '_') -> Dict[str, Any]:
    """Flatten nested dictionary for easier DataFrame conversion."""
    items = []
    for k, v in d.items():
        new_key = f"{parent_key}{sep}{k}" if parent_key else k
        if isinstance(v, dict):
            items.extend(flatten_dict(v, new_key, sep=sep).items())
        elif isinstance(v, list):
            # Convert lists to strings for CSV compatibility
            items.append((new_key, str(v)))
        else:
            items.append((new_key, v))
    return dict(items)


def aggregate_results(job_dirs: List[Path]) -> pd.DataFrame:
    """
    Aggregate metadata from all jobs into a pandas DataFrame.

    Returns:
        DataFrame with one row per job, columns for all parameters and results
    """
    records = []
    for job_dir in job_dirs:
        try:
            metadata = load_job_metadata(job_dir)
            flat_record = flatten_dict(metadata)
            flat_record['job_dir'] = str(job_dir)
            records.append(flat_record)
        except Exception as e:
            warnings.warn(f"Error loading {job_dir}: {e}")
            continue

    if not records:
        raise ValueError("No valid job results found")

    df = pd.DataFrame(records)
    return df


def generate_summary_statistics(df: pd.DataFrame) -> Dict[str, Any]:
    """Generate summary statistics from results DataFrame."""
    summary = {
        'total_jobs': len(df),
        'timing_statistics': {},
        'mesh_statistics': {},
    }

    # Timing statistics
    if 'timings_total' in df.columns:
        summary['timing_statistics'] = {
            'mean_total_time': float(df['timings_total'].mean()),
            'median_total_time': float(df['timings_total'].median()),
            'min_total_time': float(df['timings_total'].min()),
            'max_total_time': float(df['timings_total'].max()),
            'std_total_time': float(df['timings_total'].std()),
        }

    # Mesh statistics
    if 'mesh_ndof' in df.columns:
        summary['mesh_statistics'] = {
            'mean_ndof': float(df['mesh_ndof'].mean()),
            'median_ndof': float(df['mesh_ndof'].median()),
            'min_ndof': int(df['mesh_ndof'].min()),
            'max_ndof': int(df['mesh_ndof'].max()),
        }

    return summary


def save_results(df: pd.DataFrame, summary: Dict[str, Any], results_dir: Path):
    """Save aggregated results and summary."""
    # Save CSV
    csv_file = results_dir / "summary.csv"
    df.to_csv(csv_file, index=False)
    print(f"✓ Saved summary CSV: {csv_file}")

    # Save JSON summary
    summary_file = results_dir / "summary_statistics.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"✓ Saved statistics: {summary_file}")


def print_summary(df: pd.DataFrame, summary: Dict[str, Any]):
    """Print summary to console."""
    print("\n" + "=" * 70)
    print("Post-Processing Summary")
    print("=" * 70)

    print(f"\nTotal jobs processed: {summary['total_jobs']}")

    if summary['timing_statistics']:
        print("\nTiming Statistics:")
        ts = summary['timing_statistics']
        print(f"  Mean total time:   {ts['mean_total_time']:.2f} s")
        print(f"  Median total time: {ts['median_total_time']:.2f} s")
        print(f"  Min/Max:           {ts['min_total_time']:.2f} / {ts['max_total_time']:.2f} s")
        print(f"  Std deviation:     {ts['std_total_time']:.2f} s")

    if summary['mesh_statistics']:
        print("\nMesh Statistics:")
        ms = summary['mesh_statistics']
        print(f"  Mean DOFs:   {ms['mean_ndof']:.0f}")
        print(f"  Median DOFs: {ms['median_ndof']:.0f}")
        print(f"  Min/Max:     {ms['min_ndof']} / {ms['max_ndof']}")

    # Parameter ranges
    param_cols = [col for col in df.columns if col.startswith('parameters_')]
    if param_cols:
        print("\nParameter Ranges:")
        for col in param_cols:
            param_name = col.replace('parameters_', '')
            try:
                unique_vals = df[col].nunique()
                print(f"  {param_name:20s}: {unique_vals} unique values")
            except:
                pass

    print("=" * 70)


def generate_plots(df: pd.DataFrame, results_dir: Path):
    """Generate optional plots (requires matplotlib)."""
    try:
        import matplotlib.pyplot as plt
        import matplotlib
        matplotlib.use('Agg')  # Non-interactive backend
    except ImportError:
        print("\nNote: matplotlib not available, skipping plots")
        return

    # Plot 1: Timing vs wavelength
    if 'parameters_wavelength' in df.columns and 'timings_solve' in df.columns:
        try:
            fig, ax = plt.subplots(figsize=(10, 6))
            df_plot = df.sort_values('parameters_wavelength')
            ax.scatter(df_plot['parameters_wavelength'], df_plot['timings_solve'])
            ax.set_xlabel('Wavelength')
            ax.set_ylabel('Solve Time (s)')
            ax.set_title('Solve Time vs Wavelength')
            ax.grid(True, alpha=0.3)
            plot_file = results_dir / "solve_time_vs_wavelength.pdf"
            fig.savefig(plot_file, bbox_inches='tight')
            plt.close(fig)
            print(f"✓ Saved plot: {plot_file}")
        except Exception as e:
            warnings.warn(f"Error creating plot: {e}")

    # Plot 2: Timing breakdown
    timing_cols = [col for col in df.columns if col.startswith('timings_') and col != 'timings_total']
    if timing_cols:
        try:
            fig, ax = plt.subplots(figsize=(10, 6))
            timing_means = df[timing_cols].mean()
            timing_labels = [col.replace('timings_', '') for col in timing_cols]
            ax.bar(timing_labels, timing_means)
            ax.set_ylabel('Mean Time (s)')
            ax.set_title('Mean Timing Breakdown')
            ax.grid(True, alpha=0.3, axis='y')
            plt.xticks(rotation=45, ha='right')
            plot_file = results_dir / "timing_breakdown.pdf"
            fig.savefig(plot_file, bbox_inches='tight')
            plt.close(fig)
            print(f"✓ Saved plot: {plot_file}")
        except Exception as e:
            warnings.warn(f"Error creating plot: {e}")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Post-process and aggregate simulation results"
    )
    parser.add_argument(
        'results_dir',
        type=str,
        help='Results directory containing job_XXXX subdirectories'
    )
    parser.add_argument(
        'log_dir',
        type=str,
        help='Log directory containing XXXX.out and XXXX.err output files.'
    )
    parser.add_argument(
        '--no-plots',
        action='store_true',
        help='Skip plot generation'
    )

    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    if not results_dir.exists():
        print(f"Error: Results directory not found: {results_dir}", file=sys.stderr)
        return 1

    log_dir = Path(args.log_dir)

    print("=" * 70)
    print(f"Post-Processing: {results_dir}")
    print("=" * 70)

    # Find completed jobs
    print("\nScanning for completed jobs...")
    job_dirs = find_completed_jobs(results_dir)
    print(f"✓ Found {len(job_dirs)} completed jobs")

    if not job_dirs:
        print("Error: No completed jobs found", file=sys.stderr)
        return 1

    # Aggregate results
    print("\nAggregating results...")
    try:
        df = aggregate_results(job_dirs)
        print(f"✓ Aggregated {len(df)} job results")
    except Exception as e:
        print(f"Error aggregating results: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1

    # Generate summary statistics
    try:
        summary = generate_summary_statistics(df)
    except Exception as e:
        print(f"Error generating summary: {e}", file=sys.stderr)
        summary = {'total_jobs': len(df)}

    # Save results
    try:
        save_results(df, summary, results_dir)
    except Exception as e:
        print(f"Error saving results: {e}", file=sys.stderr)
        return 1

    # Generate plots (if requested)
    if not args.no_plots:
        print("\nGenerating plots...")
        generate_plots(df, results_dir)

    # Failed jobs
    error_jobs = find_error_jobs(log_dir)

    if len(error_jobs["meshing_failed"]) > 0 or len(error_jobs["error"]) > 0 or len(error_jobs["other"]) > 0:
            
        print("\n"+"=" * 70)
        print(f"Failed jobs summary:")
        print("=" * 70)

        fail_sum = 0

        print("\nFound failed jobs:")

        if len(error_jobs['meshing_failed']) > 0:
            fail_sum += len(error_jobs['meshing_failed'])
            print(f" - Meshing failed in {len(error_jobs['meshing_failed'])} jobs: {error_jobs['meshing_failed'][:100]}")
        if len(error_jobs['error']) > 0:
            fail_sum += len(error_jobs['error'])
            print(f" - Errors in {len(error_jobs['error'])} jobs: {error_jobs['error'][:100]}")
        if len(error_jobs['segment']) > 0:
            fail_sum += len(error_jobs['segment'])
            print(f" - Segmentation fault in {len(error_jobs['segment'])} jobs: {error_jobs['segment'][:100]}")
        if len(error_jobs['other']) > 0:
            fail_sum += len(error_jobs['other'])
            print(f" - Other error msg in {len(error_jobs['other'])} jobs: {error_jobs['other'][:100]}")

        print(f"\nTotal number of jobs: {fail_sum}")

    # Print summary
    print_summary(df, summary)

    print("\n✓ Post-processing completed successfully\n")
    return 0


if __name__ == '__main__':
    sys.exit(main())
