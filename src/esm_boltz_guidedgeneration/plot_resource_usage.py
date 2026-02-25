#!/usr/bin/env python3
"""
Plot GPU and CPU usage metrics from log files.

This script generates visualizations of resource utilization during
ESM3-Boltz guided generation jobs.

Usage:
    python plot_resource_usage.py GPU_LOG_FILE CPU_LOG_FILE OUTPUT_DIR

Author: Anna Su
Date: 2026-02-20
"""

import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from pathlib import Path
from datetime import datetime


def plot_gpu_usage(gpu_log_path, output_dir):
    """
    Generate GPU usage plots from log file.
    
    Args:
        gpu_log_path: Path to GPU usage CSV log
        output_dir: Directory to save plots
    
    Returns:
        Path to saved plot or None if failed
    """
    try:
        # Read GPU log
        df = pd.read_csv(gpu_log_path)
        
        if df.empty:
            print(f"Warning: GPU log file is empty: {gpu_log_path}")
            return None
        
        # Parse timestamp
        df['timestamp'] = pd.to_datetime(df['timestamp'])
        
        # Create figure with subplots
        fig, axes = plt.subplots(3, 1, figsize=(12, 10))
        fig.suptitle('GPU Resource Utilization', fontsize=16, fontweight='bold')
        
        # Plot 1: GPU Utilization
        axes[0].plot(df['timestamp'], df['gpu_util_percent'], 
                    color='#2E86AB', linewidth=2, marker='o', markersize=3)
        axes[0].set_ylabel('GPU Utilization (%)', fontsize=12, fontweight='bold')
        axes[0].set_ylim(0, 105)
        axes[0].grid(True, alpha=0.3)
        axes[0].axhline(y=df['gpu_util_percent'].mean(), 
                       color='red', linestyle='--', linewidth=1.5, 
                       label=f"Avg: {df['gpu_util_percent'].mean():.1f}%")
        axes[0].legend(loc='upper right')
        
        # Plot 2: GPU Memory Usage
        axes[1].plot(df['timestamp'], df['memory_used_mb'], 
                    color='#A23B72', linewidth=2, marker='s', markersize=3, label='Used')
        axes[1].axhline(y=df['memory_total_mb'].iloc[0], 
                       color='gray', linestyle='--', linewidth=1.5, 
                       label=f"Total: {df['memory_total_mb'].iloc[0]:.0f} MB")
        axes[1].set_ylabel('GPU Memory (MB)', fontsize=12, fontweight='bold')
        axes[1].grid(True, alpha=0.3)
        axes[1].legend(loc='upper right')
        
        # Plot 3: Temperature and Power
        ax3 = axes[2]
        ax3_twin = ax3.twinx()
        
        line1 = ax3.plot(df['timestamp'], df['temperature_c'], 
                        color='#F18F01', linewidth=2, marker='^', markersize=3, label='Temperature')
        ax3.set_ylabel('Temperature (°C)', fontsize=12, fontweight='bold', color='#F18F01')
        ax3.tick_params(axis='y', labelcolor='#F18F01')
        ax3.grid(True, alpha=0.3)
        
        line2 = ax3_twin.plot(df['timestamp'], df['power_draw_w'], 
                             color='#6A4C93', linewidth=2, marker='v', markersize=3, label='Power')
        ax3_twin.set_ylabel('Power Draw (W)', fontsize=12, fontweight='bold', color='#6A4C93')
        ax3_twin.tick_params(axis='y', labelcolor='#6A4C93')
        
        # Combine legends
        lines = line1 + line2
        labels = [l.get_label() for l in lines]
        ax3.legend(lines, labels, loc='upper right')
        
        # Format x-axis for all subplots
        for ax in axes:
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
            ax.tick_params(axis='x', rotation=45)
        
        axes[2].set_xlabel('Time', fontsize=12, fontweight='bold')
        
        # Adjust layout
        plt.tight_layout()
        
        # Save plot
        output_path = Path(output_dir) / 'gpu_usage_plot.png'
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close()
        
        print(f"✓ GPU usage plot saved to: {output_path}")
        return output_path
        
    except Exception as e:
        print(f"Error generating GPU usage plot: {e}")
        return None


def plot_cpu_usage(cpu_log_path, output_dir):
    """
    Generate CPU usage plots from log file.
    
    Args:
        cpu_log_path: Path to CPU usage CSV log
        output_dir: Directory to save plots
    
    Returns:
        Path to saved plot or None if failed
    """
    try:
        # Read CPU log
        df = pd.read_csv(cpu_log_path)
        
        if df.empty:
            print(f"Warning: CPU log file is empty: {cpu_log_path}")
            return None
        
        # Parse timestamp
        df['timestamp'] = pd.to_datetime(df['timestamp'])
        
        # Create figure with subplots
        fig, axes = plt.subplots(3, 1, figsize=(12, 10))
        fig.suptitle('CPU Resource Utilization', fontsize=16, fontweight='bold')
        
        # Plot 1: CPU Utilization
        axes[0].plot(df['timestamp'], df['cpu_util_percent'], 
                    color='#06A77D', linewidth=2, marker='o', markersize=3)
        axes[0].set_ylabel('CPU Utilization (%)', fontsize=12, fontweight='bold')
        axes[0].set_ylim(0, 105)
        axes[0].grid(True, alpha=0.3)
        axes[0].axhline(y=df['cpu_util_percent'].mean(), 
                       color='red', linestyle='--', linewidth=1.5, 
                       label=f"Avg: {df['cpu_util_percent'].mean():.1f}%")
        axes[0].legend(loc='upper right')
        
        # Plot 2: Memory Usage
        axes[1].plot(df['timestamp'], df['memory_used_mb'], 
                    color='#D62246', linewidth=2, marker='s', markersize=3, label='Used')
        axes[1].axhline(y=df['memory_total_mb'].iloc[0], 
                       color='gray', linestyle='--', linewidth=1.5, 
                       label=f"Total: {df['memory_total_mb'].iloc[0]:.0f} MB")
        axes[1].set_ylabel('Memory (MB)', fontsize=12, fontweight='bold')
        axes[1].grid(True, alpha=0.3)
        
        # Calculate memory percentage
        mem_pct = (df['memory_used_mb'] / df['memory_total_mb'] * 100).mean()
        axes[1].text(0.02, 0.95, f'Avg Usage: {mem_pct:.1f}%', 
                    transform=axes[1].transAxes, 
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
                    verticalalignment='top')
        axes[1].legend(loc='upper right')
        
        # Plot 3: Load Average
        axes[2].plot(df['timestamp'], df['load_avg_1min'], 
                    color='#F77F00', linewidth=2, marker='^', markersize=3)
        axes[2].set_ylabel('Load Average (1 min)', fontsize=12, fontweight='bold')
        axes[2].grid(True, alpha=0.3)
        axes[2].axhline(y=df['load_avg_1min'].mean(), 
                       color='red', linestyle='--', linewidth=1.5, 
                       label=f"Avg: {df['load_avg_1min'].mean():.2f}")
        axes[2].legend(loc='upper right')
        
        # Format x-axis for all subplots
        for ax in axes:
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
            ax.tick_params(axis='x', rotation=45)
        
        axes[2].set_xlabel('Time', fontsize=12, fontweight='bold')
        
        # Adjust layout
        plt.tight_layout()
        
        # Save plot
        output_path = Path(output_dir) / 'cpu_usage_plot.png'
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close()
        
        print(f"✓ CPU usage plot saved to: {output_path}")
        return output_path
        
    except Exception as e:
        print(f"Error generating CPU usage plot: {e}")
        return None


def main():
    if len(sys.argv) != 4:
        print("Usage: python plot_resource_usage.py GPU_LOG_FILE CPU_LOG_FILE OUTPUT_DIR")
        sys.exit(1)
    
    gpu_log_path = sys.argv[1]
    cpu_log_path = sys.argv[2]
    output_dir = sys.argv[3]
    
    # Create output directory if it doesn't exist
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    print("\n" + "="*60)
    print("Generating Resource Usage Plots")
    print("="*60)
    
    # Generate plots
    gpu_plot = plot_gpu_usage(gpu_log_path, output_dir) if Path(gpu_log_path).exists() else None
    cpu_plot = plot_cpu_usage(cpu_log_path, output_dir) if Path(cpu_log_path).exists() else None
    
    print("="*60)
    
    if gpu_plot or cpu_plot:
        print("✓ Resource usage plots generated successfully!")
    else:
        print("✗ Failed to generate resource usage plots")
        sys.exit(1)


if __name__ == "__main__":
    main()
