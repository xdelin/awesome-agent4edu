#!/usr/bin/env python3
"""
MLX Memory Debugging Utility

Usage (with uv):
    uv run python check_memory.py                    # Show current memory stats
    uv run python check_memory.py --watch            # Monitor memory continuously
    uv run python check_memory.py --clear            # Clear cache and show stats
    uv run python check_memory.py --watch --log memory.csv  # Log to CSV file

Or run directly if mlx is in your environment:
    python check_memory.py
"""

import argparse
import time
import os

try:
    import mlx.core as mx
except ImportError:
    print("Error: MLX not installed. Install with: uv add mlx")
    exit(1)


def format_bytes(bytes_val: int) -> str:
    """Format bytes as human-readable string."""
    if bytes_val >= 1e9:
        return f"{bytes_val / 1e9:.2f} GB"
    elif bytes_val >= 1e6:
        return f"{bytes_val / 1e6:.2f} MB"
    elif bytes_val >= 1e3:
        return f"{bytes_val / 1e3:.2f} KB"
    return f"{bytes_val} B"


def get_memory_stats() -> dict:
    """Get current MLX memory statistics."""
    return {
        "active": mx.get_active_memory(),
        "peak": mx.get_peak_memory(),
        "cache": mx.get_cache_memory(),
    }


def print_stats(stats: dict, label: str = "") -> None:
    """Print memory statistics."""
    prefix = f"[{label}] " if label else ""
    print(f"{prefix}Active: {format_bytes(stats['active'])}, "
          f"Peak: {format_bytes(stats['peak'])}, "
          f"Cache: {format_bytes(stats['cache'])}")


def main():
    parser = argparse.ArgumentParser(description="MLX Memory Debugging Utility")
    parser.add_argument("--watch", action="store_true",
                        help="Monitor memory continuously (Ctrl+C to stop)")
    parser.add_argument("--interval", type=float, default=1.0,
                        help="Watch interval in seconds (default: 1.0)")
    parser.add_argument("--clear", action="store_true",
                        help="Clear memory cache before showing stats")
    parser.add_argument("--reset-peak", action="store_true",
                        help="Reset peak memory counter")
    parser.add_argument("--log", type=str, metavar="FILE",
                        help="Log memory stats to CSV file (use with --watch)")
    args = parser.parse_args()

    if args.reset_peak:
        mx.reset_peak_memory()
        print("Peak memory counter reset.")

    if args.clear:
        before = get_memory_stats()
        mx.clear_cache()
        after = get_memory_stats()
        print("Cache cleared.")
        print(f"  Before: Cache = {format_bytes(before['cache'])}")
        print(f"  After:  Cache = {format_bytes(after['cache'])}")
        print()

    if args.watch:
        print("Monitoring MLX memory (Ctrl+C to stop)...")
        if args.log:
            print(f"Logging to: {args.log}")
            # Write CSV header if file doesn't exist
            if not os.path.exists(args.log):
                with open(args.log, "w") as f:
                    f.write("timestamp,active_bytes,peak_bytes,cache_bytes\n")
        print("-" * 60)
        try:
            while True:
                stats = get_memory_stats()
                timestamp = time.strftime("%H:%M:%S")
                print_stats(stats, timestamp)
                if args.log:
                    with open(args.log, "a") as f:
                        f.write(f"{timestamp},{stats['active']},{stats['peak']},{stats['cache']}\n")
                time.sleep(args.interval)
        except KeyboardInterrupt:
            print("\nStopped.")
    else:
        stats = get_memory_stats()
        print("MLX Memory Stats:")
        print(f"  Active:  {format_bytes(stats['active'])}")
        print(f"  Peak:    {format_bytes(stats['peak'])}")
        print(f"  Cache:   {format_bytes(stats['cache'])}")


if __name__ == "__main__":
    main()
