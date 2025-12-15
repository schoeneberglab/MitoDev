"""
Main script for MitoDev cell tracking pipeline.

This script orchestrates the cell tracking workflow:
1. Load configuration
2. Track cells across time frames
3. Fill cell volumes with smoothing
4. Save cell volumes to disk
"""

import os
import sys
import argparse
import os.path as osp
from typing import Dict, Any, List

from utils import load_config, get_paths, save_cells
from tracking import track, fill_cells_across_frames


def parse_arguments() -> argparse.Namespace:
    """
    Parse command line arguments.
    
    Returns:
        Namespace object containing parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description='MitoDev cell tracking pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python main.py config.yaml data_folder region
  python main.py /path/to/config.yaml sample1 region1
        """
    )
    
    parser.add_argument(
        'config_path',
        type=str,
        help='Path to the YAML configuration file'
    )
    
    parser.add_argument(
        'data_folder',
        type=str,
        help='Name of the data folder'
    )
    
    parser.add_argument(
        'region',
        type=str,
        help='Name of the region'
    )
    
    parser.add_argument(
        '--minimal',
        action='store_true',
        help='Only save 10 cells'
    )
    
    return parser.parse_args()


def main() -> None:
    """
    Main execution function for the cell tracking pipeline.
    
    Expects command line arguments:
        config_path: Path to YAML configuration file
        data_folder: Name of the data folder
        region: Name of the region
    """
    # Parse command line arguments
    args = parse_arguments()
    
    # Load configuration
    if not os.path.exists(args.config_path):
        raise FileNotFoundError(f"Configuration file not found: {args.config_path}")
    
    cfg = load_config(args.config_path)
    
    # Extract configuration parameters
    root_dir = cfg['data']['root_dir']
    save_root_dir = cfg['data']['save_dir']
    num_frames = cfg['data']['num_frames']
    min_z_len = cfg['tracking']['min_z_len']
    z_window = cfg['cell_smoothing']['z_window']
    sigma = cfg['cell_smoothing']['sigma']
    debug = cfg.get('debug', {'flag': False})
    resume = cfg.get('resume', {}).get('flag', False)
    minimal = args.minimal
    
    data_folder = args.data_folder
    region = args.region

    # Set up metadata directory
    metadata_dir = os.path.join(save_root_dir, data_folder, region, 'extractions')
    save_root_dir = os.path.join(save_root_dir, data_folder, region, 'single_cells')
    
    # Get file paths
    try:
        line_centers_path, mito_files_path, masks_path, center_cellmask_maps = get_paths(
            root_dir, metadata_dir, data_folder, region
        )
    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)

    # Get z-stacks for tracking
    z_stacks = cfg['tracking']['zstacks']

    # Track cells across frames
    print("=" * 60)
    print("Starting cell tracking across frames...")
    print("=" * 60)
    try:
        cell_pts, cell_labels = track(
            line_centers_path=line_centers_path,
            num_frames=num_frames, 
            z_stacks=z_stacks,
            min_z_len=min_z_len
        )
        print(f"Successfully tracked {len(cell_labels)} cells across {num_frames} frames")
    except Exception as e:
        print(f"Error during tracking: {e}")
        sys.exit(1)

    # Determine frame range
    if debug.get('flag', False):
        cell_labels = debug.get('cell_labels', cell_labels)
        frame_range = debug.get('frame_range', [0, num_frames])
    else:
        frame_range = [0, num_frames]

    # Fill cells across frames with smoothing
    print("=" * 60)
    print("Filling cell volumes across frames...")
    print("=" * 60)
    os.makedirs(save_root_dir, exist_ok=True)
    
    try:
        cell_dict = fill_cells_across_frames(
            cell_pts, 
            cell_labels, 
            frame_range,
            center_cellmask_maps, 
            save_root_dir,
            data_folder, 
            region,
            z_window=z_window, 
            sigma=sigma, 
            resume=resume,
            minimal=minimal
        )
        print(f"Successfully processed {len(cell_dict)} cells")
    except Exception as e:
        print(f"Error during cell filling: {e}")
        sys.exit(1)

    # Save cell volumes
    print("=" * 60)
    print("Saving cell volumes...")
    print("=" * 60)
    try:
        save_cells(
            cell_dict, 
            mito_files_path, 
            data_folder, 
            region, 
            save_root_dir,
            frame_range
        )
        print("Successfully saved all cell volumes")
    except Exception as e:
        print(f"Error during cell saving: {e}")
        sys.exit(1)
    
    print("=" * 60)
    print("Pipeline completed successfully!")
    print("=" * 60)


if __name__ == '__main__':
    main()
