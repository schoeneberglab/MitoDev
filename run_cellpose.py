"""
CellPose segmentation script for MitoDev.

This script processes image stacks using CellPose to segment cells and extract:
- Cell masks
- Line centers (cell centroids)
- Center-to-cell-mask mappings

Expects command line arguments:
    sys.argv[1]: data_dir - Directory containing input images
    sys.argv[2]: start_idx - Starting index for file processing
    sys.argv[3]: end_idx - Ending index for file processing (exclusive)
"""

import os
import sys
import pickle
from typing import Dict, Tuple
from pathlib import Path

import numpy as np
import skimage.io
import open3d as o3d
from tqdm import tqdm
import os.path as osp

from cellpose import models, core, utils, plot
from cellpose.io import logger_setup
from typing import Any

# Set up CellPose logging
logger_setup()

# Check GPU availability
use_GPU = core.use_gpu()
print(f'>>> GPU activated? {int(use_GPU)}')


# Configuration constants
Z_STACK_START = 20  # Starting z-stack index (remove bad z-stacks)
Z_STACK_END = 250   # Ending z-stack index (for organoid; use 90:180 for monolayer)
Z_STACK_PROCESS = 230  # Number of z-stacks to process
CELLPOSE_DIAMETER = 80  # CellPose diameter in pixels
DBSCAN_EPS = 17  # DBSCAN epsilon parameter for clustering
DBSCAN_MIN_POINTS = 7  # DBSCAN minimum points parameter


def extract_frame_number(filepath: str) -> str:
    """
    Extract frame number from file path.
    
    Args:
        filepath: Path to image file.
        
    Returns:
        Frame number as string.
    """
    filename = osp.basename(filepath)
    # Extract frame number from filename (assumes format with 'msec' in name)
    parts = filename.split('_')
    for part in parts:
        if 'msec' in part:
            return part.split('msec')[0]
    # Fallback: try to extract from end of filename
    return filename.split('_')[-2].split('msec')[0] if 'msec' in filename else '0'


def process_cellpose_segmentation(
    imgs: np.ndarray,
    model: Any,  # CellPose model instance (varies by version)
    diameter: int
) -> np.ndarray:
    """
    Run CellPose segmentation on image stack.
    
    Args:
        imgs: Image array of shape (z, y, x).
        model: CellPose model instance.
        diameter: Cell diameter in pixels.
        
    Returns:
        Stitched masks array.
    """
    # z_axis=0 means the first dimension is the z-axis
    # channels parameter is deprecated in v4.0.1+, so we remove it
    out = model.eval(
        imgs[0],  # Input image of shape (z, y, x)
        diameter=diameter,
        do_3D=False,  # Use 2D algorithm with stitching
        stitch_threshold=0.5,
        z_axis=0  # Specify that first axis (axis 0) is the z-axis
    )
    masks_stitched = out[0]
    return masks_stitched


def extract_centers_and_masks(
    masks_stitched: np.ndarray,
    z_start: int = 0,
    z_end: int = Z_STACK_PROCESS
) -> Tuple[np.ndarray, Dict[Tuple[float, float], np.ndarray]]:
    """
    Extract cell centers and create center-to-mask mapping.
    
    Args:
        masks_stitched: Stitched masks array from CellPose.
        z_start: Starting z-stack index. Default is 0.
        z_end: Ending z-stack index (exclusive). Default is Z_STACK_PROCESS.
        
    Returns:
        Tuple of (centers, center_cellmask_map) where:
        - centers: Array of shape (n, 4) with columns [x, y, z, label]
        - center_cellmask_map: Dictionary mapping (cx, cy) to array of (x, y, z) points
    """
    center = []
    center_cellmask_map: Dict[Tuple[float, float], np.ndarray] = {}
    
    for i in tqdm(range(z_start, z_end), desc="Extracting centers"):
        if i >= masks_stitched.shape[0]:
            break
            
        labels = np.unique(masks_stitched[i])
        
        # Remove background label (0)
        labels = labels[labels != 0]

        for label in labels:
            y_coords, x_coords = np.where(masks_stitched[i] == label)
            
            if len(x_coords) == 0:
                continue
                
            # Calculate centroid
            cx = np.mean(x_coords).astype(np.float16)
            cy = np.mean(y_coords).astype(np.float16)
            center.append([cx, cy, i, label])

            # Store mask points for this center
            repeated_z = np.repeat(i, len(x_coords))
            center_cellmask_map[(cx, cy)] = np.column_stack([
                x_coords, y_coords, repeated_z
            ]).astype(np.int16)

    center = np.array(center).astype(np.float16)
    return center, center_cellmask_map


def cluster_centers(center: np.ndarray) -> np.ndarray:
    """
    Cluster cell centers using DBSCAN to group centers from same cell.
    
    Args:
        center: Array of shape (n, 4) with columns [x, y, z, label].
        
    Returns:
        Array of cluster labels for each center point.
    """
    if len(center) == 0:
        return np.array([])
    
    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(center[:, :3])

    with o3d.utility.VerbosityContextManager(
            o3d.utility.VerbosityLevel.Error) as cm:
        labels = np.array(
            pcd.cluster_dbscan(
                eps=DBSCAN_EPS, 
                min_points=DBSCAN_MIN_POINTS, 
                print_progress=False
            )
        )

    return labels


def main() -> None:
    """
    Main execution function for CellPose processing.
    
    Expects command line arguments:
        sys.argv[1]: data_dir - Directory containing input images
        sys.argv[2]: start_idx - Starting index for file processing
        sys.argv[3]: end_idx - Ending index for file processing (exclusive)
    """
    if len(sys.argv) < 4:
        raise ValueError(
            "Usage: python run_cellpose.py <data_dir> <start_idx> <end_idx>\n"
            "  data_dir: Directory containing input images\n"
            "  start_idx: Starting index for file processing\n"
            "  end_idx: Ending index for file processing (exclusive)"
        )
    
    data_dir = sys.argv[1]
    start_idx = int(sys.argv[2])
    end_idx = int(sys.argv[3])
    
    if not os.path.exists(data_dir):
        raise FileNotFoundError(f"Data directory not found: {data_dir}")

    # Set up save directory
    # Extract suffix from data_dir (assumes format: .../date/sample/region)
    path_parts = Path(data_dir).parts
    if len(path_parts) >= 4:
        save_dir_suffix = '/'.join(path_parts[-4:-1])
    else:
        save_dir_suffix = osp.basename(data_dir)
    
    # TODO: Make save_dir_root configurable via environment variable or config file
    save_dir_root = os.environ.get(
        'MITODEV_SAVE_ROOT',
        './data'
    )
    save_data_dir = osp.join(save_dir_root, save_dir_suffix, 'extractions')
    os.makedirs(save_data_dir, exist_ok=True)

    # Get image files
    files = sorted([
        f for f in os.listdir(data_dir) 
        if f.endswith('.tif') and '560nm' in f
    ])
    
    if start_idx >= len(files) or end_idx > len(files):
        raise IndexError(
            f"File indices [{start_idx}, {end_idx}) out of range "
            f"for {len(files)} files"
        )
    
    files = files[start_idx:end_idx]
    fpaths = [os.path.join(data_dir, f) for f in files]

    # Extract frame number and check if already processed
    frame_no = extract_frame_number(fpaths[0])
    line_centers_path = osp.join(save_data_dir, 'line_centers', f"frame_{frame_no}.npy")
    
    if osp.exists(line_centers_path):
        print(f"Frame {frame_no} already processed. Skipping...")
        sys.exit(0)

    # Load images
    print(f"Loading {len(fpaths)} image files...")
    imgs = [skimage.io.imread(f) for f in fpaths]
    imgs = np.array(imgs)
    
    # Remove bad z-stacks
    imgs = imgs[:, Z_STACK_START:Z_STACK_END]

    # Initialize CellPose model
    print("Initializing CellPose model...")
    # In CellPose 4.x, the model class is CellposeModel
    # Try CellposeModel first (CellPose 4.x), then fallback to Cellpose (older versions)
    if hasattr(models, 'CellposeModel'):
        model = models.CellposeModel(gpu=use_GPU, model_type='cyto2')
    elif hasattr(models, 'Cellpose'):
        model = models.Cellpose(gpu=use_GPU, model_type='cyto2')
    else:
        # Last resort: try to get model from models module
        model = models.CellposeModel(gpu=use_GPU, model_type='cyto2')

    # Run segmentation
    print("Running CellPose segmentation...")
    masks_stitched = process_cellpose_segmentation(
        imgs, model, CELLPOSE_DIAMETER
    )

    # Process each image in the batch
    for img_no in range(len(imgs)):
        print(f"Processing image {img_no + 1}/{len(imgs)}...")
        
        # Extract centers and create mapping
        center, center_cellmask_map = extract_centers_and_masks(
            masks_stitched, z_start=0, z_end=Z_STACK_PROCESS
        )

        if len(center) == 0:
            print(f"Warning: No centers found in image {img_no}")
            continue

        # Cluster centers to group centers from same cell
        print("Clustering centers...")
        labels = cluster_centers(center)
        
        # Update center labels with cluster assignments
        center[:, 3] = labels

        # Remove noise points (label -1)
        valid_mask = labels != -1
        center = center[valid_mask]

        # Create output directories
        os.makedirs(osp.join(save_data_dir, 'line_centers'), exist_ok=True)
        os.makedirs(osp.join(save_data_dir, 'cell_masks'), exist_ok=True)
        os.makedirs(osp.join(save_data_dir, 'center_mask_map'), exist_ok=True)

        # Save results
        np.save(
            osp.join(save_data_dir, 'line_centers', f"frame_{frame_no}.npy"), 
            center
        )
        np.save(
            osp.join(save_data_dir, 'cell_masks', f"frame_{frame_no}.npy"), 
            masks_stitched
        )

        # Save center-to-mask mapping
        center_mask_map_path = osp.join(
            save_data_dir, 'center_mask_map', 
            f'center_cellmask_map_{frame_no}.pkl'
        )
        with open(center_mask_map_path, 'wb') as f:
            pickle.dump(center_cellmask_map, f)
        
        print(f"Successfully saved results for frame {frame_no}")


if __name__ == '__main__':
    main()
