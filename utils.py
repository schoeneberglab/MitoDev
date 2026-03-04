"""
Utility functions for cell tracking and processing in MitoDev.

This module provides functions for:
- Cell filtering and tracking
- Point cloud processing and visualization
- Cell smoothing and filtering
- Configuration and file path management
- Cell data saving and loading
"""

import pickle
from typing import Dict, List, Tuple, Optional, Union, Any
from pathlib import Path

import numpy as np
import skimage.io
import yaml
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import os
from cellpose import plot
import os.path as osp
import pandas as pd
from scipy.ndimage import gaussian_filter
import open3d as o3d
from tqdm import tqdm


def visualize_line_centers(
    points_with_labels: np.ndarray, 
    colormap: str = 'tab10'
) -> None:
    """
    Visualize a point cloud with labels using Open3D.

    Args:
        points_with_labels: Array of shape (n, 4), where first 3 columns are 
                           x, y, z coordinates and last column is integer label.
        colormap: Matplotlib colormap name for label colors. Default is 'tab10'.
    
    Raises:
        AssertionError: If input array does not have shape (n, 4).
    """
    assert points_with_labels.shape[1] == 4, "Input array must have shape (n, 4)"

    # Separate points and labels
    points = points_with_labels[:, :3]
    labels = points_with_labels[:, 3].astype(int)

    # Unique labels and color mapping
    unique_labels = np.unique(labels)
    cmap = plt.get_cmap(colormap, len(unique_labels))
    label_to_color = {lab: cmap(i)[:3] for i, lab in enumerate(unique_labels)}

    # Assign colors to each point
    colors = np.array([label_to_color[lab] for lab in labels])

    # Create Open3D point cloud
    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(points)
    pcd.colors = o3d.utility.Vector3dVector(colors)

    # Visualize
    o3d.visualization.draw_geometries([pcd])


def filter_cells_by_z_length(centers: np.ndarray, min_z_stacks: int = 20) -> List[float]:
    """
    Filter cells that span at least a minimum number of z-stacks.

    Removes cells with labels -1.0 and 0.0 (background/noise labels) and cells
    that don't meet the minimum z-stack length requirement.

    Args:
        centers: Array of shape (n, 4) containing cell centers with columns 
                [x, y, z, label].
        min_z_stacks: Minimum number of z-stacks a cell must span to be kept.
                     Default is 20.

    Returns:
        List of cell labels that meet the filtering criteria.
    """
    label_list = np.unique(centers[:, -1])

    # Count occurrences of each label
    label_counts: Dict[float, int] = {}
    for label in label_list:
        label_counts[label] = 0

    for point in centers:
        label_counts[point[3]] += 1

    # Filter labels: keep those with enough z-stacks and exclude background labels
    keep = [
        label for label in label_list
        if label_counts[label] >= min_z_stacks and label not in (-1.0, 0.0)
    ]
    discard = [
        label for label in label_list
        if label_counts[label] < min_z_stacks or label in (-1.0, 0.0)
    ]

    print(f"Keeping {len(keep)} cells and discarding {len(discard)} cells")

    return keep


def make_xy_list(
    centers: np.ndarray, 
    label_list: List[float], 
    zstack: int, 
    window: int = 1
) -> Optional[np.ndarray]:
    """
    Select cell centers within a z-stack window and filter by label list.

    Args:
        centers: Array of shape (n, 4) containing cell centers with columns 
                [x, y, z, label].
        label_list: List of cell labels to include.
        zstack: Central z-stack index.
        window: Half-width of z-stack window. Centers within [zstack-window, 
               zstack+window] are included. Default is 1.

    Returns:
        Array of shape (m, 3) with columns [x, y, label] for selected centers,
        or None if no centers are found.
    """
    selected_centers = []
    z_min = zstack - window
    z_max = zstack + window + 1

    for point in centers:
        if (z_min <= point[2] < z_max) and (point[3] in label_list):
            selected_centers.append(point)

    if len(selected_centers) == 0:
        print(f'No centers found in z-stack {zstack} with window {window}')
        return None

    selected_centers = np.array(selected_centers)
    # Extract (x, y, label) from selected centers
    xy_list = np.concatenate([selected_centers[:, :2], selected_centers[:, 3:]], axis=1)
    
    return xy_list


def find_closest_cell(
    frame0_pts: np.ndarray, 
    frame1_pts: np.ndarray, 
    distance_threshold: float = 5.0
) -> np.ndarray:
    """
    Find closest matching cells between two frames based on Euclidean distance.

    For each cell in frame0, finds the closest cell in frame1. Only pairs
    with distance below the threshold are returned.

    Args:
        frame0_pts: Array of shape (n, 3) with columns [x, y, label] for frame 0.
        frame1_pts: Array of shape (m, 3) with columns [x, y, label] for frame 1.
        distance_threshold: Maximum distance for a valid cell pair. Default is 5.0.

    Returns:
        Array of shape (k, 2) with columns [label_frame0, label_frame1] for 
        matched cell pairs.
    """
    label_pairs = []

    unique_labels_frame0 = np.unique(frame0_pts[:, -1])
    unique_labels_frame1 = np.unique(frame1_pts[:, -1])

    for label_0 in unique_labels_frame0:
        min_dist = float('inf')
        label_assignment = -1

        frame0_pts_filtered = frame0_pts[frame0_pts[:, -1] == label_0]
        
        if len(frame0_pts_filtered) == 0:
            continue

        for label_1 in unique_labels_frame1:
            frame1_pts_filtered = frame1_pts[frame1_pts[:, -1] == label_1]

            if len(frame1_pts_filtered) == 0:
                continue

            # Compute pairwise distances between all points
            distances = []
            for pt0 in frame0_pts_filtered:
                x0, y0, _ = pt0
                for pt1 in frame1_pts_filtered:
                    x1, y1, _ = pt1
                    dist = np.sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2)
                    distances.append(dist)

            closest_dist = np.min(distances)

            if closest_dist < min_dist:
                min_dist = closest_dist
                label_assignment = label_1

        if min_dist < distance_threshold:
            label_pairs.append([label_0, label_assignment, min_dist])

    if len(label_pairs) == 0:
        return np.array([]).reshape(0, 2)

    label_pairs = np.array(label_pairs)
    return label_pairs[:, :2]  # Return only [label_0, label_1]


def change_labels_f1(
    pairs: np.ndarray, 
    centers0: np.ndarray, 
    centers1: np.ndarray
) -> Tuple[List[List[float]], List[List[float]]]:
    """
    Update cell labels in frame 1 to match frame 0 based on tracking pairs.

    This function propagates labels from frame 0 to frame 1 for tracked cells,
    ensuring consistent labeling across frames.

    Args:
        pairs: Array of shape (n, 2) with columns [label_frame0, label_frame1] 
              representing tracked cell pairs.
        centers0: Array of shape (m, 4) with columns [x, y, z, label] for frame 0.
        centers1: Array of shape (k, 4) with columns [x, y, z, label] for frame 1.

    Returns:
        Tuple of (new_frame0_pts, new_frame1_pts) where both are lists of 
        [x, y, z, label] arrays with updated labels.
    """
    new_frame0_pts = []
    new_frame1_pts = []

    for pair in pairs:
        label_frame0, label_frame1 = pair[0], pair[1]
        
        # Add cells from frame 0 that have been tracked
        for point in centers0:
            if point[3] == label_frame0:
                new_frame0_pts.append([point[0], point[1], point[2], label_frame0])
        
        # Add cells from frame 1 with labels updated to match frame 0
        for point in centers1:
            if point[3] == label_frame1:
                new_frame1_pts.append([point[0], point[1], point[2], label_frame0])

    return new_frame0_pts, new_frame1_pts


def generate_distinct_colors(n: int) -> np.ndarray:
    """
    Generate n visually distinct colors using KMeans clustering.

    Args:
        n: Number of distinct colors to generate.

    Returns:
        Array of shape (n, 3) with RGB values in range [0, 1].
    """
    # Generate random RGB values
    rgb_values = np.random.random((n, 3))

    # Apply KMeans clustering to find n distinct colors
    kmeans = KMeans(n_clusters=n, random_state=0, n_init=10).fit(rgb_values)
    distinct_colors = kmeans.cluster_centers_

    return distinct_colors


def vis_cell_mask(imgs, frame0_linecenters, frame1_linecenters, pairs_of_interest, z_stacks):
    zstack = z_stacks

    img_f0 = imgs[0]
    img_f1 = imgs[1]

    img0 = plot.image_to_rgb(img_f0[zstack, :, :].copy(), channels=[0])

    for pair in pairs_of_interest:
        print("distance:", pair[2])
        plt.subplot(121)
        pair0 = pair[0]
        # crop image data to where there's signal
        start_idx = -1
        end_idx = -1
        for idx in range(img0.shape[1] - 1):
            delta = img0[1000, idx] - img0[1000, idx + 1]
            if abs(delta[0]) > 0:
                start_idx = idx
                break
        for idx in range(img0.shape[1] - 2, -1, -1):
            delta = img0[1000, idx] - img0[1000, idx + 1]
            if abs(delta[0]) > 0:
                end_idx = idx
                break

        # collect labels
        to_plot = []
        pt_labels = []
        for val in frame0_linecenters:
            if val[2] == zstack - 20:
                to_plot.append([val[0], val[1]])
                pt_labels.append(val[3])
        to_plot_df = pd.DataFrame(to_plot)
        to_plot_df.index = pt_labels

        # plot
        imgout = img0.copy()[:, start_idx:end_idx]  # [vertical, horizontal]
        plt.imshow(imgout, origin="lower")
        # plt.gca().invert_yaxis()

        plt.scatter(to_plot_df[1], to_plot_df[0], marker='.', c='gray')
        for i, txt in enumerate(to_plot_df.index):
            if str(txt) == str(pair0):
                plt.annotate(txt, (to_plot[i][1], to_plot[i][0]), fontsize=4, c='red')
            else:
                plt.annotate(txt, (to_plot[i][1], to_plot[i][0]), fontsize=4, c='green')

        # plt.show()

        ############## frame 1
        pair1 = pair[1]
        plt.subplot(122)
        img1 = plot.image_to_rgb(img_f1[zstack, :, :].copy(), channels=[0])

        # crop image data to where there's signal
        start_idx = -1
        end_idx = -1
        for idx in range(img1.shape[1] - 1):
            delta = img1[1000, idx] - img1[1000, idx + 1]
            if abs(delta[0]) > 0:
                start_idx = idx
                break
        for idx in range(img1.shape[1] - 2, -1, -1):
            delta = img1[1000, idx] - img1[1000, idx + 1]
            if abs(delta[0]) > 0:
                end_idx = idx
                break

        # collect labels
        to_plot = []
        pt_labels = []
        for val in frame1_linecenters:
            if val[2] == zstack - 20:
                to_plot.append([val[0], val[1]])
                pt_labels.append(val[3])
        to_plot_df = pd.DataFrame(to_plot)
        to_plot_df.index = pt_labels

        # plot
        imgout = img1.copy()[:, start_idx:end_idx]  # [vertical, horizontal]

        plt.imshow(imgout, origin="lower")
        # plt.gca().invert_yaxis()

        plt.scatter(to_plot_df[1], to_plot_df[0], marker='.', c='gray')
        for i, txt in enumerate(to_plot_df.index):
            if str(txt) == str(pair1):
                plt.annotate(txt, (to_plot[i][1], to_plot[i][0]), fontsize=4, c='red')
            else:
                plt.annotate(txt, (to_plot[i][1], to_plot[i][0]), fontsize=4, c='green')

        plt.show()


def find_faulty_assignments(assignments: Union[set, list]) -> List[Tuple[float, float]]:
    """
    Find conflicting cell label assignments.

    Detects cases where multiple cells from frame 0 are assigned to the same
    cell in frame 1, which indicates tracking conflicts.

    Args:
        assignments: Set or list of tuples (label_frame0, label_frame1) 
                    representing cell assignments.

    Returns:
        List of tuples (label_frame0, assigned_label) representing conflicts.
    """
    assigned_labels = set()
    conflicts = []

    for label_0, assigned_label in assignments:
        if assigned_label in assigned_labels:
            conflicts.append((label_0, assigned_label))
        else:
            assigned_labels.add(assigned_label)

    return conflicts


def make_cell(cell_pts: np.ndarray) -> Dict[int, np.ndarray]:
    """
    Organize cell points by z-stack into a dictionary.

    Args:
        cell_pts: Array of shape (n, 3) with columns [x, y, z].

    Returns:
        Dictionary mapping z-stack indices to arrays of shape (m, 2) with 
        columns [x, y] for points at that z-stack.
    """
    min_z = int(np.min(cell_pts[:, 2]))
    max_z = int(np.max(cell_pts[:, 2]))
    cell: Dict[int, np.ndarray] = {}
    
    for z in range(min_z, max_z + 1):
        z_points = cell_pts[cell_pts[:, 2] == z]
        if len(z_points) > 0:
            cell[z] = z_points[:, :2]

    return cell


def cell_to_cellpoints(cell: Dict[int, np.ndarray]) -> np.ndarray:
    """
    Convert cell dictionary back to point array format.

    Args:
        cell: Dictionary mapping z-stack indices to arrays of shape (m, 2) 
              with columns [x, y].

    Returns:
        Array of shape (n, 3) with columns [x, y, z].
    """
    cell_points = []
    for z in sorted(cell.keys()):
        xy_points = cell[z]
        for xy in xy_points:
            cell_points.append([xy[0], xy[1], z])
    
    if len(cell_points) == 0:
        return np.array([]).reshape(0, 3)
    
    cell_points = np.vstack(cell_points)
    return cell_points


def smooth_cells(cell_points: np.ndarray, window: int = 5) -> np.ndarray:
    """
    Smooth cell points across z-stacks by filling gaps.

    For empty z-stacks, copies points from the previous z-stack. This helps
    maintain continuity in cell tracking across z-stacks.

    Args:
        cell_points: Array of shape (n, 3) with columns [x, y, z].
        window: Window size for area calculation (currently unused but kept for 
               future use). Default is 5.

    Returns:
        Array of shape (m, 3) with columns [x, y, z] after smoothing.
    """
    if len(cell_points) == 0:
        return cell_points

    # Sort the cell points by z
    cell_points = cell_points[cell_points[:, 2].argsort()]
    min_z = int(np.min(cell_points[:, 2]))
    max_z = int(np.max(cell_points[:, 2]))

    unsmooth_cell = make_cell(cell_points)

    for z in range(min_z, max_z + 1):
        if z not in unsmooth_cell or len(unsmooth_cell[z]) == 0:
            # Fill empty z-stacks with previous z-stack data
            if z - 1 >= min_z and (z - 1) in unsmooth_cell:
                unsmooth_cell[z] = unsmooth_cell[z - 1].copy()
            continue

        # Calculate average area in the previous window z-stacks
        areas = []
        z_prev = z - 1
        while z_prev >= min_z and z_prev >= z - window:
            if z_prev in unsmooth_cell:
                area_z = len(unsmooth_cell[z_prev])
                areas.append(area_z)
            z_prev -= 1

        # Note: Area-based filtering is currently commented out but could be
        # re-enabled in the future for more sophisticated smoothing

    cell_points = cell_to_cellpoints(unsmooth_cell)
    return cell_points


def apply_gaussian_filter(
    cell_points: np.ndarray, 
    sigma: float = 10.0, 
    threshold: float = 0.5
) -> np.ndarray:
    """
    Apply Gaussian smoothing filter to cell points.

    Converts cell points to a 3D binary mask, applies Gaussian filtering,
    and extracts smoothed points above a threshold.

    Args:
        cell_points: Array of shape (n, 3) with columns [x, y, z].
        sigma: Standard deviation for Gaussian filter. Default is 10.0.
        threshold: Threshold value for extracting smoothed points. Points with
                  filter values above this threshold are kept. Default is 0.5.

    Returns:
        Array of shape (m, 3) with columns [x, y, z] after Gaussian smoothing.
    """
    if len(cell_points) == 0:
        return cell_points

    # Compute the min and max coordinates along each axis
    min_x, max_x = cell_points[:, 0].min(), cell_points[:, 0].max()
    min_y, max_y = cell_points[:, 1].min(), cell_points[:, 1].max()
    min_z, max_z = cell_points[:, 2].min(), cell_points[:, 2].max()

    # Compute the range for each axis
    x_range = int(max_x - min_x) + 1
    y_range = int(max_y - min_y) + 1
    z_range = int(max_z - min_z) + 1

    # Initialize the cell mask and mark the cell points
    cell_mask = np.zeros((x_range, y_range, z_range), dtype=np.float32)
    cell_indices = np.rint(cell_points - [min_x, min_y, min_z]).astype(int)
    
    # Ensure indices are within bounds
    valid_mask = (
        (cell_indices[:, 0] >= 0) & (cell_indices[:, 0] < x_range) &
        (cell_indices[:, 1] >= 0) & (cell_indices[:, 1] < y_range) &
        (cell_indices[:, 2] >= 0) & (cell_indices[:, 2] < z_range)
    )
    valid_indices = cell_indices[valid_mask]
    
    if len(valid_indices) > 0:
        cell_mask[valid_indices[:, 0], valid_indices[:, 1], valid_indices[:, 2]] = 1.0

    # Apply Gaussian filter
    cell_mask = gaussian_filter(cell_mask, sigma=sigma)

    # Vectorized approach to find points where the mask exceeds the threshold
    smoothed_indices = np.argwhere(cell_mask > threshold)

    # Convert the indices back to original coordinate space
    smoothed_cell_points = smoothed_indices + [min_x, min_y, min_z]

    return smoothed_cell_points


def vis_cells(
    cell_points: np.ndarray, 
    save: bool = False, 
    name: str = 'cell_points',
    save_path: Optional[str] = None
) -> None:
    """
    Visualize or save cell points as a 3D point cloud.

    Args:
        cell_points: Array of shape (n, 3) with columns [x, y, z].
        save: If True, save the point cloud to file instead of displaying.
              Default is False.
        name: Base name for saved file (used if save_path is None).
              Default is 'cell_points'.
        save_path: Full path for saving point cloud. If None and save=True,
                  uses default path. Default is None.
    """
    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(cell_points)

    if not save:
        o3d.visualization.draw_geometries(
            [pcd],
            window_name='3D Points Visualization',
            width=800, 
            height=600
        )
    else:
        if save_path is None:
            save_path = f"/home/schoeneberglab/Desktop/{name}.pcd"
        o3d.io.write_point_cloud(save_path, pcd)


def load_config(config_path: Union[str, Path]) -> Dict[str, Any]:
    """
    Load configuration from a YAML file.

    Args:
        config_path: Path to the YAML configuration file.

    Returns:
        Dictionary containing configuration parameters.

    Raises:
        FileNotFoundError: If the config file doesn't exist.
        yaml.YAMLError: If the YAML file is malformed.
    """
    config_path = str(config_path)
    
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    
    with open(config_path, 'r') as file:
        try:
            cfg = yaml.safe_load(file)
            if cfg is None:
                raise ValueError("Configuration file is empty")
            return cfg
        except yaml.YAMLError as exc:
            raise yaml.YAMLError(f"Error parsing YAML file {config_path}: {exc}") from exc


def get_paths(
    data_dir: str, 
    metadata_dir: str, 
    data_folder: str, 
    region: str
) -> Tuple[List[str], List[str], List[str], List[str]]:
    """
    Get file paths for line centers, masks, mappings, and mitochondrial images.

    Args:
        data_dir: Root directory containing data.
        metadata_dir: Directory containing metadata (line centers, masks, etc.).
        data_folder: Name of the data folder.
        region: Name of the region.

    Returns:
        Tuple of (line_centers_paths, mito_files_paths, masks_paths, 
                 center_cellmask_maps_paths), each a sorted list of file paths.

    Raises:
        FileNotFoundError: If required directories don't exist.
    """
    line_centers_dir = osp.join(metadata_dir, 'line_centers')
    if not os.path.exists(line_centers_dir):
        raise FileNotFoundError(f"Line centers directory not found: {line_centers_dir}")
    all_sample_path = sorted([
        osp.join(line_centers_dir, f) 
        for f in os.listdir(line_centers_dir)
        if os.path.isfile(osp.join(line_centers_dir, f))
    ])

    cell_masks_dir = osp.join(metadata_dir, 'cell_masks')
    if not os.path.exists(cell_masks_dir):
        raise FileNotFoundError(f"Cell masks directory not found: {cell_masks_dir}")
    all_masks = sorted([
        osp.join(cell_masks_dir, f) 
        for f in os.listdir(cell_masks_dir)
        if os.path.isfile(osp.join(cell_masks_dir, f))
    ])

    center_cellmask_mapping_dir = osp.join(metadata_dir, 'center_mask_map')
    if not os.path.exists(center_cellmask_mapping_dir):
        raise FileNotFoundError(f"Center mask mapping directory not found: {center_cellmask_mapping_dir}")
    center_cellmask_maps = sorted([
        osp.join(center_cellmask_mapping_dir, f) 
        for f in os.listdir(center_cellmask_mapping_dir)
        if os.path.isfile(osp.join(center_cellmask_mapping_dir, f))
    ])

    mito_img_dir = osp.join(data_dir, data_folder, region, "Processed_Data")
    if not os.path.exists(mito_img_dir):
        raise FileNotFoundError(f"Mitochondrial image directory not found: {mito_img_dir}")
    all_sample_path_mito = sorted([
        osp.join(mito_img_dir, f) 
        for f in os.listdir(mito_img_dir) 
        if '488nm' in f and 'processed' in f and os.path.isfile(osp.join(mito_img_dir, f))
    ])

    return all_sample_path, all_sample_path_mito, all_masks, center_cellmask_maps


def load_center_cell_masks(
    center_cell_masks_fpaths: List[str], 
    frame_range: List[int]
) -> List[Dict[Tuple[float, float], np.ndarray]]:
    """
    Load center-to-cell-mask mappings for specified frame range.

    Args:
        center_cell_masks_fpaths: List of file paths to pickle files containing
                                 center-to-cell-mask mappings.
        frame_range: List of [start_frame, end_frame] indices (end_frame exclusive).

    Returns:
        List of dictionaries mapping (x, y) center coordinates to arrays of
        shape (n, 3) with columns [x, y, z] for cell mask points.

    Raises:
        IndexError: If frame_range is out of bounds.
        FileNotFoundError: If any mask file doesn't exist.
    """
    start_frame, end_frame = frame_range[0], frame_range[1]
    
    if end_frame > len(center_cell_masks_fpaths):
        raise IndexError(
            f"Frame range [{start_frame}, {end_frame}) exceeds available "
            f"frames ({len(center_cell_masks_fpaths)})"
        )
    
    center_cellmask_maps = []
    print(f"Loading cell masks for frames {start_frame} to {end_frame - 1}: ")
    
    for idx, path in tqdm(
        enumerate(center_cell_masks_fpaths[start_frame:end_frame]),
        total=end_frame - start_frame
    ):
        if not os.path.exists(path):
            raise FileNotFoundError(f"Mask file not found: {path}")
        
        with open(path, 'rb') as f:
            center_cellmask_maps.append(pickle.load(f))
    
    return center_cellmask_maps


def get_cell_bounds(cell_list: List[np.ndarray]) -> List[float]:
    """
    Compute bounding box for a list of cell point arrays.

    Args:
        cell_list: List of arrays, each of shape (n, 3) with columns [x, y, z].

    Returns:
        List of [min_x, max_x, min_y, max_y, min_z, max_z] bounding box coordinates.

    Raises:
        ValueError: If cell_list is empty.
    """
    if len(cell_list) == 0:
        raise ValueError("cell_list cannot be empty")

    min_x, min_y, min_z = float('inf'), float('inf'), float('inf')
    max_x, max_y, max_z = float('-inf'), float('-inf'), float('-inf')
    
    for cell in cell_list:
        if len(cell) == 0:
            continue
        x, y, z = cell[:, 0], cell[:, 1], cell[:, 2]
        min_x = min(min_x, x.min())
        min_y = min(min_y, y.min())
        min_z = min(min_z, z.min())
        max_x = max(max_x, x.max())
        max_y = max(max_y, y.max())
        max_z = max(max_z, z.max())

    bounds = [min_x, max_x, min_y, max_y, min_z, max_z]
    return bounds


def save_cells(
    cell_dict: Dict[float, Dict[str, Any]], 
    all_sample_path_mito: List[str], 
    data_folder: str, 
    region: str,
    save_root_dir: str, 
    frame_range: List[int],
    z_offset: int = 20
) -> None:
    """
    Save cell volumes extracted from mitochondrial images.

    For each cell and frame, extracts the cell volume from the mitochondrial
    image, normalizes it, and saves as a TIFF file.

    Args:
        cell_dict: Dictionary mapping cell labels to dictionaries with keys:
                  'cell' (list of point arrays) and 'bounds' (bounding box).
        all_sample_path_mito: List of file paths to mitochondrial images.
        data_folder: Name of the data folder.
        region: Name of the region.
        save_root_dir: Root directory for saving cell volumes.
        frame_range: List of [start_frame, end_frame] indices (end_frame exclusive).
        z_offset: Z-stack offset for accessing mitochondrial images. Default is 20.

    Raises:
        IndexError: If frame indices are out of bounds.
        FileNotFoundError: If mitochondrial image files don't exist.
    """
    start_frame, end_frame = frame_range[0], frame_range[1]
    
    if end_frame > len(all_sample_path_mito):
        raise IndexError(
            f"Frame range [{start_frame}, {end_frame}) exceeds available "
            f"mitochondrial images ({len(all_sample_path_mito)})"
        )

    for frame in range(start_frame, end_frame):
        # Load mitochondrial image file
        img_mito_file = all_sample_path_mito[frame]
        if not os.path.exists(img_mito_file):
            raise FileNotFoundError(f"Mitochondrial image not found: {img_mito_file}")
        
        print(f'Loading mitochondrial file for frame {frame}: {img_mito_file}')
        img_mito = skimage.io.imread(img_mito_file)

        for cell_label, cell_data in cell_dict.items():
            cell_over_time = cell_data['cell']  # List of cells over time
            boundaries = cell_data['bounds']

            min_x, max_x, min_y, max_y, min_z, max_z = boundaries

            # Create cell volume array
            z_size = int(max_z - min_z) + 1
            x_size = int(max_x - min_x) + 1
            y_size = int(max_y - min_y) + 1
            # Use (z, y, x) ordering so indexing matches img_mito[z, y, x]
            cell_volume = np.zeros((z_size, y_size, x_size), dtype=np.float16)
            
            cell_points = cell_over_time[frame]

            # Fill cell volume with mitochondrial image values
            for point_idx, point in enumerate(cell_points):
                xi, yi, zi = point[0], point[1], point[2]
                
                # Calculate indices
                z_idx = int(zi - min_z)
                x_idx = int(xi - min_x)
                y_idx = int(yi - min_y)
                img_z_idx = int(zi + z_offset)
                
                # Check bounds
                if (0 <= z_idx < z_size and 
                    0 <= x_idx < x_size and 
                    0 <= y_idx < y_size and
                    img_z_idx < img_mito.shape[0] and
                    int(xi) < img_mito.shape[1] and
                    int(yi) < img_mito.shape[2]):
                    try:
                        cell_volume[z_idx, y_idx, x_idx] = img_mito[
                            img_z_idx, int(yi), int(xi)
                        ]
                    except (IndexError, ValueError) as e:
                        print(f'Warning: Error accessing point {point_idx} in frame {frame}: {e}')
                        continue

            # Normalize cell volume
            sc = cell_volume.copy()
            # Find the 99.95th percentile for normalization
            max_val = np.percentile(sc.astype(np.float32), 99.95)
            if max_val > 0:
                sc = np.clip(sc, 0, max_val)
                sc = (sc / max_val) * 255.0
            else:
                sc = np.zeros_like(sc)
            
            # Save the cell volume
            save_dir = osp.join(
                save_root_dir, 'tracked_cells_smooth',
                f'cell_{cell_label}', 'mitograph', f'frame_{frame}'
            )
  
            # Convert to uint8 first as you intended
            image_to_save = sc.astype(np.uint8)

            # Check contrast: returns True if the image is low contrast
            if not skimage.exposure.is_low_contrast(image_to_save):
                os.makedirs(save_dir, exist_ok=True)
                save_path = osp.join(save_dir, f'frame_{frame}.tif')
                skimage.io.imsave(save_path, image_to_save)
            else:
                print(f"Skipping {cell_label}: Low contrast detected.")


# Backward compatibility: alias for renamed function
filter = filter_cells_by_z_length
