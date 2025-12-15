"""
Cell tracking functions for MitoDev.

This module provides functions for:
- Tracking cells across time frames
- Filling cell data across frames with smoothing
"""

import pickle
import os
from typing import Dict, List, Tuple, Any, Optional
import numpy as np
from tqdm import tqdm
import os.path as osp

from utils import (
    make_xy_list, 
    find_closest_cell, 
    find_faulty_assignments, 
    change_labels_f1, 
    filter_cells_by_z_length, 
    load_center_cell_masks, 
    smooth_cells, 
    apply_gaussian_filter, 
    get_cell_bounds
)


def track(
    line_centers_path: List[str], 
    num_frames: int, 
    z_stacks: List[int], 
    min_z_len: int
) -> Tuple[List[np.ndarray], np.ndarray]:
    """
    Track cells across multiple time frames.

    For each pair of consecutive frames, matches cells based on spatial proximity
    across multiple z-stacks. Resolves conflicts where multiple cells are assigned
    to the same target.

    Args:
        line_centers_path: List of file paths to line center arrays for each frame.
        num_frames: Total number of frames to process.
        z_stacks: List of z-stack indices to use for tracking.
        min_z_len: Minimum number of z-stacks a cell must span to be tracked.

    Returns:
        Tuple of (cell_pts, final_tracked_labels) where:
        - cell_pts: List of arrays, each of shape (n, 4) with columns [x, y, z, label]
                  for each frame.
        - final_tracked_labels: Array of unique cell labels that were successfully
                               tracked across all frames.

    Raises:
        FileNotFoundError: If any line center file doesn't exist.
        ValueError: If num_frames doesn't match the number of line_centers_path files.
    """
    if len(line_centers_path) < num_frames:
        raise ValueError(
            f"Number of line center paths ({len(line_centers_path)}) is less than "
            f"num_frames ({num_frames})"
        )

    cell_pts: List[np.ndarray] = []
    
    # Load first frame
    if not os.path.exists(line_centers_path[0]):
        raise FileNotFoundError(f"Line center file not found: {line_centers_path[0]}")
    line_centers_t0 = np.load(line_centers_path[0])

    # Handle single frame case
    if num_frames == 1:
        cell_pts.append(line_centers_t0)
        return cell_pts, np.unique(line_centers_t0[:, -1])

    # Track across frames
    for time_f in range(1, num_frames):
        # Filter cells by minimum z-stack length
        keep_t0 = filter_cells_by_z_length(line_centers_t0, min_z_stacks=min_z_len)
        print(f'Tracking across frames {time_f - 1} and {time_f}')

        # Load and filter next frame
        if not os.path.exists(line_centers_path[time_f]):
            raise FileNotFoundError(f"Line center file not found: {line_centers_path[time_f]}")
        print(f'Processing frame {time_f}')
        line_centers_t1 = np.load(line_centers_path[time_f])
        keep_t1 = filter_cells_by_z_length(line_centers_t1, min_z_stacks=min_z_len)

        # Find cell pairs across all z-stacks
        all_pairs = set()
        for zstack in z_stacks:
            frame_t0_z = make_xy_list(line_centers_t0, keep_t0, zstack, window=2)
            frame_t1_z = make_xy_list(line_centers_t1, keep_t1, zstack, window=2)

            if frame_t0_z is None or frame_t1_z is None:
                continue

            pairs_z = find_closest_cell(frame_t0_z, frame_t1_z)
            pairs_z_set = set(tuple(pair) for pair in pairs_z)
            all_pairs = all_pairs.union(pairs_z_set)

        # Resolve conflicts (multiple cells assigned to same target)
        conflicts = find_faulty_assignments(all_pairs)
        if len(conflicts) > 0:
            conflict_set = set(conflicts)
            all_pairs = [pair for pair in all_pairs if pair not in conflict_set]
        
        print(f'{len(all_pairs)} cells tracked across frames {time_f - 1} and {time_f}')

        # Update labels in frame 1 to match frame 0
        all_pairs_array = np.array(list(all_pairs)) if len(all_pairs) > 0 else np.array([]).reshape(0, 2)
        new_frame_t0_pts, new_frame_t1_pts = change_labels_f1(
            all_pairs_array, line_centers_t0, line_centers_t1
        )

        # Frame t1 becomes frame t0 for next iteration
        line_centers_t0 = np.array(new_frame_t1_pts)

        # Save tracked points for frame t0
        cell_pts.append(np.array(new_frame_t0_pts).astype(np.float16))

    # Add final frame
    cell_pts.append(np.array(new_frame_t1_pts).astype(np.float16))

    # Filter to only keep cells tracked across all frames
    final_tracked_labels = np.unique(cell_pts[-1][:, -1])

    for i in range(len(cell_pts)):
        cell_pts[i] = cell_pts[i][np.isin(cell_pts[i][:, -1], final_tracked_labels)]
        unique_labels = np.unique(cell_pts[i][:, -1])
        assert len(unique_labels) == len(final_tracked_labels), (
            f"Frame {i}: Expected {len(final_tracked_labels)} labels, "
            f"got {len(unique_labels)}"
        )

    return cell_pts, final_tracked_labels


def fill_cells_across_frames(
    cell_pts: List[np.ndarray], 
    final_tracked_labels: np.ndarray, 
    frame_range: List[int], 
    center_cellmask_maps: List[str], 
    save_root_dir: str,
    data_folder: str, 
    region: str,
    z_window: int = 5, 
    sigma: float = 10.0, 
    resume: bool = False,
    minimal: bool = False,
) -> Dict[float, Dict[str, Any]]:
    """
    Fill cell volumes across frames with smoothing and save incrementally.

    For each tracked cell, extracts cell points from center-to-mask mappings,
    applies smoothing filters, and saves the results incrementally to disk.

    Args:
        cell_pts: List of arrays, each of shape (n, 4) with columns [x, y, z, label]
                 for each frame.
        final_tracked_labels: Array of unique cell labels to process.
        frame_range: List of [start_frame, end_frame] indices (end_frame exclusive).
        center_cellmask_maps: List of file paths to center-to-mask mapping dictionaries.
        save_root_dir: Root directory for saving cell dictionaries.
        data_folder: Name of the data folder.
        region: Name of the region.
        z_window: Window size for z-stack smoothing. Default is 5.
        sigma: Standard deviation for Gaussian smoothing. Default is 10.0.
        resume: If True, load existing cell_dict and continue. Default is False.
        minimal: If True, only save 10 cells. Default is False.
    Returns:
        Dictionary mapping cell labels to dictionaries with keys:
        - 'cell': List of point arrays (one per frame) of shape (n, 3) with [x, y, z]
        - 'bounds': List of [min_x, max_x, min_y, max_y, min_z, max_z] bounding box

    Raises:
        FileNotFoundError: If center_cellmask_maps files don't exist.
        KeyError: If center coordinates are not found in mapping dictionaries.
    """
    save_dir = save_root_dir
    cell_dict_path = osp.join(save_dir, 'cell_dict.pkl')
    cell_dict_temp_path = osp.join(save_dir, 'cell_dict_temp.pkl')

    if resume:
        if os.path.exists(cell_dict_path):
            with open(cell_dict_path, 'rb') as f:
                cell_dict_old = pickle.load(f)
            cell_dict = {}
        else:
            print(f"Warning: Resume flag set but cell_dict.pkl not found. Starting fresh.")
            cell_dict = {}
    else:
        cell_dict: Dict[float, Dict[str, Any]] = {}

    # Load center-to-cell-mask mappings for all frames
    center_cellmask_maps_loaded = load_center_cell_masks(center_cellmask_maps, frame_range)

    start_frame, end_frame = frame_range[0], frame_range[1]

    if minimal:
        print(f"Saving only 10 cells since minimal flag is set")
        final_tracked_labels = final_tracked_labels[:10]

    for lbl_idx, cell_lbl in enumerate(tqdm(final_tracked_labels, desc="Processing cells")):
        cell_dict[cell_lbl] = {'cell': [], 'bounds': []}

        for frame in tqdm(
            range(start_frame, end_frame),
            desc=f'Cell Label - {cell_lbl} ({lbl_idx + 1}/{len(final_tracked_labels)})'
        ):
            # Get line centers for this cell in this frame
            idx = np.where(cell_pts[frame][:, -1] == cell_lbl)[0]
            if len(idx) == 0:
                # No points for this cell in this frame
                cell_dict[cell_lbl]['cell'].append(np.array([]).reshape(0, 3).astype(np.uint16))
                continue

            line_centers_lb = cell_pts[frame][idx]
            cell_points = []

            # Extract cell points from center-to-mask mappings
            for line_center in line_centers_lb:
                cx, cy, z, lab = line_center
                center_key = (float(cx), float(cy))
                
                if center_key not in center_cellmask_maps_loaded[frame - start_frame]:
                    print(f"Warning: Center ({cx}, {cy}) not found in frame {frame} mapping")
                    continue
                
                cell_points_for_z = center_cellmask_maps_loaded[frame - start_frame][center_key]
                cell_points.extend(cell_points_for_z)

            if len(cell_points) == 0:
                cell_dict[cell_lbl]['cell'].append(np.array([]).reshape(0, 3).astype(np.uint16))
                continue

            # Apply smoothing
            cell_points = np.array(cell_points)
            cell_points = smooth_cells(cell_points, window=z_window)
            cell_points = apply_gaussian_filter(cell_points, sigma=sigma)

            cell_dict[cell_lbl]['cell'].append(np.array(cell_points).astype(np.uint16))

        # Compute bounding box for this cell
        cell_bounds = get_cell_bounds(cell_dict[cell_lbl]['cell'])
        cell_dict[cell_lbl]['bounds'] = cell_bounds

        # Save incrementally after each cell
        os.makedirs(save_dir, exist_ok=True)
        print(f"Saving cell_dict after processing cell {cell_lbl}...")
        with open(cell_dict_temp_path, 'wb') as f:
            pickle.dump(cell_dict, f)

        # Atomic rename for safe saving
        os.rename(cell_dict_temp_path, cell_dict_path)

    return cell_dict
