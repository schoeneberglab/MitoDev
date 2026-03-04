import os
import os.path as osp
from glob import glob
import argparse
import base64

import tifffile
import imageio
import numpy as np
from tqdm import tqdm


# --------------------------------------------------
# CLI
# --------------------------------------------------
def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate embedded HTML with MIP videos for tracked cells"
    )
    parser.add_argument(
        "--root",
        required=True,
        help="Root directory containing tracked_cells_smooth/",
    )
    parser.add_argument("--fps", type=int, default=10)
    parser.add_argument("--out", default="cell_gallery")
    return parser.parse_args()


# --------------------------------------------------
# Data loading
# --------------------------------------------------
def load_cell_frames(cell_dir):
    """
    Expected:
    cell_dir/
        mitograph/
            frame_0/frame_0.tif
            frame_1/frame_1.tif
            ...
    """
    mito_dir = osp.join(cell_dir, "mitograph")
    frame_dirs = sorted(glob(osp.join(mito_dir, "frame_*")))

    imgs = []
    for fd in frame_dirs:
        frame_id = osp.basename(fd)
        tif_path = osp.join(fd, f"{frame_id}.tif")
        if osp.exists(tif_path):
            imgs.append(tifffile.imread(tif_path))

    return imgs


def make_mip_video(imgs):
    """Max intensity projection over Z"""
    return np.stack([np.max(img, axis=0) for img in imgs], axis=0)


# --------------------------------------------------
# Video writing
# --------------------------------------------------
def normalize_to_uint8(video):
    vmin = video.min()
    vmax = video.max()
    if vmax > vmin:
        video = (video - vmin) / (vmax - vmin)
    else:
        video = np.zeros_like(video)
    video_u8 = (video * 255).astype(np.uint8)

    # Pad to even dimensions (libx264 requirement)
    h, w = video_u8.shape[1], video_u8.shape[2]
    h_padded = h if h % 2 == 0 else h + 1
    w_padded = w if w % 2 == 0 else w + 1
    
    if h_padded != h or w_padded != w:
        video_padded = np.zeros((video_u8.shape[0], h_padded, w_padded), dtype=np.uint8)
        video_padded[:, :h, :w] = video_u8
        return video_padded
    
    return video_u8


def save_mp4(video, path, fps):
    video_u8 = normalize_to_uint8(video)

    writer = imageio.get_writer(
        path,
        fps=fps,
        format="ffmpeg",
        codec="libx264",
        pixelformat="yuv420p",
        macro_block_size=1,  # Allow any frame size
    )

    for frame in video_u8:
        writer.append_data(frame)

    writer.close()


# --------------------------------------------------
# HTML embedding
# --------------------------------------------------
def file_to_base64(path):
    with open(path, "rb") as f:
        return base64.b64encode(f.read()).decode("utf-8")


def write_embedded_html(entries, out_dir):
    """
    entries: list of (cell_id, mp4_path)
    """
    html_path = osp.join(out_dir, "index_embedded.html")

    with open(html_path, "w") as f:
        f.write(
            """
<!DOCTYPE html>
<html>
<head>
<title>Tracked Cell MIP Gallery</title>
<style>
body {
    font-family: Arial, sans-serif;
    background: #121212;  /* dark page background */
    color: #ffffff;       /* default text color for contrast */
    margin: 20px;
}
h1 {
    text-align: center;
    color: #ffffff;
}
.grid {
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(260px, 1fr));
    gap: 20px;
}
.card {
    background: #1e1e1e;  /* slightly lighter gray for cards */
    padding: 12px;
    border-radius: 12px;
    box-shadow: 0 4px 12px rgba(0,0,0,0.5);
    text-align: center;
    transition: transform 0.2s, box-shadow 0.2s;
}
.card:hover {
    transform: translateY(-5px);
    box-shadow: 0 8px 20px rgba(0,0,0,0.7);
}
video {
    width: 100%;
    border-radius: 6px;
    cursor: pointer;
}
h3 {
    font-size: 14px;
    margin: 8px 0;
    word-break: break-word;
    color: #dddddd;  /* light text on gray cards */
}
</style>
</head>
<body>

<h1>Tracked Cell MIP Videos</h1>

<div class="grid" id="gridContainer">
"""
        )

        for label, mp4_path in entries:
            b64 = file_to_base64(mp4_path)
            f.write(
                f"""
<div class="card">
  <h3>cell {label}</h3>
  <video loop muted playsinline preload="metadata">
    <source src="data:video/mp4;base64,{b64}" type="video/mp4">
    Your browser does not support the video tag.
  </video>
</div>
"""
            )

        f.write(
            """
</div>

<script>
// Hover to play video
document.querySelectorAll('video').forEach(v => {
    v.addEventListener('mouseenter', () => v.play());
    v.addEventListener('mouseleave', () => v.pause());
});
</script>

</body>
</html>
"""
        )

    return html_path





# --------------------------------------------------
# Main
# --------------------------------------------------
def main():
    args = parse_args()

    root_dir = osp.abspath(args.root)
    tracked_dir = osp.join(root_dir, "tracked_cells_smooth")

    print(f"\nRoot: {root_dir}")
    print(f"Scanning: {tracked_dir}")

    if not osp.isdir(tracked_dir):
        raise FileNotFoundError("tracked_cells_smooth not found")

    out_dir = osp.abspath(args.out)
    mp4_dir = osp.join(out_dir, "mp4")
    os.makedirs(mp4_dir, exist_ok=True)

    cell_dirs = sorted(glob(osp.join(tracked_dir, "cell_*")))
    if len(cell_dirs) == 0:
        raise RuntimeError("No cell_* directories found")
    
    # Sort cell_dirs numerically by the label
    cell_dirs = sorted(cell_dirs, key=lambda x: float(osp.basename(x).split('_')[1]))

    print(f"Found {len(cell_dirs)} cells\n")

    entries = []
    processed = 0
    skipped = 0

    for cell_dir in tqdm(cell_dirs, desc="Processing cells", unit="cell"):
        cell_id = osp.basename(cell_dir)
        label = int(float(cell_id.split('_')[1]))  # Extract as integer, e.g., 10.0 -> 10

        imgs = load_cell_frames(cell_dir)
        if len(imgs) == 0:
            tqdm.write(f"{cell_id}: no frames, skipping")
            skipped += 1
            continue

        mip_vid = make_mip_video(imgs)

        mp4_path = osp.join(mp4_dir, f"{cell_id}.mp4")
        save_mp4(mip_vid, mp4_path, args.fps)

        entries.append((label, mp4_path))
        processed += 1

    html_path = write_embedded_html(entries, out_dir)

    print("\nDone!")
    print(f"Cells processed: {processed}")
    print(f"Cells skipped: {skipped}")
    print(f"Open anywhere: {html_path}")


if __name__ == "__main__":
    main()
