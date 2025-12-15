# MitoDev – 4D Cell Segmentation and Tracking

This repository contains the code for the paper **production-ready pipeline**:

* Downloading microscopy data
* Segmenting cells using **Cellpose**
* Running mitochondrial tracking
* Generating **per-cell MIP videos**
* Visualizing all cells in a **HTML dashboard**

Some stages are **computationally expensive** and may take **hours to complete**, depending on dataset size and hardware.

---

## 📁 Repository Structure

```
.
├── config.yaml
├── main.py
├── visualise_in_html.py
├── scripts/
│   ├── download.sh
│   ├── run_cellpose.sh
├── data/
├── README.md
```

---

## ⚙️ Environment Setup

### 1️⃣ Install `uv`

`uv` is used for fast Python environment and dependency management.

```bash
pip install uv
```

---

### 2️⃣ Create Virtual Environment

```bash
make venv
```

This will:

* Install Python 3.10
* Create a virtual environment at `.mitodev/`
* Install all required dependencies

---

### 3️⃣ Activate Environment

```bash
source .mitodev/bin/activate
```

---

## 📥 Data Download

Download the raw dataset into the `data/` directory:

```bash
bash scripts/download.sh ./data
```

⏱ **Runtime**:

* Depends on network speed
* Can take several minutes for large datasets

---

## 🧠 Cell Segmentation (Cellpose)

Run Cellpose on the processed data directory:

```bash
bash scripts/run_cellpose.sh "./data/20231221 Gillian Lung Organoid/Sample 1/1/Processed_Data" 0 10
```

Arguments:

* **Processed data path**
* **Start index**
* **End index (exclusive)**

⏱ **Runtime**:

* **Can take hours** for large datasets
* Strongly recommended to run on a machine with GPU support

---

## 🧬 Mitochondrial Tracking

Run the main processing pipeline:

```bash
python -m main ./config.yaml "20231221 Gillian Lung Organoid/Sample 2" 1 --minimal
```

Arguments:

* `config.yaml`: Pipeline configuration
* Sample name
* Sample index
* `--minimal`: Runs a reduced output version (recommended for testing)

⏱ **Runtime**:

* **Several hours** depending on:

  * Number of cells
  * Timepoints

---

## 🎥 Cell Visualization (MIP Videos → HTML)

Generate **per-cell MIP videos** and embed them into a **self-contained HTML page**:

```bash
python visualise_in_html.py --root "./data/20231221 Gillian Lung Organoid/Sample 1/1/single_cells/"
```

This will:

* Traverse all tracked cells
* Generate MIP MP4 videos
* Embed all videos directly into an HTML file

📄 Output:

```
index_embedded.html
```

You can open this file **on any system**, even without access to the data or Python.

---

## 🖥️ Hardware Recommendations

| Task          | Recommendation           |
| ------------- | ------------------------ |
| Cellpose      | GPU strongly recommended |
| Main pipeline | GPU + high RAM           |
| Visualization | CPU sufficient           |

---

## ⚠️ Notes

* Some pipeline stages are **long-running** by design

---

## 📬 Contact

For issues, improvements, or extensions, please open a GitHub issue or contact the maintainers.

---