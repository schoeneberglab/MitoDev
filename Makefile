VENV := .mitodev
# Use the system's python 3.10 or fallback to python3
PYTHON := $(shell which python3.10 || which python3)

.PHONY: venv clean

venv:
	@echo "Using Python from: $(PYTHON)"
	# Create the venv using the existing system python
	uv venv $(VENV) --python $(PYTHON)
	# Install packages
	. $(VENV)/bin/activate && \
	uv pip install \
		numpy \
		scikit-image \
		pyyaml \
		scikit-learn \
		matplotlib \
		pandas \
		cellpose==2.2.3 \
		open3d \
		imageio[ffmpeg] \
		gdown

clean:
	rm -rf $(VENV)
