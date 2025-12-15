VENV := .mitodev

.PHONY: venv clean

venv:
	uv python install 3.10
	uv venv $(VENV) --python 3.10
	. $(VENV)/bin/activate && \
	uv pip install \
		numpy \
		scikit-image \
		pyyaml \
		scikit-learn \
		matplotlib \
		pandas \
		cellpose \
		open3d \
		imageio[ffmpeg] \
		gdown

clean:
	rm -rf $(VENV)
