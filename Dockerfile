# Use an official Miniconda base image
FROM continuumio/miniconda3

# Set the working directory in the container
WORKDIR /app

# Copy the current directory contents into the container at /streamlit_app
COPY . /app

# Optionally update Conda
RUN conda update -n base -c defaults conda

# Create the Conda environment
RUN conda env create -f env.yml

# Make RUN commands use the new environment
SHELL ["conda", "run", "-n", "env310", "/bin/bash", "-c"]

# Ensure the environment is activated
RUN echo "conda activate env310" >> ~/.bashrc

# Install vim
RUN apt-get update && apt-get install -y vim

# Install the package in editable mode
RUN pip install -e .

# Start an interactive shell
CMD ["/bin/bash"]
