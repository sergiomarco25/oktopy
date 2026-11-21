# oktopy with everything it runs: Python (analysis extras) and R (edgeR, limma, DESeq2, data.table), from
# envs/oktopy.yml. Laptop, workstation, cloud VM or (through Apptainer) an HPC cluster:
#   docker build -t oktopy .
#   docker run --rm -u "$(id -u):$(id -g)" -v "$PWD":/work oktopy example demo        # writes ./demo
#   docker run --rm -u "$(id -u):$(id -g)" -v "$PWD/demo":/work oktopy run            # the worked example, end to end
#   apptainer build oktopy.sif docker-daemon://oktopy:latest     # or docker://<registry>/oktopy:<tag>
#   apptainer exec --bind "$PWD" oktopy.sif oktopy run
FROM mambaorg/micromamba:1.5.10
COPY --chown=$MAMBA_USER:$MAMBA_USER envs/oktopy.yml /tmp/oktopy.yml
RUN micromamba install -y -n base -f /tmp/oktopy.yml && micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1
COPY --chown=$MAMBA_USER:$MAMBA_USER . /opt/oktopy
RUN python -m pip install --no-deps /opt/oktopy && oktopy --version \
    && Rscript -e 'suppressPackageStartupMessages({library(edgeR); library(limma); library(DESeq2); library(data.table)})'
WORKDIR /work
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "oktopy"]
CMD ["--help"]
