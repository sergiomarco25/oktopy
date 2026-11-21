# Running oktopy: laptop, cluster, container, cloud

The steps are the same everywhere (`oktopy run` lists them: `oktopy run --list`); what changes is who schedules them.

## Laptop or workstation

```bash
cd my-field && oktopy run
```

One step after another for the datasets of `[pipeline] datasets`. A failing step stops the run and says how to
resume (`oktopy run --from STEP`); `--datasets A,B` limits the per-dataset steps. The worked example takes about five
minutes on four cores and 8 GB of memory. A real dataset of a few hundred thousand nuclei needs tens of GB for
harmonisation; set `data_root` to a disk with room.

With Snakemake (installed by `envs/oktopy.yml`), only what changed is rebuilt and independent datasets run in
parallel:

```bash
snakemake -n
```

```bash
snakemake --cores 4
```

`workflow/profiles/default/config.yaml` is the laptop profile (4 cores, re-runs triggered by file times; after a
code or config change, `--forcerun <rule>`). `snakemake robustness` adds the shuffled-label calibration and power.

## SLURM cluster

1. Install the environment on the cluster (`mamba env create -f envs/oktopy.yml`, then `pip install` oktopy at the
   version the project pins), or use the container (below).
2. Put the project folder in your home or project space and `data_root` (or `<PREFIX>_DATA_ROOT`) on a file system
   every compute node reads. Job files and temporary data should not go to a node-local `/tmp`.
3. In `workflow/profiles/slurm/config.yaml` set `slurm_partition` and `qos` (`sinfo` and `sacctmgr show qos` list
   them) and in `config/project.toml` the memory per dataset (`[pipeline] mem_gb`).
4. Run the workflow controller where it may run for hours (a login node with `tmux`, or a small job):

```bash
snakemake --profile workflow/profiles/slurm -n
```

```bash
snakemake --profile workflow/profiles/slurm
```

One job per dataset and step; `harmonise` and `tests` get the dataset's memory, `de_cache` 16 GB. If preempted jobs
are cancelled rather than requeued on your cluster, use a QOS without preemption: the controller otherwise waits for
a job that will never finish.

## Containers

The image has Python, R (edgeR, limma, DESeq2, data.table) and oktopy; `oktopy` is its entry point.

```bash
docker build -t oktopy .
```

```bash
docker run --rm -u "$(id -u):$(id -g)" -v "$PWD/my-field":/work oktopy run
```

On a cluster without Docker, build an Apptainer (Singularity) image from `apptainer.def` (the same steps as the
Dockerfile; add `--fakeroot` if your site requires it), or convert a Docker image you built elsewhere
(`apptainer build oktopy.sif docker-daemon://oktopy:latest`), and run commands through it:

```bash
apptainer build oktopy.sif apptainer.def
```

```bash
apptainer exec --bind "$PWD" oktopy.sif oktopy run
```

The project's own package (dataset loaders) is read from the mounted project folder, so the image does not need to
change when a project adds a dataset.

## Cloud

Any VM with Docker runs the image as above; size memory to the largest dataset. A built store needs nothing but
Python's standard library to serve: `oktopy http` (MCP Streamable HTTP plus REST and OpenAPI) with
`<PREFIX>_HTTP_HOST=0.0.0.0`, `PORT` and `<PREFIX>_STORE_DB` set, behind any container service. It is read-only
unless `<PREFIX>_MCP_ALLOW_WRITE=1`; `<PREFIX>_MCP_TOKEN` adds a shared secret.

## Serving a store to AI tools

In Claude Code, add the stdio server to a project's `.mcp.json`:

```json
{"mcpServers": {"my-field": {"command": "oktopy", "args": ["--project", "/path/to/my-field", "serve"]}}}
```

The server answers from the store and cites claim ids; `config/server_instructions.md` tells clients how to use it.
