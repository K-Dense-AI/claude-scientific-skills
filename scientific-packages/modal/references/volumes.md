# Modal Volumes

## Overview

Modal Volumes provide high-performance distributed file systems for Modal applications. Designed for write-once, read-many workloads like ML model weights and distributed data processing.

## Creating Volumes

### Via CLI

```bash
modal volume create my-volume
```

For Volumes v2 (beta):
```bash
modal volume create --version=2 my-volume
```

### From Code

```python
vol = modal.Volume.from_name("my-volume", create_if_missing=True)

# For v2
vol = modal.Volume.from_name("my-volume", create_if_missing=True, version=2)
```

## Using Volumes

Attach to functions via mount points:

```python
vol = modal.Volume.from_name("my-volume")

@app.function(volumes={"/data": vol})
def run():
    with open("/data/xyz.txt", "w") as f:
        f.write("hello")
    vol.commit()  # Persist changes
```

## Commits and Reloads

### Commits

Persist changes to Volume:

```python
@app.function(volumes={"/data": vol})
def write_data():
    with open("/data/file.txt", "w") as f:
        f.write("data")
    vol.commit()  # Make changes visible to other containers
```

**Background commits**: Modal automatically commits Volume changes every few seconds and on container shutdown.

### Reloads

Fetch latest changes from other containers:

```python
@app.function(volumes={"/data": vol})
def read_data():
    vol.reload()  # Fetch latest changes
    with open("/data/file.txt", "r") as f:
        content = f.read()
```

At container creation, latest Volume state is mounted. Reload needed to see subsequent commits from other containers.

## Uploading Files

### Batch Upload (Efficient)

```python
vol = modal.Volume.from_name("my-volume")

with vol.batch_upload() as batch:
    batch.put_file("local-path.txt", "/remote-path.txt")
    batch.put_directory("/local/directory/", "/remote/directory")
    batch.put_file(io.BytesIO(b"some data"), "/foobar")
```

### Via Image

```python
image = modal.Image.debian_slim().add_local_dir(
    local_path="/home/user/my_dir",
    remote_path="/app"
)

@app.function(image=image)
def process():
    # Files available at /app
    ...
```

## Downloading Files

### Via CLI

```bash
modal volume get my-volume remote.txt local.txt
```

Max file size via CLI: No limit
Max file size via dashboard: 16 MB

### Via Python SDK

```python
vol = modal.Volume.from_name("my-volume")

for data in vol.read_file("path.txt"):
    print(data)
```

## Volume Performance

### Volumes v1

Best for:
- <50,000 files (recommended)
- <500,000 files (hard limit)
- Sequential access patterns
- <5 concurrent writers

### Volumes v2 (Beta)

Improved for:
- Unlimited files
- Hundreds of concurrent writers
- Random access patterns
- Large files (up to 1 TiB)

Current v2 limits:
- Max file size: 1 TiB
- Max files per directory: 32,768
- Unlimited directory depth

## Model Storage

### Saving Model Weights

```python
volume = modal.Volume.from_name("model-weights", create_if_missing=True)
MODEL_DIR = "/models"

@app.function(volumes={MODEL_DIR: volume})
def train():
    model = train_model()
    save_model(f"{MODEL_DIR}/my_model.pt", model)
    volume.commit()
```

### Loading Model Weights

```python
@app.function(volumes={MODEL_DIR: volume})
def inference(model_id: str):
    try:
        model = load_model(f"{MODEL_DIR}/{model_id}")
    except NotFound:
        volume.reload()  # Fetch latest models
        model = load_model(f"{MODEL_DIR}/{model_id}")
    return model.run(request)
```

## Model Checkpointing

Save checkpoints during long training jobs:

```python
volume = modal.Volume.from_name("checkpoints")
VOL_PATH = "/vol"

@app.function(
    gpu="A10G",
    timeout=2*60*60,  # 2 hours
    volumes={VOL_PATH: volume}
)
def finetune():
    from transformers import Seq2SeqTrainer, Seq2SeqTrainingArguments

    training_args = Seq2SeqTrainingArguments(
        output_dir=str(VOL_PATH / "model"),  # Checkpoints saved to Volume
        save_steps=100,
        # ... more args
    )

    trainer = Seq2SeqTrainer(model=model, args=training_args, ...)
    trainer.train()
```

Background commits ensure checkpoints persist even if training is interrupted.

## CLI Commands

```bash
# List files
modal volume ls my-volume

# Upload
modal volume put my-volume local.txt remote.txt

# Download
modal volume get my-volume remote.txt local.txt

# Copy within Volume
modal volume cp my-volume src.txt dst.txt

# Delete
modal volume rm my-volume file.txt

# List all volumes
modal volume list

# Delete volume
modal volume delete my-volume
```

## Ephemeral Volumes

Create temporary volumes that are garbage collected:

```python
with modal.Volume.ephemeral() as vol:
    sb = modal.Sandbox.create(
        volumes={"/cache": vol},
        app=my_app,
    )
    # Use volume
    # Automatically cleaned up when context exits
```

## Concurrent Access

### Concurrent Reads

Multiple containers can read simultaneously without issues.

### Concurrent Writes

Supported but:
- Avoid modifying same files concurrently
- Last write wins (data loss possible)
- v1: Limit to ~5 concurrent writers
- v2: Hundreds of concurrent writers supported

## Volume Errors

### "Volume Busy"

Cannot reload when files are open:

```python
# WRONG
f = open("/vol/data.txt", "r")
volume.reload()  # ERROR: volume busy
```

```python
# CORRECT
with open("/vol/data.txt", "r") as f:
    data = f.read()
# File closed before reload
volume.reload()
```

### "File Not Found"

Remember to use mount point:

```python
# WRONG - file saved to local disk
with open("/xyz.txt", "w") as f:
    f.write("data")

# CORRECT - file saved to Volume
with open("/data/xyz.txt", "w") as f:
    f.write("data")
```

## Upgrading from v1 to v2

No automated migration currently. Manual steps:

1. Create new v2 Volume
2. Copy data using `cp` or `rsync`
3. Update app to use new Volume

```bash
modal volume create --version=2 my-volume-v2
modal shell --volume my-volume --volume my-volume-v2

# In shell:
cp -rp /mnt/my-volume/. /mnt/my-volume-v2/.
sync /mnt/my-volume-v2
```

Warning: Deployed apps reference Volumes by ID. Re-deploy after creating new Volume.
