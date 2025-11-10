# Common Patterns for Scientific Computing

## Machine Learning Model Inference

### Basic Model Serving

```python
import modal

app = modal.App("ml-inference")

image = (
    modal.Image.debian_slim()
    .uv_pip_install("torch", "transformers")
)

@app.cls(
    image=image,
    gpu="L40S",
)
class Model:
    @modal.enter()
    def load_model(self):
        from transformers import AutoModel, AutoTokenizer
        self.model = AutoModel.from_pretrained("bert-base-uncased")
        self.tokenizer = AutoTokenizer.from_pretrained("bert-base-uncased")

    @modal.method()
    def predict(self, text: str):
        inputs = self.tokenizer(text, return_tensors="pt")
        outputs = self.model(**inputs)
        return outputs.last_hidden_state.mean(dim=1).tolist()

@app.local_entrypoint()
def main():
    model = Model()
    result = model.predict.remote("Hello world")
    print(result)
```

### Model Serving with Volume

```python
volume = modal.Volume.from_name("models", create_if_missing=True)
MODEL_PATH = "/models"

@app.cls(
    image=image,
    gpu="A100",
    volumes={MODEL_PATH: volume}
)
class ModelServer:
    @modal.enter()
    def load(self):
        import torch
        self.model = torch.load(f"{MODEL_PATH}/model.pt")
        self.model.eval()

    @modal.method()
    def infer(self, data):
        import torch
        with torch.no_grad():
            return self.model(torch.tensor(data)).tolist()
```

## Batch Processing

### Parallel Data Processing

```python
@app.function(
    image=modal.Image.debian_slim().uv_pip_install("pandas", "numpy"),
    cpu=2.0,
    memory=8192
)
def process_batch(batch_id: int):
    import pandas as pd

    # Load batch
    df = pd.read_csv(f"s3://bucket/batch_{batch_id}.csv")

    # Process
    result = df.apply(lambda row: complex_calculation(row), axis=1)

    # Save result
    result.to_csv(f"s3://bucket/results_{batch_id}.csv")

    return batch_id

@app.local_entrypoint()
def main():
    # Process 100 batches in parallel
    results = list(process_batch.map(range(100)))
    print(f"Processed {len(results)} batches")
```

### Batch Processing with Progress

```python
@app.function()
def process_item(item_id: int):
    # Expensive processing
    result = compute_something(item_id)
    return result

@app.local_entrypoint()
def main():
    items = list(range(1000))

    print(f"Processing {len(items)} items...")
    results = []
    for i, result in enumerate(process_item.map(items)):
        results.append(result)
        if (i + 1) % 100 == 0:
            print(f"Completed {i + 1}/{len(items)}")

    print("All items processed!")
```

## Data Analysis Pipeline

### ETL Pipeline

```python
volume = modal.Volume.from_name("data-pipeline")
DATA_PATH = "/data"

@app.function(
    image=modal.Image.debian_slim().uv_pip_install("pandas", "polars"),
    volumes={DATA_PATH: volume},
    cpu=4.0,
    memory=16384
)
def extract_transform_load():
    import polars as pl

    # Extract
    raw_data = pl.read_csv(f"{DATA_PATH}/raw/*.csv")

    # Transform
    transformed = (
        raw_data
        .filter(pl.col("value") > 0)
        .group_by("category")
        .agg([
            pl.col("value").mean().alias("avg_value"),
            pl.col("value").sum().alias("total_value")
        ])
    )

    # Load
    transformed.write_parquet(f"{DATA_PATH}/processed/data.parquet")
    volume.commit()

    return transformed.shape

@app.function(schedule=modal.Cron("0 2 * * *"))
def daily_pipeline():
    result = extract_transform_load.remote()
    print(f"Processed data shape: {result}")
```

## GPU-Accelerated Computing

### Distributed Training

```python
@app.function(
    gpu="A100:2",
    image=modal.Image.debian_slim().uv_pip_install("torch", "accelerate"),
    timeout=7200,
)
def train_model():
    import torch
    from torch.nn.parallel import DataParallel

    # Load data
    train_loader = get_data_loader()

    # Initialize model
    model = MyModel()
    model = DataParallel(model)
    model = model.cuda()

    # Train
    optimizer = torch.optim.Adam(model.parameters())
    for epoch in range(10):
        for batch in train_loader:
            loss = train_step(model, batch, optimizer)
            print(f"Epoch {epoch}, Loss: {loss}")

    return "Training complete"
```

### GPU Batch Inference

```python
@app.function(
    gpu="L40S",
    image=modal.Image.debian_slim().uv_pip_install("torch", "transformers")
)
def batch_inference(texts: list[str]):
    from transformers import pipeline

    classifier = pipeline("sentiment-analysis", device=0)
    results = classifier(texts, batch_size=32)

    return results

@app.local_entrypoint()
def main():
    # Process 10,000 texts
    texts = load_texts()

    # Split into chunks of 100
    chunks = [texts[i:i+100] for i in range(0, len(texts), 100)]

    # Process in parallel on multiple GPUs
    all_results = []
    for results in batch_inference.map(chunks):
        all_results.extend(results)

    print(f"Processed {len(all_results)} texts")
```

## Scientific Computing

### Molecular Dynamics Simulation

```python
@app.function(
    image=modal.Image.debian_slim().apt_install("openmpi-bin").uv_pip_install("mpi4py", "numpy"),
    cpu=16.0,
    memory=65536,
    timeout=7200,
)
def run_simulation(config: dict):
    import numpy as np

    # Initialize system
    positions = initialize_positions(config["n_particles"])
    velocities = initialize_velocities(config["temperature"])

    # Run MD steps
    for step in range(config["n_steps"]):
        forces = compute_forces(positions)
        velocities += forces * config["dt"]
        positions += velocities * config["dt"]

        if step % 1000 == 0:
            energy = compute_energy(positions, velocities)
            print(f"Step {step}, Energy: {energy}")

    return positions, velocities
```

### Distributed Monte Carlo

```python
@app.function(cpu=2.0)
def monte_carlo_trial(trial_id: int, n_samples: int):
    import random

    count = sum(1 for _ in range(n_samples)
                if random.random()**2 + random.random()**2 <= 1)

    return count

@app.local_entrypoint()
def estimate_pi():
    n_trials = 100
    n_samples_per_trial = 1_000_000

    # Run trials in parallel
    results = list(monte_carlo_trial.map(
        range(n_trials),
        [n_samples_per_trial] * n_trials
    ))

    total_count = sum(results)
    total_samples = n_trials * n_samples_per_trial

    pi_estimate = 4 * total_count / total_samples
    print(f"Estimated π = {pi_estimate}")
```

## Data Processing with Volumes

### Image Processing Pipeline

```python
volume = modal.Volume.from_name("images")
IMAGE_PATH = "/images"

@app.function(
    image=modal.Image.debian_slim().uv_pip_install("Pillow", "numpy"),
    volumes={IMAGE_PATH: volume}
)
def process_image(filename: str):
    from PIL import Image
    import numpy as np

    # Load image
    img = Image.open(f"{IMAGE_PATH}/raw/{filename}")

    # Process
    img_array = np.array(img)
    processed = apply_filters(img_array)

    # Save
    result_img = Image.fromarray(processed)
    result_img.save(f"{IMAGE_PATH}/processed/{filename}")

    return filename

@app.function(volumes={IMAGE_PATH: volume})
def process_all_images():
    import os

    # Get all images
    filenames = os.listdir(f"{IMAGE_PATH}/raw")

    # Process in parallel
    results = list(process_image.map(filenames))

    volume.commit()
    return f"Processed {len(results)} images"
```

## Web API for Scientific Computing

```python
image = modal.Image.debian_slim().uv_pip_install("fastapi[standard]", "numpy", "scipy")

@app.function(image=image)
@modal.fastapi_endpoint(method="POST")
def compute_statistics(data: dict):
    import numpy as np
    from scipy import stats

    values = np.array(data["values"])

    return {
        "mean": float(np.mean(values)),
        "median": float(np.median(values)),
        "std": float(np.std(values)),
        "skewness": float(stats.skew(values)),
        "kurtosis": float(stats.kurtosis(values))
    }
```

## Scheduled Data Collection

```python
@app.function(
    schedule=modal.Cron("*/30 * * * *"),  # Every 30 minutes
    secrets=[modal.Secret.from_name("api-keys")],
    volumes={"/data": modal.Volume.from_name("sensor-data")}
)
def collect_sensor_data():
    import requests
    import json
    from datetime import datetime

    # Fetch from API
    response = requests.get(
        "https://api.example.com/sensors",
        headers={"Authorization": f"Bearer {os.environ['API_KEY']}"}
    )

    data = response.json()

    # Save with timestamp
    timestamp = datetime.now().isoformat()
    with open(f"/data/{timestamp}.json", "w") as f:
        json.dump(data, f)

    volume.commit()

    return f"Collected {len(data)} sensor readings"
```

## Best Practices

### Use Classes for Stateful Workloads

```python
@app.cls(gpu="A100")
class ModelService:
    @modal.enter()
    def setup(self):
        # Load once, reuse across requests
        self.model = load_heavy_model()

    @modal.method()
    def predict(self, x):
        return self.model(x)
```

### Batch Similar Workloads

```python
@app.function()
def process_many(items: list):
    # More efficient than processing one at a time
    return [process(item) for item in items]
```

### Use Volumes for Large Datasets

```python
# Store large datasets in volumes, not in image
volume = modal.Volume.from_name("dataset")

@app.function(volumes={"/data": volume})
def train():
    data = load_from_volume("/data/training.parquet")
    model = train_model(data)
```

### Profile Before Scaling to GPUs

```python
# Test on CPU first
@app.function(cpu=4.0)
def test_pipeline():
    ...

# Then scale to GPU if needed
@app.function(gpu="A100")
def gpu_pipeline():
    ...
```
