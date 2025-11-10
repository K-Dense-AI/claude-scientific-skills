# Scheduled Jobs and Cron

## Basic Scheduling

Schedule functions to run automatically at regular intervals or specific times.

### Simple Daily Schedule

```python
import modal

app = modal.App()

@app.function(schedule=modal.Period(days=1))
def daily_task():
    print("Running daily task")
    # Process data, send reports, etc.
```

Deploy to activate:
```bash
modal deploy script.py
```

Function runs every 24 hours from deployment time.

## Schedule Types

### Period Schedules

Run at fixed intervals from deployment time:

```python
# Every 5 hours
@app.function(schedule=modal.Period(hours=5))
def every_5_hours():
    ...

# Every 30 minutes
@app.function(schedule=modal.Period(minutes=30))
def every_30_minutes():
    ...

# Every day
@app.function(schedule=modal.Period(days=1))
def daily():
    ...
```

**Note**: Redeploying resets the period timer.

### Cron Schedules

Run at specific times using cron syntax:

```python
# Every Monday at 8 AM UTC
@app.function(schedule=modal.Cron("0 8 * * 1"))
def weekly_report():
    ...

# Daily at 6 AM New York time
@app.function(schedule=modal.Cron("0 6 * * *", timezone="America/New_York"))
def morning_report():
    ...

# Every hour on the hour
@app.function(schedule=modal.Cron("0 * * * *"))
def hourly():
    ...

# Every 15 minutes
@app.function(schedule=modal.Cron("*/15 * * * *"))
def quarter_hourly():
    ...
```

**Cron syntax**: `minute hour day month day_of_week`
- Minute: 0-59
- Hour: 0-23
- Day: 1-31
- Month: 1-12
- Day of week: 0-6 (0 = Sunday)

### Timezone Support

Specify timezone for cron schedules:

```python
@app.function(schedule=modal.Cron("0 9 * * *", timezone="Europe/London"))
def uk_morning_task():
    ...

@app.function(schedule=modal.Cron("0 17 * * 5", timezone="Asia/Tokyo"))
def friday_evening_jp():
    ...
```

## Deployment

### Deploy Scheduled Functions

```bash
modal deploy script.py
```

Scheduled functions persist until explicitly stopped.

### Programmatic Deployment

```python
if __name__ == "__main__":
    app.deploy()
```

## Monitoring

### View Execution Logs

Check https://modal.com/apps for:
- Past execution logs
- Execution history
- Failure notifications

### Run Manually

Trigger scheduled function immediately via dashboard "Run now" button.

## Schedule Management

### Pausing Schedules

Schedules cannot be paused. To stop:
1. Remove `schedule` parameter
2. Redeploy app

### Updating Schedules

Change schedule parameters and redeploy:

```python
# Update from daily to weekly
@app.function(schedule=modal.Period(days=7))
def task():
    ...
```

```bash
modal deploy script.py
```

## Common Patterns

### Data Pipeline

```python
@app.function(
    schedule=modal.Cron("0 2 * * *"),  # 2 AM daily
    timeout=3600,                       # 1 hour timeout
)
def etl_pipeline():
    # Extract data from sources
    data = extract_data()

    # Transform data
    transformed = transform_data(data)

    # Load to warehouse
    load_to_warehouse(transformed)
```

### Model Retraining

```python
volume = modal.Volume.from_name("models")

@app.function(
    schedule=modal.Cron("0 0 * * 0"),  # Weekly on Sunday midnight
    gpu="A100",
    timeout=7200,                       # 2 hours
    volumes={"/models": volume}
)
def retrain_model():
    # Load latest data
    data = load_training_data()

    # Train model
    model = train(data)

    # Save new model
    save_model(model, "/models/latest.pt")
    volume.commit()
```

### Report Generation

```python
@app.function(
    schedule=modal.Cron("0 9 * * 1"),  # Monday 9 AM
    secrets=[modal.Secret.from_name("email-creds")]
)
def weekly_report():
    # Generate report
    report = generate_analytics_report()

    # Send email
    send_email(
        to="team@company.com",
        subject="Weekly Analytics Report",
        body=report
    )
```

### Data Cleanup

```python
@app.function(schedule=modal.Period(hours=6))
def cleanup_old_data():
    # Remove data older than 30 days
    cutoff = datetime.now() - timedelta(days=30)
    delete_old_records(cutoff)
```

## Configuration with Secrets and Volumes

Scheduled functions support all function parameters:

```python
vol = modal.Volume.from_name("data")
secret = modal.Secret.from_name("api-keys")

@app.function(
    schedule=modal.Cron("0 */6 * * *"),  # Every 6 hours
    secrets=[secret],
    volumes={"/data": vol},
    cpu=4.0,
    memory=16384,
)
def sync_data():
    import os

    api_key = os.environ["API_KEY"]

    # Fetch from external API
    data = fetch_external_data(api_key)

    # Save to volume
    with open("/data/latest.json", "w") as f:
        json.dump(data, f)

    vol.commit()
```

## Dynamic Scheduling

Update schedules programmatically:

```python
@app.function()
def main_task():
    ...

@app.function(schedule=modal.Cron("0 6 * * *", timezone="America/New_York"))
def enable_high_traffic_mode():
    main_task.update_autoscaler(min_containers=5)

@app.function(schedule=modal.Cron("0 22 * * *", timezone="America/New_York"))
def disable_high_traffic_mode():
    main_task.update_autoscaler(min_containers=0)
```

## Error Handling

Scheduled functions that fail will:
- Show failure in dashboard
- Send notifications (configurable)
- Retry on next scheduled run

```python
@app.function(
    schedule=modal.Cron("0 * * * *"),
    retries=3,  # Retry failed runs
    timeout=1800
)
def robust_task():
    try:
        perform_task()
    except Exception as e:
        # Log error
        print(f"Task failed: {e}")
        # Optionally send alert
        send_alert(f"Scheduled task failed: {e}")
        raise
```

## Best Practices

1. **Set timeouts**: Always specify timeout for scheduled functions
2. **Use appropriate schedules**: Period for relative timing, Cron for absolute
3. **Monitor failures**: Check dashboard regularly for failed runs
4. **Idempotent operations**: Design tasks to handle reruns safely
5. **Resource limits**: Set appropriate CPU/memory for scheduled workloads
6. **Timezone awareness**: Specify timezone for cron schedules
