# Secrets and Environment Variables

## Creating Secrets

### Via Dashboard

Create secrets at https://modal.com/secrets

Templates available for:
- Database credentials (Postgres, MongoDB)
- Cloud providers (AWS, GCP, Azure)
- ML platforms (Weights & Biases, Hugging Face)
- And more

### Via CLI

```bash
# Create secret with key-value pairs
modal secret create my-secret KEY1=value1 KEY2=value2

# Use environment variables
modal secret create db-secret PGHOST=uri PGPASSWORD="$PGPASSWORD"

# List secrets
modal secret list

# Delete secret
modal secret delete my-secret
```

### Programmatically

From dictionary:

```python
if modal.is_local():
    local_secret = modal.Secret.from_dict({"FOO": os.environ["LOCAL_FOO"]})
else:
    local_secret = modal.Secret.from_dict({})

@app.function(secrets=[local_secret])
def some_function():
    import os
    print(os.environ["FOO"])
```

From .env file:

```python
@app.function(secrets=[modal.Secret.from_dotenv()])
def some_function():
    import os
    print(os.environ["USERNAME"])
```

## Using Secrets

Inject secrets into functions:

```python
@app.function(secrets=[modal.Secret.from_name("my-secret")])
def some_function():
    import os
    secret_key = os.environ["MY_PASSWORD"]
    # Use secret
    ...
```

### Multiple Secrets

```python
@app.function(secrets=[
    modal.Secret.from_name("database-creds"),
    modal.Secret.from_name("api-keys"),
])
def other_function():
    # All keys from both secrets available
    ...
```

Later secrets override earlier ones if keys clash.

## Environment Variables

### Reserved Runtime Variables

**All Containers**:
- `MODAL_CLOUD_PROVIDER` - Cloud provider (AWS/GCP/OCI)
- `MODAL_IMAGE_ID` - Image ID
- `MODAL_REGION` - Region identifier (e.g., us-east-1)
- `MODAL_TASK_ID` - Container task ID

**Function Containers**:
- `MODAL_ENVIRONMENT` - Modal Environment name
- `MODAL_IS_REMOTE` - Set to '1' in remote containers
- `MODAL_IDENTITY_TOKEN` - OIDC token for function identity

**Sandbox Containers**:
- `MODAL_SANDBOX_ID` - Sandbox ID

### Setting Environment Variables

Via Image:

```python
image = modal.Image.debian_slim().env({"PORT": "6443"})

@app.function(image=image)
def my_function():
    import os
    port = os.environ["PORT"]
```

Via Secrets:

```python
secret = modal.Secret.from_dict({"API_KEY": "secret-value"})

@app.function(secrets=[secret])
def my_function():
    import os
    api_key = os.environ["API_KEY"]
```

## Common Secret Patterns

### AWS Credentials

```python
aws_secret = modal.Secret.from_name("my-aws-secret")

@app.function(secrets=[aws_secret])
def use_aws():
    import boto3
    s3 = boto3.client('s3')
    # AWS_ACCESS_KEY_ID and AWS_SECRET_ACCESS_KEY automatically used
```

### Hugging Face Token

```python
hf_secret = modal.Secret.from_name("huggingface")

@app.function(secrets=[hf_secret])
def download_model():
    from transformers import AutoModel
    # HF_TOKEN automatically used for authentication
    model = AutoModel.from_pretrained("private-model")
```

### Database Credentials

```python
db_secret = modal.Secret.from_name("postgres-creds")

@app.function(secrets=[db_secret])
def query_db():
    import psycopg2
    conn = psycopg2.connect(
        host=os.environ["PGHOST"],
        port=os.environ["PGPORT"],
        user=os.environ["PGUSER"],
        password=os.environ["PGPASSWORD"],
    )
```

## Best Practices

1. **Never hardcode secrets** - Always use Modal Secrets
2. **Use specific secrets** - Create separate secrets for different purposes
3. **Rotate secrets regularly** - Update secrets periodically
4. **Minimal scope** - Only attach secrets to functions that need them
5. **Environment-specific** - Use different secrets for dev/staging/prod

## Security Notes

- Secrets are encrypted at rest
- Only available to functions that explicitly request them
- Not logged or exposed in dashboards
- Can be scoped to specific environments
