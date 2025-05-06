## method 1: Submitting the Job with AWS CLI

```
aws batch submit-job \
    --job-name nextflow-pipeline-job \
    --job-queue <your-job-queue-name> \
    --job-definition <your-job-definition-name> \
    --container-overrides '{
        "command": [
            "/bin/bash", "-c",
            "aws s3 cp s3://mybucket/scripts /app/scripts/ --recursive && nextflow run /app/scripts/main.nf -profile \"aws,Fargate\" -c /app/scripts/nextflow.config -bucket-dir s3://mybucket/demo/work"
        ]
    }'
```

## method 2: Submit the job with environment variables

- update Dockerfile
    <Details>
    
    ```
        FROM python:3.9-slim
        
        # Install system dependencies (procps provides ps command for NF to track process)
        RUN apt-get update && apt-get install -y libglib2.0-0 libsm6 libxext6 libxrender-dev procps curl unzip openjdk-17-jre \
            && rm -rf /var/lib/apt/lists/*
        
        # Install AWS CLI
        RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" \
            && unzip awscliv2.zip \
            && ./aws/install \
            && rm -rf awscliv2.zip aws
        
        # Install Nextflow
        RUN curl -s https://get.nextflow.io | bash \
            && mv nextflow /usr/local/bin/ \
            && chmod +x /usr/local/bin/nextflow
        
        # Install Python packages (if needed for kaleido or other packages)
        RUN pip install --no-cache-dir -U pandas scikit-learn plotly kaleido
        
        # Set working directory
        WORKDIR /app
        
        # Create a shell script to handle S3 copy and Nextflow execution
        RUN echo '#!/bin/bash\n\
        aws s3 cp "${S3_SCRIPTS_PATH:-s3://mybucket/scripts}" /app/scripts/ --recursive\n\
        nextflow run /app/scripts/main.nf -profile "${NF_PROFILE:-aws,Fargate}" -c /app/scripts/nextflow.config -bucket-dir "${NF_BUCKET_DIR:-s3://mybucket/demo/work}"\n\
        ' > /app/run.sh \
            && chmod +x /app/run.sh
        
        # Set the entrypoint to the shell script
        ENTRYPOINT ["/app/run.sh"]    
    ```
    
    </Details>

- submit job

```
aws batch submit-job \
    --job-name nextflow-pipeline-job \
    --job-queue <your-job-queue-name> \
    --job-definition <your-job-definition-name> \
    --container-overrides '{
        "command": ["/app/run.sh"],
        "environment": [
            {"name": "S3_SCRIPTS_PATH", "value": "s3://mybucket/scripts"},
            {"name": "NF_PROFILE", "value": "aws,Fargate"},
            {"name": "NF_BUCKET_DIR", "value": "s3://mybucket/demo/work"}
        ]
    }'
```

