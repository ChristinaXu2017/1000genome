# Nextflow Pipeline on AWS Fargate
A Nextflow pipeline is proposed to Examine genomic variation across populations with AWS. Here, both Nextflow head and task jobs are running on AWS Fargate, leveraging Spot instances, and automating infrastructure with Terraform, our approach overcomes the scalability, cost, and speed limitations of prior methods.

## AWS Cloud Infrasture
The folder named terraform provides the template to set up cloud infrastructure. 
- step 1: set up AWS cloud credential in your local PC. 
  ```
  # eg. less ~/.aws/config 
  [default]
    region = us-east-1 
    output = json

  # eg. less ~/.aws/credentials 
  [default]
    aws_access_key_id=ASIAQT...X
    aws_secret_access_key=1mSq9i0+Vcr...A/u
    aws_session_token=IQoJ/...=
  ```
- step 2: terraform deployment
  ```
  cd ./aws.run/terraform
  terraform init
  terraform apply
  ```
- step 3: now you can get aws resource through AWS console or terraform command. eg.
  - job_definition arn to launch docker container running PCA and bcftool, refer to "scripts/nextflow.config"
  - name of your private bucket, job_queue and Nextflow head job_definition, refer to "submit.sh"
    
## Nextflow Pipeline


 

