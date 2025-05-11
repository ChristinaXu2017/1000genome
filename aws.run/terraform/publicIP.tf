# Data source to fetch the default VPC
data "aws_vpc" "default" {
  default = true
}

# Data source to fetch public subnets in the default VPC
data "aws_subnets" "default_public" {
  filter {
    name   = "vpc-id"
    values = [data.aws_vpc.default.id]
  }
  filter {
    name   = "default-for-az"
    values = ["true"]
  }
}

# Security group for Fargate tasks
resource "aws_security_group" "fargate_sg" {
  name        = "batch-fargate-sg"
  description = "Security group for AWS Batch Fargate tasks"
  vpc_id      = data.aws_vpc.default.id

  # Allow outbound traffic to access ECR and other AWS services
  egress {
    from_port   = 0
    to_port     = 0
    protocol    = "-1"
    cidr_blocks = ["0.0.0.0/0"]
  }
}


