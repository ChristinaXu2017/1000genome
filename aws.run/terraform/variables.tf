variable "region" {
  default = "us-east-1"
}

variable "environment" {
  default = "production"
}

variable "project" {
  description = "The name of the project"
  type        = string
  default     = "1000genome" 
}

variable "bucket_output" {
  description = "bucket for results and scripts"
  type        = string
  default     = "1000genome" 

}



