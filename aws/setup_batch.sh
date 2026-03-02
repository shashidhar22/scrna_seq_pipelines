#!/usr/bin/env bash
# =============================================================================
# AWS Batch Infrastructure Setup — scRNA-seq Pipeline
# =============================================================================
#
# Creates: ECR repos, IAM roles, launch templates, compute environments,
# job queues for running the scRNA-seq pipeline on AWS Batch with Spot.
#
# Prerequisites:
#   - AWS CLI v2 configured with appropriate permissions
#   - jq installed
#
# Usage:
#   bash aws/setup_batch.sh
#
# To tear down:
#   bash aws/setup_batch.sh --teardown
#
# SSO/admin-managed accounts (no IAM permissions):
#   Have your admin create the IAM resources per aws/iam_admin_instructions.md,
#   then pass the ARNs as environment variables:
#
#   export BATCH_SERVICE_ROLE_ARN="arn:aws:iam::123456789012:role/scrna-batch-service-role"
#   export ECS_INSTANCE_PROFILE_ARN="arn:aws:iam::123456789012:instance-profile/scrna-ecs-instance-profile"
#   export SPOT_FLEET_ROLE_ARN="arn:aws:iam::123456789012:role/scrna-spot-fleet-role"
#   export S3_BUCKET="my-bucket-name"
#   bash aws/setup_batch.sh
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Configuration — update these for your account
# ---------------------------------------------------------------------------
AWS_REGION="${AWS_REGION:-us-east-1}"
ACCOUNT_ID=$(aws sts get-caller-identity --query Account --output text)
PROJECT="scrna"
S3_BUCKET="${S3_BUCKET:-kstme-scrna}"

# Networking — update with your VPC/subnet IDs
VPC_SUBNETS="${VPC_SUBNETS:-}"            # comma-separated subnet IDs
SECURITY_GROUPS="${SECURITY_GROUPS:-}"     # comma-separated security group IDs

# Pre-existing IAM ARNs (for SSO/admin-managed accounts without IAM permissions)
BATCH_SERVICE_ROLE_ARN="${BATCH_SERVICE_ROLE_ARN:-}"
ECS_INSTANCE_PROFILE_ARN="${ECS_INSTANCE_PROFILE_ARN:-}"
SPOT_FLEET_ROLE_ARN="${SPOT_FLEET_ROLE_ARN:-}"

# Determine whether to create IAM resources or use pre-existing ARNs
_iam_vars_set=0
_iam_vars_names=("BATCH_SERVICE_ROLE_ARN" "ECS_INSTANCE_PROFILE_ARN" "SPOT_FLEET_ROLE_ARN")
for _var in "${_iam_vars_names[@]}"; do
    [[ -n "${!_var}" ]] && (( _iam_vars_set++ ))
done

if [[ "${_iam_vars_set}" -eq 0 ]]; then
    CREATE_IAM=true
elif [[ "${_iam_vars_set}" -eq 3 ]]; then
    CREATE_IAM=false
    # Validate ARN formats
    for _var in "${_iam_vars_names[@]}"; do
        _val="${!_var}"
        if [[ ! "${_val}" =~ ^arn:aws:iam::[0-9]{12}:(role|instance-profile)/ ]]; then
            echo "ERROR: ${_var} does not look like a valid IAM ARN."
            echo "  Got:      ${_val}"
            echo "  Expected: arn:aws:iam::<12-digit-account>:role/<name>"
            echo "            or arn:aws:iam::<12-digit-account>:instance-profile/<name>"
            exit 1
        fi
    done
else
    echo "ERROR: Either set all three IAM ARN variables or none of them."
    for _var in "${_iam_vars_names[@]}"; do
        _val="${!_var}"
        if [[ -n "${_val}" ]]; then
            echo "  ${_var} = ${_val}  (SET)"
        else
            echo "  ${_var}  (MISSING)"
        fi
    done
    echo ""
    echo "See aws/iam_admin_instructions.md for details."
    exit 1
fi

if [[ "${CREATE_IAM}" == "true" ]]; then
    _iam_mode="Creating roles (local account)"
else
    _iam_mode="Using pre-existing ARNs (SSO/admin-managed)"
fi

echo "============================================"
echo "AWS Batch Setup for scRNA-seq Pipeline"
echo "============================================"
echo "Account:  ${ACCOUNT_ID}"
echo "Region:   ${AWS_REGION}"
echo "Project:  ${PROJECT}"
echo "S3:       s3://${S3_BUCKET}/"
echo "IAM:      ${_iam_mode}"
echo "============================================"

if [[ -z "${VPC_SUBNETS}" || -z "${SECURITY_GROUPS}" ]]; then
    echo ""
    echo "ERROR: Set VPC_SUBNETS and SECURITY_GROUPS environment variables."
    echo "  export VPC_SUBNETS='subnet-abc123,subnet-def456'"
    echo "  export SECURITY_GROUPS='sg-abc123'"
    exit 1
fi

# ---------------------------------------------------------------------------
# Teardown mode
# ---------------------------------------------------------------------------
if [[ "${1:-}" == "--teardown" ]]; then
    echo ""
    echo "Tearing down AWS Batch resources..."
    echo ""

    echo "Disabling job queues..."
    aws batch update-job-queue --job-queue "${PROJECT}-cpu-queue" --state DISABLED --region "${AWS_REGION}" 2>/dev/null || true
    aws batch update-job-queue --job-queue "${PROJECT}-gpu-queue" --state DISABLED --region "${AWS_REGION}" 2>/dev/null || true
    sleep 10

    echo "Deleting job queues..."
    aws batch delete-job-queue --job-queue "${PROJECT}-cpu-queue" --region "${AWS_REGION}" 2>/dev/null || true
    aws batch delete-job-queue --job-queue "${PROJECT}-gpu-queue" --region "${AWS_REGION}" 2>/dev/null || true
    sleep 15

    echo "Disabling compute environments..."
    aws batch update-compute-environment --compute-environment "${PROJECT}-cpu-ce" --state DISABLED --region "${AWS_REGION}" 2>/dev/null || true
    aws batch update-compute-environment --compute-environment "${PROJECT}-gpu-ce" --state DISABLED --region "${AWS_REGION}" 2>/dev/null || true
    sleep 15

    echo "Deleting compute environments..."
    aws batch delete-compute-environment --compute-environment "${PROJECT}-cpu-ce" --region "${AWS_REGION}" 2>/dev/null || true
    aws batch delete-compute-environment --compute-environment "${PROJECT}-gpu-ce" --region "${AWS_REGION}" 2>/dev/null || true

    echo "Teardown complete. IAM roles and launch templates preserved."
    exit 0
fi

# ---------------------------------------------------------------------------
# 1. ECR Repositories
# ---------------------------------------------------------------------------
echo ""
echo "--- Creating ECR repositories ---"

for repo in "${PROJECT}-cpu" "${PROJECT}-gpu"; do
    if aws ecr describe-repositories --repository-names "${repo}" --region "${AWS_REGION}" &>/dev/null; then
        echo "ECR repo '${repo}' already exists"
    else
        aws ecr create-repository \
            --repository-name "${repo}" \
            --region "${AWS_REGION}" \
            --image-scanning-configuration scanOnPush=true \
            --encryption-configuration encryptionType=AES256
        echo "Created ECR repo: ${repo}"
    fi
done

# ---------------------------------------------------------------------------
# 2. IAM Roles
# ---------------------------------------------------------------------------
echo ""

if [[ "${CREATE_IAM}" == "true" ]]; then
    echo "--- Creating IAM roles ---"

    # 2a. Batch service role (allows Batch to manage EC2 instances)
    BATCH_SERVICE_ROLE="${PROJECT}-batch-service-role"
    if aws iam get-role --role-name "${BATCH_SERVICE_ROLE}" &>/dev/null; then
        echo "Role '${BATCH_SERVICE_ROLE}' already exists"
    else
        aws iam create-role \
            --role-name "${BATCH_SERVICE_ROLE}" \
            --assume-role-policy-document '{
                "Version": "2012-10-17",
                "Statement": [{
                    "Effect": "Allow",
                    "Principal": {"Service": "batch.amazonaws.com"},
                    "Action": "sts:AssumeRole"
                }]
            }'
        aws iam attach-role-policy \
            --role-name "${BATCH_SERVICE_ROLE}" \
            --policy-arn "arn:aws:iam::aws:policy/service-role/AWSBatchServiceRole"
        echo "Created role: ${BATCH_SERVICE_ROLE}"
    fi

    # 2b. ECS instance role (EC2 instances that run containers)
    ECS_INSTANCE_ROLE="${PROJECT}-ecs-instance-role"
    ECS_INSTANCE_PROFILE="${PROJECT}-ecs-instance-profile"
    if aws iam get-role --role-name "${ECS_INSTANCE_ROLE}" &>/dev/null; then
        echo "Role '${ECS_INSTANCE_ROLE}' already exists"
    else
        aws iam create-role \
            --role-name "${ECS_INSTANCE_ROLE}" \
            --assume-role-policy-document '{
                "Version": "2012-10-17",
                "Statement": [{
                    "Effect": "Allow",
                    "Principal": {"Service": "ec2.amazonaws.com"},
                    "Action": "sts:AssumeRole"
                }]
            }'
        aws iam attach-role-policy \
            --role-name "${ECS_INSTANCE_ROLE}" \
            --policy-arn "arn:aws:iam::aws:policy/service-role/AmazonEC2ContainerServiceforEC2Role"
        echo "Created role: ${ECS_INSTANCE_ROLE}"
    fi

    # Attach S3 access policy for the pipeline data bucket
    S3_POLICY_NAME="${PROJECT}-s3-access"
    S3_POLICY_ARN="arn:aws:iam::${ACCOUNT_ID}:policy/${S3_POLICY_NAME}"
    if aws iam get-policy --policy-arn "${S3_POLICY_ARN}" &>/dev/null; then
        echo "Policy '${S3_POLICY_NAME}' already exists"
    else
        aws iam create-policy \
            --policy-name "${S3_POLICY_NAME}" \
            --policy-document "{
                \"Version\": \"2012-10-17\",
                \"Statement\": [{
                    \"Effect\": \"Allow\",
                    \"Action\": [
                        \"s3:GetObject\",
                        \"s3:PutObject\",
                        \"s3:DeleteObject\",
                        \"s3:ListBucket\"
                    ],
                    \"Resource\": [
                        \"arn:aws:s3:::${S3_BUCKET}\",
                        \"arn:aws:s3:::${S3_BUCKET}/*\"
                    ]
                }]
            }"
        echo "Created S3 policy: ${S3_POLICY_NAME}"
    fi
    aws iam attach-role-policy \
        --role-name "${ECS_INSTANCE_ROLE}" \
        --policy-arn "${S3_POLICY_ARN}" 2>/dev/null || true

    # Create instance profile if needed
    if aws iam get-instance-profile --instance-profile-name "${ECS_INSTANCE_PROFILE}" &>/dev/null; then
        echo "Instance profile '${ECS_INSTANCE_PROFILE}' already exists"
    else
        aws iam create-instance-profile \
            --instance-profile-name "${ECS_INSTANCE_PROFILE}"
        aws iam add-role-to-instance-profile \
            --instance-profile-name "${ECS_INSTANCE_PROFILE}" \
            --role-name "${ECS_INSTANCE_ROLE}"
        echo "Created instance profile: ${ECS_INSTANCE_PROFILE}"
        echo "Waiting for instance profile propagation..."
        sleep 15
    fi

    # 2c. Spot fleet role
    SPOT_FLEET_ROLE="${PROJECT}-spot-fleet-role"
    if aws iam get-role --role-name "${SPOT_FLEET_ROLE}" &>/dev/null; then
        echo "Role '${SPOT_FLEET_ROLE}' already exists"
    else
        aws iam create-role \
            --role-name "${SPOT_FLEET_ROLE}" \
            --assume-role-policy-document '{
                "Version": "2012-10-17",
                "Statement": [{
                    "Effect": "Allow",
                    "Principal": {"Service": "spotfleet.amazonaws.com"},
                    "Action": "sts:AssumeRole"
                }]
            }'
        aws iam attach-role-policy \
            --role-name "${SPOT_FLEET_ROLE}" \
            --policy-arn "arn:aws:iam::aws:policy/service-role/AmazonEC2SpotFleetTaggingRole"
        echo "Created role: ${SPOT_FLEET_ROLE}"
    fi

    # Set ARN variables from created resource names (single code path downstream)
    BATCH_SERVICE_ROLE_ARN="arn:aws:iam::${ACCOUNT_ID}:role/${BATCH_SERVICE_ROLE}"
    ECS_INSTANCE_PROFILE_ARN="arn:aws:iam::${ACCOUNT_ID}:instance-profile/${ECS_INSTANCE_PROFILE}"
    SPOT_FLEET_ROLE_ARN="arn:aws:iam::${ACCOUNT_ID}:role/${SPOT_FLEET_ROLE}"

else
    echo "--- Skipping IAM creation (using pre-existing ARNs) ---"
    echo "  Batch service role:     ${BATCH_SERVICE_ROLE_ARN}"
    echo "  ECS instance profile:   ${ECS_INSTANCE_PROFILE_ARN}"
    echo "  Spot fleet role:        ${SPOT_FLEET_ROLE_ARN}"
fi

# ---------------------------------------------------------------------------
# 3. Launch Templates (EBS sizing)
# ---------------------------------------------------------------------------
echo ""
echo "--- Creating launch templates ---"

# CPU launch template: 500 GB gp3 (Cell Ranger needs space for FASTQs + ref)
CPU_LT_NAME="${PROJECT}-cpu-lt"
if aws ec2 describe-launch-templates --launch-template-names "${CPU_LT_NAME}" --region "${AWS_REGION}" &>/dev/null; then
    echo "Launch template '${CPU_LT_NAME}' already exists"
else
    aws ec2 create-launch-template \
        --launch-template-name "${CPU_LT_NAME}" \
        --region "${AWS_REGION}" \
        --launch-template-data '{
            "BlockDeviceMappings": [{
                "DeviceName": "/dev/xvda",
                "Ebs": {
                    "VolumeSize": 500,
                    "VolumeType": "gp3",
                    "Throughput": 250,
                    "Iops": 3000,
                    "Encrypted": true
                }
            }]
        }'
    echo "Created launch template: ${CPU_LT_NAME}"
fi

# GPU launch template: 200 GB gp3
GPU_LT_NAME="${PROJECT}-gpu-lt"
if aws ec2 describe-launch-templates --launch-template-names "${GPU_LT_NAME}" --region "${AWS_REGION}" &>/dev/null; then
    echo "Launch template '${GPU_LT_NAME}' already exists"
else
    aws ec2 create-launch-template \
        --launch-template-name "${GPU_LT_NAME}" \
        --region "${AWS_REGION}" \
        --launch-template-data '{
            "BlockDeviceMappings": [{
                "DeviceName": "/dev/xvda",
                "Ebs": {
                    "VolumeSize": 200,
                    "VolumeType": "gp3",
                    "Throughput": 250,
                    "Iops": 3000,
                    "Encrypted": true
                }
            }]
        }'
    echo "Created launch template: ${GPU_LT_NAME}"
fi

# Get launch template IDs
CPU_LT_ID=$(aws ec2 describe-launch-templates \
    --launch-template-names "${CPU_LT_NAME}" \
    --region "${AWS_REGION}" \
    --query 'LaunchTemplates[0].LaunchTemplateId' --output text)
GPU_LT_ID=$(aws ec2 describe-launch-templates \
    --launch-template-names "${GPU_LT_NAME}" \
    --region "${AWS_REGION}" \
    --query 'LaunchTemplates[0].LaunchTemplateId' --output text)

# ---------------------------------------------------------------------------
# 4. Batch Compute Environments
# ---------------------------------------------------------------------------
echo ""
echo "--- Creating Batch compute environments ---"

# CPU Compute Environment (Spot)
CPU_CE="${PROJECT}-cpu-ce"
if aws batch describe-compute-environments --compute-environments "${CPU_CE}" --region "${AWS_REGION}" \
    --query 'computeEnvironments[0].status' --output text 2>/dev/null | grep -q VALID; then
    echo "Compute environment '${CPU_CE}' already exists"
else
    aws batch create-compute-environment \
        --compute-environment-name "${CPU_CE}" \
        --region "${AWS_REGION}" \
        --type MANAGED \
        --state ENABLED \
        --service-role "${BATCH_SERVICE_ROLE_ARN}" \
        --compute-resources "{
            \"type\": \"SPOT\",
            \"allocationStrategy\": \"SPOT_CAPACITY_OPTIMIZED\",
            \"minvCpus\": 0,
            \"maxvCpus\": 512,
            \"instanceTypes\": [
                \"m6a.xlarge\",
                \"m6a.2xlarge\",
                \"r6a.2xlarge\",
                \"r6a.4xlarge\",
                \"r6a.8xlarge\"
            ],
            \"subnets\": [$(echo "${VPC_SUBNETS}" | sed 's/,/\",\"/g' | sed 's/^/\"/' | sed 's/$/\"/')],
            \"securityGroupIds\": [$(echo "${SECURITY_GROUPS}" | sed 's/,/\",\"/g' | sed 's/^/\"/' | sed 's/$/\"/')],
            \"instanceRole\": \"${ECS_INSTANCE_PROFILE_ARN}\",
            \"spotIamFleetRole\": \"${SPOT_FLEET_ROLE_ARN}\",
            \"launchTemplate\": {
                \"launchTemplateId\": \"${CPU_LT_ID}\",
                \"version\": \"\$Latest\"
            },
            \"tags\": {
                \"Project\": \"scrna-pipeline\",
                \"Environment\": \"batch\"
            }
        }"
    echo "Created compute environment: ${CPU_CE}"
fi

# GPU Compute Environment (Spot)
GPU_CE="${PROJECT}-gpu-ce"
if aws batch describe-compute-environments --compute-environments "${GPU_CE}" --region "${AWS_REGION}" \
    --query 'computeEnvironments[0].status' --output text 2>/dev/null | grep -q VALID; then
    echo "Compute environment '${GPU_CE}' already exists"
else
    aws batch create-compute-environment \
        --compute-environment-name "${GPU_CE}" \
        --region "${AWS_REGION}" \
        --type MANAGED \
        --state ENABLED \
        --service-role "${BATCH_SERVICE_ROLE_ARN}" \
        --compute-resources "{
            \"type\": \"SPOT\",
            \"allocationStrategy\": \"SPOT_CAPACITY_OPTIMIZED\",
            \"minvCpus\": 0,
            \"maxvCpus\": 128,
            \"instanceTypes\": [
                \"g4dn.xlarge\",
                \"g4dn.2xlarge\",
                \"g5.xlarge\",
                \"g5.2xlarge\"
            ],
            \"subnets\": [$(echo "${VPC_SUBNETS}" | sed 's/,/\",\"/g' | sed 's/^/\"/' | sed 's/$/\"/')],
            \"securityGroupIds\": [$(echo "${SECURITY_GROUPS}" | sed 's/,/\",\"/g' | sed 's/^/\"/' | sed 's/$/\"/')],
            \"instanceRole\": \"${ECS_INSTANCE_PROFILE_ARN}\",
            \"spotIamFleetRole\": \"${SPOT_FLEET_ROLE_ARN}\",
            \"launchTemplate\": {
                \"launchTemplateId\": \"${GPU_LT_ID}\",
                \"version\": \"\$Latest\"
            },
            \"tags\": {
                \"Project\": \"scrna-pipeline\",
                \"Environment\": \"batch\"
            }
        }"
    echo "Created compute environment: ${GPU_CE}"
fi

# Wait for compute environments to become VALID
echo ""
echo "Waiting for compute environments to become VALID..."
for ce in "${CPU_CE}" "${GPU_CE}"; do
    for i in $(seq 1 30); do
        status=$(aws batch describe-compute-environments \
            --compute-environments "${ce}" \
            --region "${AWS_REGION}" \
            --query 'computeEnvironments[0].status' --output text 2>/dev/null || echo "CREATING")
        if [[ "${status}" == "VALID" ]]; then
            echo "  ${ce}: VALID"
            break
        fi
        echo "  ${ce}: ${status} (attempt ${i}/30)..."
        sleep 10
    done
done

# ---------------------------------------------------------------------------
# 5. Job Queues
# ---------------------------------------------------------------------------
echo ""
echo "--- Creating job queues ---"

# CPU queue
CPU_QUEUE="${PROJECT}-cpu-queue"
if aws batch describe-job-queues --job-queues "${CPU_QUEUE}" --region "${AWS_REGION}" \
    --query 'jobQueues[0].status' --output text 2>/dev/null | grep -q VALID; then
    echo "Job queue '${CPU_QUEUE}' already exists"
else
    aws batch create-job-queue \
        --job-queue-name "${CPU_QUEUE}" \
        --region "${AWS_REGION}" \
        --state ENABLED \
        --priority 1 \
        --compute-environment-order "order=1,computeEnvironment=arn:aws:batch:${AWS_REGION}:${ACCOUNT_ID}:compute-environment/${CPU_CE}"
    echo "Created job queue: ${CPU_QUEUE}"
fi

# GPU queue
GPU_QUEUE="${PROJECT}-gpu-queue"
if aws batch describe-job-queues --job-queues "${GPU_QUEUE}" --region "${AWS_REGION}" \
    --query 'jobQueues[0].status' --output text 2>/dev/null | grep -q VALID; then
    echo "Job queue '${GPU_QUEUE}' already exists"
else
    aws batch create-job-queue \
        --job-queue-name "${GPU_QUEUE}" \
        --region "${AWS_REGION}" \
        --state ENABLED \
        --priority 1 \
        --compute-environment-order "order=1,computeEnvironment=arn:aws:batch:${AWS_REGION}:${ACCOUNT_ID}:compute-environment/${GPU_CE}"
    echo "Created job queue: ${GPU_QUEUE}"
fi

# ---------------------------------------------------------------------------
# 6. Summary
# ---------------------------------------------------------------------------
echo ""
echo "============================================"
echo "Setup complete!"
echo "============================================"
echo ""
echo "ECR repositories:"
echo "  ${ACCOUNT_ID}.dkr.ecr.${AWS_REGION}.amazonaws.com/${PROJECT}-cpu"
echo "  ${ACCOUNT_ID}.dkr.ecr.${AWS_REGION}.amazonaws.com/${PROJECT}-gpu"
echo ""
echo "Compute environments:"
echo "  ${CPU_CE} (CPU Spot, max 512 vCPUs)"
echo "  ${GPU_CE} (GPU Spot, max 128 vCPUs)"
echo ""
echo "Job queues:"
echo "  ${CPU_QUEUE}"
echo "  ${GPU_QUEUE}"
echo ""
echo "Next steps:"
echo "  1. Build and push Docker images:  bash aws/push_images.sh"
echo "  2. Upload data to S3:             aws s3 sync data/ s3://${S3_BUCKET}/config/"
echo "  3. Upload FASTQs:                 aws s3 sync /path/to/fastq/ s3://${S3_BUCKET}/fastq/"
echo "  4. Update nextflow.aws.config with your account ID (${ACCOUNT_ID})"
echo "  5. Run pipeline:                  nextflow run python_pipeline.nf -c nextflow.aws.config"
echo ""
echo "Tip: Set S3_BUCKET to change the data bucket (default: kstme-scrna)."
echo "     For SSO accounts, see aws/iam_admin_instructions.md."
