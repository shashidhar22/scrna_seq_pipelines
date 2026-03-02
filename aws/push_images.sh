#!/usr/bin/env bash
# =============================================================================
# Build and Push Docker Images to ECR — scRNA-seq Pipeline
# =============================================================================
#
# Builds CPU and GPU Docker images and pushes them to ECR.
#
# Prerequisites:
#   - AWS CLI v2 configured
#   - Docker installed and running
#   - Cell Ranger tarball placed in aws/ directory (for CPU image)
#   - Reference genomes in aws/ directory (for CPU image):
#       aws/refdata-gex-GRCh38-2024-A/
#       aws/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0/
#
# Usage:
#   bash aws/push_images.sh            # Build and push both images
#   bash aws/push_images.sh --cpu      # CPU image only
#   bash aws/push_images.sh --gpu      # GPU image only
# =============================================================================

set -euo pipefail

AWS_REGION="${AWS_REGION:-us-east-1}"
ACCOUNT_ID=$(aws sts get-caller-identity --query Account --output text)
PROJECT="scrna"
REPO_BASE="${ACCOUNT_ID}.dkr.ecr.${AWS_REGION}.amazonaws.com"

# Navigate to project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "${SCRIPT_DIR}")"
cd "${PROJECT_ROOT}"

echo "============================================"
echo "Docker Image Build & Push"
echo "============================================"
echo "Account:  ${ACCOUNT_ID}"
echo "Region:   ${AWS_REGION}"
echo "Registry: ${REPO_BASE}"
echo "============================================"

# Authenticate Docker with ECR
echo ""
echo "Authenticating Docker with ECR..."
aws ecr get-login-password --region "${AWS_REGION}" \
    | docker login --username AWS --password-stdin "${REPO_BASE}"

BUILD_CPU=true
BUILD_GPU=true
if [[ "${1:-}" == "--cpu" ]]; then
    BUILD_GPU=false
elif [[ "${1:-}" == "--gpu" ]]; then
    BUILD_CPU=false
fi

# ---------------------------------------------------------------------------
# CPU Image
# ---------------------------------------------------------------------------
if [[ "${BUILD_CPU}" == "true" ]]; then
    echo ""
    echo "--- Building CPU image ---"
    echo ""

    # Verify required build context files
    if [[ ! -f "aws/cellranger-8.0.1.tar.gz" ]]; then
        echo "ERROR: Cell Ranger tarball not found at aws/cellranger-8.0.1.tar.gz"
        echo "Download from: https://www.10xgenomics.com/support/software/cell-ranger/downloads"
        exit 1
    fi
    if [[ ! -d "aws/refdata-gex-GRCh38-2024-A" ]]; then
        echo "ERROR: GEX reference not found at aws/refdata-gex-GRCh38-2024-A/"
        echo "Download from: https://www.10xgenomics.com/support/software/cell-ranger/downloads"
        exit 1
    fi
    if [[ ! -d "aws/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0" ]]; then
        echo "ERROR: VDJ reference not found at aws/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0/"
        echo "Download from: https://www.10xgenomics.com/support/software/cell-ranger/downloads"
        exit 1
    fi

    # Copy shared files into build context
    cp requirements_cpu.txt aws/requirements_cpu.txt
    cp -r pipeline/ aws/pipeline/

    docker build \
        -f aws/Dockerfile.cpu \
        -t "${PROJECT}-cpu:latest" \
        aws/

    # Clean up build context copies
    rm -f aws/requirements_cpu.txt
    rm -rf aws/pipeline/

    # Tag and push
    docker tag "${PROJECT}-cpu:latest" "${REPO_BASE}/${PROJECT}-cpu:latest"
    docker push "${REPO_BASE}/${PROJECT}-cpu:latest"
    echo "Pushed: ${REPO_BASE}/${PROJECT}-cpu:latest"
fi

# ---------------------------------------------------------------------------
# GPU Image
# ---------------------------------------------------------------------------
if [[ "${BUILD_GPU}" == "true" ]]; then
    echo ""
    echo "--- Building GPU image ---"
    echo ""

    # Copy shared files into build context
    cp requirements.txt aws/requirements.txt
    cp -r pipeline/ aws/pipeline/

    docker build \
        -f aws/Dockerfile.gpu \
        -t "${PROJECT}-gpu:latest" \
        aws/

    # Clean up build context copies
    rm -f aws/requirements.txt
    rm -rf aws/pipeline/

    # Tag and push
    docker tag "${PROJECT}-gpu:latest" "${REPO_BASE}/${PROJECT}-gpu:latest"
    docker push "${REPO_BASE}/${PROJECT}-gpu:latest"
    echo "Pushed: ${REPO_BASE}/${PROJECT}-gpu:latest"
fi

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo ""
echo "============================================"
echo "Done!"
echo "============================================"
echo ""
if [[ "${BUILD_CPU}" == "true" ]]; then
    echo "CPU: ${REPO_BASE}/${PROJECT}-cpu:latest"
fi
if [[ "${BUILD_GPU}" == "true" ]]; then
    echo "GPU: ${REPO_BASE}/${PROJECT}-gpu:latest"
fi
echo ""
echo "Update nextflow.aws.config with these image URIs."
