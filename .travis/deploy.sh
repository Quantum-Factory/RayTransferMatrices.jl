#!/bin/sh
mv docs/build docs/RayTransferMatrices.jl
echo 'put -r docs/RayTransferMatrices.jl' |sftp -P $SSH_PORT $SSH_USER@$SSH_HOST:docs/github.com/Quantum-Factory
