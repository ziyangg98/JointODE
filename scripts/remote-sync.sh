#!/bin/bash
# Sync local code to server and run command
# Usage:
#   ./scripts/remote-sync.sh              # sync only
#   ./scripts/remote-sync.sh install      # sync + install package
#   ./scripts/remote-sync.sh test         # sync + run tests (no reinstall)
#   ./scripts/remote-sync.sh check        # sync + R CMD check
#   ./scripts/remote-sync.sh restore      # sync + renv::restore()
#   ./scripts/remote-sync.sh run FILE.R   # sync + run arbitrary R script

SERVER="sky18"
REMOTE_DIR="~/JointODE"

# --- Sync ---
rsync -az --delete \
  --exclude='.git' \
  --exclude='src/*.o' \
  --exclude='src/*.so' \
  --exclude='src/*.dll' \
  --exclude='Meta' \
  --exclude='*.Rcheck' \
  --exclude='renv/library' \
  --exclude='renv/staging' \
  --exclude='renv/sandbox' \
  ./ "${SERVER}:${REMOTE_DIR}/"

echo "Synced to ${SERVER}:${REMOTE_DIR}"
[ -z "${1}" ] && exit 0

# --- Helper: run R command with renv library ---
remote_r() {
  ssh "${SERVER}" "cd ${REMOTE_DIR} && \
    export R_LIBS=\$(Rscript --vanilla -e 'source(\"renv/activate.R\"); cat(paste(.libPaths(), collapse=\":\"))' 2>/dev/null) && \
    $*"
}

case "${1}" in
  restore)
    ssh "${SERVER}" "cd ${REMOTE_DIR} && \
      Rscript --vanilla -e 'source(\"renv/activate.R\"); renv::restore(prompt = FALSE)'"
    ;;
  install)
    remote_r R CMD INSTALL --no-docs --no-multiarch .
    ;;
  test)
    remote_r R CMD INSTALL --no-docs --no-multiarch . '&&' \
      Rscript -e "'library(JointODE); testthat::test_local()'"
    ;;
  check)
    remote_r R CMD build --no-manual . '&&' \
      R CMD check --no-manual JointODE_*.tar.gz
    ;;
  run)
    [ -z "${2}" ] && echo "Usage: $0 run FILE.R" && exit 1
    remote_r Rscript "${2}"
    ;;
  *)
    echo "Unknown: ${1}. Use: install|test|check|bench|restore|run"
    ;;
esac
