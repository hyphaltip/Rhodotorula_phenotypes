#!/bin/bash
# plink_arch.sh — resolve the correct plink2 binary for the CURRENT node.
# UCR HPCC ships three AVX2 variants under /opt/linux/centos/8.x/x86_64/pkgs/plink2/...:
#   amd_avx2 (AMD EPYC), intel_avx2 (Intel), x86_64 (portable).
# Usage:  source plink_arch.sh   ->  sets $PLINK to a valid binary (echoes it).
# Fallback order: vendor-matched avx2 > portable x86_64 > any other existing binary.
PLINK_DIR=/opt/linux/centos/8.x/x86_64/pkgs/plink2/2.0.0-a.7.3/bin
if [ ! -d "$PLINK_DIR" ]; then
  PLINK_DIR=/opt/linux/rocky/8.x/x86_64/pkgs/plink2/2.0.0-a.7.3/bin
fi

VENDOR=$(grep -m1 "^vendor_id" /proc/cpuinfo | awk '{print $3}')
case "$VENDOR" in
  GenuineIntel) ARCH=intel_avx2 ;;
  AuthenticAMD) ARCH=amd_avx2 ;;
  *)            ARCH= ;;
esac
PLINK=""
for cand in "${ARCH:+$PLINK_DIR/$ARCH/plink2}" "$PLINK_DIR/x86_64/plink2"; do
  if [ -n "$cand" ] && [ -x "$cand" ]; then PLINK="$cand"; break; fi
done
# last resort: any plink2 binary
if [ -z "$PLINK" ] && [ -d "$PLINK_DIR" ]; then
  PLINK=$(find "$PLINK_DIR" -maxdepth 2 -name plink2 -type f | head -1)
fi
if [ -z "$PLINK" ] || [ ! -x "$PLINK" ]; then
  echo "ERROR: no plink2 binary found under $PLINK_DIR" >&2
  return 1 2>/dev/null || exit 1
fi
echo "PLINK=$PLINK (vendor=${VENDOR:-unknown})" >&2
export PLINK
