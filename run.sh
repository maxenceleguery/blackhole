#!/usr/bin/env bash
# Launch on the discrete Nvidia GPU. NixOS hybrid laptop, PRIME offload, Wayland session.
#
# Root cause chain (diagnosed step by step):
#  1. browsers render via EGL; libglvnd EGL enumerates both GPUs and picks Intel.
#     -> force the Nvidia-only EGL vendor json below.
#  2. on a Wayland session Chromium picks wayland-ozone + Vulkan, which crash-loops
#     on Nvidia. -> force x11 ozone + GL ANGLE backend (no Vulkan). Firefox: XWayland.
here="$(cd "$(dirname "$0")" && pwd)"
nv_egl=$(ls /run/opengl-driver/share/glvnd/egl_vendor.d/*nvidia*.json 2>/dev/null | head -1)
[ -z "$nv_egl" ] && { echo "no nvidia EGL vendor json under /run/opengl-driver"; exit 1; }
export __NV_PRIME_RENDER_OFFLOAD=1 __GLX_VENDOR_LIBRARY_NAME=nvidia
export __EGL_VENDOR_LIBRARY_FILENAMES="$nv_egl"
export MOZ_ENABLE_WAYLAND=0                       # Firefox -> XWayland

B="${BROWSER:-firefox}"
case "$B" in
  *chrom*) exec "$B" --ozone-platform=x11 --use-angle=gl --ignore-gpu-blocklist "$here/index.html" ;;
  *)       exec "$B" "$here/index.html" ;;
esac
