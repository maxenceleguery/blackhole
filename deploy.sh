#!/usr/bin/env bash
# Push the sim + SEO files to the VPS docroot. Run after the vhost + cert exist.
set -e
here="$(cd "$(dirname "$0")" && pwd)"
files=(index.html qnm.html gw.html binary.html robots.txt sitemap.xml)
[ -f "$here/preview.png" ] && files+=(preview.png)   # 1200x630 OG/social card, if you add one
scp "${files[@]/#/$here/}" ovh:/var/www/blackhole/
echo "deployed -> https://blackhole.maxenceleguery.net"
