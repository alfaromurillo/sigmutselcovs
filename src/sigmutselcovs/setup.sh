#!/bin/bash
# Deprecated: use `sigmutselcovs download-gtex` (or, from Python,
# sigmutselcovs.download.ensure_gtex_gct). The GCT now lives in the
# user cache, never inside the installed package.
echo "setup.sh is deprecated; running sigmutselcovs download-gtex..."
exec python -m sigmutselcovs.cli download-gtex "$@"
