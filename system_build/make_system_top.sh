#!/bin/bash
set -euo pipefail

HEADER="system_header.top"
INSANE="insane.top"
OUT="system.top"

USER_NAME="${1:-}"   # Optional argument

# Extract only the include lines from system_header.top
grep '^#include' "$HEADER" > "$OUT"

# Write the [ system ] block
echo "" >> "$OUT"
echo "[ system ]" >> "$OUT"
echo "; name" >> "$OUT"

if [[ -n "$USER_NAME" ]]; then
    # User provided a custom name
    echo "$USER_NAME" >> "$OUT"
else
    # Extract the name from insane.top (line after '; name')
    awk '
        found_name && $0 !~ /^;/ { print; exit }
        /^; name/ { found_name=1 }
    ' "$INSANE" >> "$OUT"
fi

echo "" >> "$OUT"

# Append everything from insane.top starting at [ molecules ]
awk '
    /^\[ molecules \]/ { printing=1 }
    printing
' "$INSANE" >> "$OUT"

echo "Wrote merged topology to $OUT"
