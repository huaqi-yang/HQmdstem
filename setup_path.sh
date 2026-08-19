#!/usr/bin/env bash
set -u
HOME_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$HOME/.local/bin"
ln -sf "$HOME_DIR/HQmdstemkit.sh" "$HOME/.local/bin/HQmdstemkit.sh"
chmod +x "$HOME_DIR/HQmdstemkit.sh" "$HOME/.local/bin/HQmdstemkit.sh"
if ! grep -q '.local/bin' "$HOME/.bashrc" 2>/dev/null; then
  echo 'export PATH="$HOME/.local/bin:$PATH"' >> "$HOME/.bashrc"
fi
export PATH="$HOME/.local/bin:$PATH"
echo "HQmdstemkit.sh -> $(command -v HQmdstemkit.sh)"
HQmdstemkit.sh help | head -30