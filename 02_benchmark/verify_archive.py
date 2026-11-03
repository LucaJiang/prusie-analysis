#!/usr/bin/env python3
"""Authenticate exact execution bytes, including the documented restored copies."""
from pathlib import Path
from common import read, sha


def main():
    archive = Path(__file__).resolve().parent / 'archive'
    restored = []
    total = 0
    for manifest in sorted(archive.rglob('checksums.json')):
        checksums = read(manifest)
        for relative, expected in checksums.items():
            total += 1
            published = manifest.parent / relative
            if sha(published) == expected:
                continue
            original = manifest.parent / 'as_run' / relative
            if not original.is_file() or sha(original) != expected:
                raise ValueError(f'Missing or changed execution bytes: {relative}')
            restored.append(relative)
    print(f'Authenticated {total} original archive files; {len(restored)} use restored copies.')


if __name__ == '__main__':
    main()
