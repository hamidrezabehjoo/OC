# Panel registration record (unpinned on purpose)

This file is **not** covered by `panel/MANIFEST.md5`. It records the
external timestamp of the frozen protocol without invalidating the
manifest (filling it after archival must not change any pinned file).

Fill in before running `./panel/run_panel.sh screen`:

- git tag: `panel-v2.5`
- commit hash: `2c64ee3` (sandbox repo /mnt/agents/output/paper;
  re-derive after moving to the real repo: commit same content, verify
  manifest, then tag — hash there will differ and takes precedence)
- md5 of `panel/MANIFEST.md5`: `74944cebc45f2c424e79af8384bb0eec`
- repository URL: ____________________
- GitHub release: ____________________
- Zenodo DOI (or OSF registration): ____________________
- date archived (UTC): ____________________

Procedure:
1. Move this tree into the real project repo (same content).
2. `python3 panel/make_manifest.py --verify` — must pass.
3. Tag `panel-v2.5`, push, create a GitHub release, let Zenodo mint the
   DOI (or create an OSF registration).
4. Fill this file in. Do not re-tag afterwards.
