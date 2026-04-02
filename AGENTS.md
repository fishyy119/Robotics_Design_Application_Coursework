# AGENTS

## Scope

- This repository is a robotics coursework project.
- `MATLAB_workspace/` contains MATLAB simulation source code.
- `LaTeX_workspace/` contains theory and derivation documents.
- Coursework content is organized with `hwNN_*` or `hwNN` prefixes. Preserve this classification unless the user explicitly requests restructuring.

## General

- Do not introduce conventions that are not stated here or already established in the repository.
- When a style or structure decision is unclear, follow nearby existing project files.
- Keep edits local and minimal. Do not perform broad rewrites unless requested.

## MATLAB

- `reference/` stores teacher-provided example code.
- Files under `reference/` may be used as content or algorithm reference, but not as code style reference.
- `MATLAB_workspace/+utils/` is the shared location for reusable plotting utilities.
- When creating plots, prefer existing functions in `+utils/` for figure initialization and common graphics settings.
- Use clear script sectioning. Separate parameter setup, data computation, and plotting with `%%` blocks where appropriate.
- Do not couple plotting logic with data computation.
- Write MATLAB code comments in Chinese.
- Keep all figure text in English, including titles, axis labels, legends, and annotations.
- When function signatures become too long or parameter passing becomes cumbersome, prefer `classdef` for parameter encapsulation.

## LaTeX

- The LaTeX preamble should strictly follow existing documents in the workspace.
- The build target is a single-page document.
- Adjusting page dimensions is acceptable when needed to achieve the intended single-page layout.
- Keep new homework documents aligned with the existing `hwNN.tex` organization.
- Preserve the existing single-document entry style under `LaTeX_workspace/`.
- Treat `LaTeX_workspace/build/` as build output, not source content.

## Validation

- Do not automatically execute MATLAB scripts in this repository.
- If MATLAB execution is needed, ask the user to run it manually.
- For LaTeX changes, preserve successful single-page compilation as the primary output goal.
