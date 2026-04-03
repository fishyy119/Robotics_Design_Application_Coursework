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
- `AGENTS.md` is the authoritative machine-readable style guide for this repository.

## MATLAB

- `reference/` stores teacher-provided example code. Use it as content or algorithm reference only, not as a code style reference.
- Third-party or auto-generated files may be kept for functionality, but do not use them as the primary style baseline for new edits.
- `MATLAB_workspace/+utils/` is the shared location for reusable plotting utilities. When creating plots, prefer existing functions there for figure initialization and common graphics settings.
- Use clear script sectioning. Separate parameter setup, data computation, and plotting with concise `%%` blocks such as `%% 参数`, `%% 初始条件`, `%% 数值积分`, `%% 后处理`, and `%% 绘图`.
- Avoid decorative separators such as repeated `=` or `-` in comments and section headers.
- Do not couple plotting logic with data computation.
- Write MATLAB code comments in Chinese.
- Keep identifiers in English or conventional symbolic form such as `q1`, `u2`, `KE`, and `params`.
- Keep all figure text in English, including titles, axis labels, legends, and annotations.
- Keep comments concise and informative. Explain physical meaning, assumptions, units, or non-obvious implementation choices; do not keep stale example comments or commented-out scratch code.
- For reusable functions, add a short Chinese summary comment immediately below the `function` line, plus concise input/output notes only when they are helpful.
- Use consistent whitespace around operators, commas, and matrix delimiters to keep formulas readable.
- When function signatures become too long or parameter passing becomes cumbersome, prefer `classdef` for parameter encapsulation.

## LaTeX

- The LaTeX preamble should strictly follow existing documents in the workspace.
- The build target is a single-page document.
- Adjusting page dimensions is acceptable when needed to achieve the intended single-page layout.
- Keep new homework documents aligned with the existing `hwNN.tex` organization and single-document entry style under `LaTeX_workspace/`.
- Treat `LaTeX_workspace/build/` as build output, not source content.
- Write LaTeX source comments in Chinese, using plain `% 注释` style comments. Avoid decorative dashed separators in TikZ or equation blocks unless they already convey structure that would otherwise be unclear.

## Validation

- Do not automatically execute MATLAB scripts in this repository. If MATLAB execution is needed, ask the user to run it manually.
- For LaTeX changes, preserve successful single-page compilation as the primary output goal.
