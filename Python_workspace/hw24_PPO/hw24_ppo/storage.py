import csv
import json
from pathlib import Path
from typing import Mapping, Sequence, cast

import torch


def write_json(path: Path, payload: Mapping[str, object]) -> None:
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def initialize_csv(path: Path, fieldnames: Sequence[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()


def append_csv_row(
    path: Path,
    fieldnames: Sequence[str],
    row: Mapping[str, object],
) -> None:
    with path.open("a", encoding="utf-8", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writerow(row)


def save_checkpoint(path: Path, checkpoint: Mapping[str, object]) -> None:
    torch.save(dict(checkpoint), path)


def load_checkpoint(path: Path, device: torch.device) -> dict[str, object]:
    if not path.exists():
        raise FileNotFoundError(f"Checkpoint not found: {path}")

    try:
        checkpoint = torch.load(path, map_location=device, weights_only=False)
    except TypeError:
        checkpoint = torch.load(path, map_location=device)

    if not isinstance(checkpoint, dict):
        raise ValueError("Checkpoint payload is not a dictionary.")

    return cast(dict[str, object], checkpoint)
