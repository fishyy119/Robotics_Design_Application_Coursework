# hw24 PPO

面向 Gymnasium `LunarLander-v3` 的 PPO 作业工程，使用 Python 3.12、PyTorch、`pyproject.toml` 和严格模式 Pyright。

## 项目结构

```text
hw24_PPO/
├─ pyproject.toml
├─ README.md
├─ reference/
├─ hw24_ppo/
│  ├─ __init__.py
│  ├─ __main__.py
│  ├─ buffer.py
│  ├─ config.py
│  ├─ curriculum.py
│  ├─ env.py
│  ├─ evaluate.py
│  ├─ networks.py
│  ├─ paths.py
│  ├─ plot.py
│  ├─ storage.py
│  ├─ train.py
│  ├─ types.py
│  └─ video_grid.py
└─ artifacts/
```

## 当前实现

- 默认环境为 `LunarLander-v3`。
- 仅支持一维向量观测和离散动作空间。
- 提供单一 CLI：`train`、`eval`、`plot`、`stitch-videos`。
- 训练默认启用高度里程碑奖励：当高度低于目标阈值、水平位置与倾角都落在随高度线性收紧的倒三角容差内，且横向速度低于阈值时，首次达成该里程碑会给一次进展奖励。
- 使用 `jaxtyping` 为 NumPy / Torch 张量补尺寸注解。
- checkpoint 采用 `state_dict` 字典保存，不保存整个模型对象。

## 运行方式

在 `hw24_PPO/` 目录下执行：

```bash
python -m hw24_ppo train --run-name smoke --epochs 2 --steps-per-epoch 512
python -m hw24_ppo eval --checkpoint artifacts/smoke/checkpoint_latest.pt --episodes 3
python -m hw24_ppo plot --metrics-csv artifacts/smoke/metrics.csv
python -m hw24_ppo stitch-videos --prefix artifacts/class/eval_videos/checkpoint_best2500_deterministic --grid-size 2 5
```

如果使用 `conda`，推荐环境名为 `py312`。

## 训练输出

每次训练会在 `artifacts/<run_name>/` 下生成：

- `config.json`
- `metrics.csv`
- `checkpoint_latest.pt`
- `checkpoint_best.pt`
- `eval.csv`
- `eval_videos/*.mp4`
- `training_curves.png`
- `training_curves_return_only.png`

训练时控制台会额外打印奖励分解：位置、速度、姿态、腿接触、发动机、高度里程碑奖励、终局奖励、超时惩罚

评估时会为每个 episode 额外保存一个 mp4 视频，文件名包含 checkpoint 名、策略模式（`deterministic` 或 `stochastic`）以及回合序号。

`stitch-videos` 会自动收集指定前缀的视频文件，按最长视频补齐其他视频的末帧，再按给定网格拼接为一个总视频。该命令依赖系统中的 `ffmpeg` / `ffprobe` 可执行文件。

## 说明

- `reference/` 中的脚本只作为算法思路参考，不直接作为工程实现基线。
- 当前版本重点是工程闭环、类型检查和可运行性，不保证默认超参数已经稳定达到 solved 水平。
